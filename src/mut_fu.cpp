

#include <string.h>
#include <time.h>

#include "cmain.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
	extern FILE *tsld;
#endif

extern PROGARGS ArgKit;

extern vector < XROSOMA > vecDNK;
vector < SAMPLE >  vecSAMPL;

//const char mutaFiles[NO_CANCER_ID][32] = {
//    "snv_mnv_BLCA-US.tsv", "snv_mnv_BRCA-US.tsv", "snv_mnv_CESC-US.tsv",
//    "snv_mnv_HNSC-US.tsv", "snv_mnv_LUAD-US.tsv", "snv_mnv_LUSC-US.tsv"
//};

void parsMutRec( char *pB, int *xrNum, int *Xpos, char *chREF, char *chALT,
                string &samp );

////////////////////////////////////////////////////////////////////////

bool lesser_MUT ( const MUTANT &x1, const MUTANT &x2 )
{
    return x1.nucPos < x2.nucPos;
}
/////////////////////////////////////////////////////////////////////////////////

int  xtrSamplName ( const char  *vcf_Fname, char *SamplName)
{

    int n=0;
    const char *p;
    for ( p=vcf_Fname; *p;  p++, n++) {
        if ( *p == '.' )
            break;
        *(SamplName+n) = *p;
    }
    if ( *p != '.' )   {
        printf ("\nCann't extract SAMPLE_NAME from '%s'\n", vcf_Fname);
        return 0;
    }
    *(SamplName+n) = '\0';
    
    return 1;
}

////////////////////////////////////////////////////////////////////////
/*
int  xtrSamplName ( const char  *vcf_Fpath, char *SamplName, char *Folder)
{
    const char *pFp;
    
    pFp = vcf_Fpath + strlen(vcf_Fpath) - 1;
    
    while ( *pFp != '/' ) {
        if ( pFp == vcf_Fpath ) {
            printf ( "INV.MutationFilePath (2-nd parm) : '%s'\n", vcf_Fpath );
            return 0;
        }
        pFp--;
    }
    pFp++;  //  '/'
    
    const char *p;
    int n=0;
    if ( Folder )   {
        for ( p=vcf_Fpath; p < pFp;  p++, n++)
            *(Folder+n) = *p;
        *(Folder+n) = '\0';
    }
    
    n=0;
    for ( p=pFp; *p;  p++, n++) {
        if ( *p == '.' )
            break;
        *(SamplName+n) = *p;
    }
    if ( *p != '.' )   {
        printf ("\nCann't extract SAMPLE_NAME from '%s'\n", vcf_Fpath);
        return 0;
    }
    *(SamplName+n) = '\0';
    
    return 1;
}
*/
/////////////////////////////////////////////////////////////////////////////////

int loadVCFMutation( const char *vcf_Fname )
{
    char cBuffer[4096];
    int XroID, nXro, NucPos;
    char cREF, cALT, cAPO;
    char cXr[16];
    int cntAllMut=0;
    
    printf("\nLoadMutations from '%s' ......\n", vcf_Fname);
    
    sprintf(cBuffer, "%s%s", ArgKit.MUTdir.c_str(), vcf_Fname );
    FILE *f_Mutas=fopen( cBuffer, "rb");
    if ( f_Mutas==NULL ) {
        printf ("\nNOT_FOUND Mutations_File ='%s'\n", cBuffer );
        return -1;
    }

    if ( ! xtrSamplName ( vcf_Fname, cBuffer) ) {
        fclose(f_Mutas);
        return -1;
    }
    
    if  ( vecSAMPL.empty() )
        vecSAMPL.push_back(SAMPLE(cBuffer));
    else {
        vecSAMPL[0].clearSampl( );
        vecSAMPL[0].SampName = cBuffer;
    }
    
    while ( 1 ) {
        memset( (void *)cBuffer, '\0', sizeof(cBuffer));
        if ( fgets_ShortRec(cBuffer, sizeof(cBuffer)-1, f_Mutas )==0 )
            break;
        if ( cBuffer[0] == '#' )
            continue;
        
        NucPos = -1; cREF = '?'; cALT = '?';
        //                  Pos  Sam  Ref Alt Qu  Filt  Info Clu  C_typ Xr
        if ( sscanf(cBuffer, "%s\t%d\t%*s\t%c\t%c\t%*s\n", cXr, &NucPos, &cREF, &cALT) != 4 ) {
            printf("File '%s': ERR_Rec_FRMT:: cntMut=%d : [X=%s, Pos=%d, Ref=%c, Alt=%c\n\t'%s'\n",
                   vcf_Fname, cntAllMut,  cXr, NucPos, cREF, cALT, cBuffer);
            fclose(f_Mutas);
            return -1;
        }
        if ( cBuffer[0] == 'X' )
            XroID = 23;
        else
            if ( cBuffer[0] == 'Y' )
                XroID = 24;
            else
                sscanf(cBuffer, "%d", &XroID);
        nXro = findXroByID(XroID);
        if ( nXro  < 0 ) {
            printf("\tINV.chrID=%d  cntMut=%d : '%s'\n", XroID, cntAllMut, cBuffer);
            continue;
        }
        vecDNK[nXro].testValidDNK(NucPos, cREF, cALT );
        cAPO = ( vecDNK[nXro].APOtest(NucPos, cREF, cALT) > 0 ) ? 1 : 0;
        if ( vecSAMPL[0].vMutAPO[nXro].empty() || vecSAMPL[0].vMutAPO[nXro].back().nucPos <= NucPos )
            vecSAMPL[0].vMutAPO[nXro].push_back( MUTANT(NucPos, cREF, cALT, cAPO ));   //chAPO
        else {
            vector <MUTANT>::iterator Iter =
            lower_bound( vecSAMPL[0].vMutAPO[nXro].begin( ), vecSAMPL[0].vMutAPO[nXro].end( ),
                        MUTANT(NucPos),  lesser_MUT);
            vecSAMPL[0].vMutAPO[nXro].insert( Iter, MUTANT(NucPos, cREF, cALT, cAPO ) );   //chAPO
        }
        cntAllMut += 1;
    }
    
    fclose(f_Mutas);
    
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ )
        if ( vecSAMPL[0].maxMUTsize < vecSAMPL[0].vMutAPO[nX].size() )
            vecSAMPL[0].maxMUTsize = vecSAMPL[0].vMutAPO[nX].size() + 1;
    
    return cntAllMut;
}
/////////////////////////////////////////////////////////////////////////////////

void SAMPLE::clearSampl( )
{
    int n;
    for ( n=0; n<NO_HUMAN_XRO; n++ )    {
        vMutAPO[n].clear();
        PorogMut[n] = 0;   //ClustPorog
    
        vClust[n].clear();
        vAggrC[n].clear();
        PorogClust[n] = 0;
    }
}
/////////////////////////////////////////////////////////////////////////////////

void tst_loadMut()
{
    FILE *Ftest;
    char buffr[1024];
    
    sprintf(buffr,"%stestLoadMut.txt", ArgKit.OUTdir.c_str() );
    Ftest =fopen(buffr, "w");
    fprintf(Ftest, "nS\tXRO\tPOS\tREF\tALT\tSampl\n");
    for ( int nS=0; nS<vecSAMPL.size(); nS++) {
        for (int nX=0; nX<NO_HUMAN_XRO; nX++)
            for (int nM=0; nM<vecSAMPL[nS].vMutAPO[nX].size(); nM++)
                fprintf(Ftest, "%d\t%d\t%ld\t%c\t%c\t%s\n", nS, nX,
                        vecSAMPL[nS].vMutAPO[nX][nM].nucPos,
                        vecSAMPL[nS].vMutAPO[nX][nM].nucREF,
                        vecSAMPL[nS].vMutAPO[nX][nM].nucALT,
                        vecSAMPL[nS].SampName.c_str() );
    }
    fclose(Ftest);
    return;
}
///////////////////////////////////////////////////////

void print_MutSize()
{
    FILE *Ftest;
    char buffr[1024];
    int nX, nS;
    
    sprintf(buffr,"%sMutSizes.txt", ArgKit.OUTdir.c_str() );
    Ftest =fopen(buffr, "w");
    
    fprintf(Ftest, "XRO\tXsize");
    for ( nS=0; nS<vecSAMPL.size(); nS++)
        fprintf(Ftest, "\tS.%d", nS);
    fprintf(Ftest, "\n");
    
    for ( nX=0; nX<NO_HUMAN_XRO; nX++)  {
        fprintf(Ftest, "%d\t%ld", vecDNK[nX].chrNum, vecDNK[nX].Xsize );
        for ( nS=0; nS<vecSAMPL.size(); nS++)
            fprintf(Ftest, "\t%ld", vecSAMPL[nS].vMutAPO[nX].size() );
        fprintf(Ftest, "\n");
    }
    fclose(Ftest);
    return;
}
///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

void parsMutRec( char *pB, int *xrNum, int *Xpos, char *chREF, char *chALT, string &samp )
{

	if ( *pB == 'X' || *pB=='x' )	
		*xrNum = 23;
	else
	if ( *pB == 'Y' || *pB=='y' )
		*xrNum = 24;
	else
		*xrNum = atoi( pB );

	while ( *pB && *pB!='\t' ) pB++;	//skip CHROM
	if (*pB)	pB++;

	*Xpos = atoi( pB );								//POS
	while ( *pB && *pB!='\t' ) pB++;	//skip POS
	if (*pB)	pB++;

	while ( *pB && *pB!='\t' ) pB++;	//skip ID
	if (*pB)	pB++;

	*chREF = *pB;											//REF
	while ( *pB && *pB!='\t' ) pB++;	//skip REF
	if (*pB)	pB++;

	*chALT = *pB;											//ALT
	while ( *pB && *pB!='\t' ) pB++;	//skip ALT
	if (*pB)	pB++;

	while ( *pB && *pB!='\t' ) pB++;	//skip QUAL
	if (*pB)	pB++;

	while ( *pB && *pB!='\t' ) pB++;	//skip FILTER
	if (*pB)	pB++;
	
	while ( *pB && *pB!='\t' ) pB++;	//skip INFO
	if (*pB)	pB++;

	char *pBend = pB;
	while ( *pBend && *pBend!='\t' ) pBend ++;
	*pBend = '\0';
	samp = pB;

	return;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
