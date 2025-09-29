

#include <string.h>
#include <time.h>

#include "cmain.h"
#include "xrosoma.h"

extern PROGARGS ArgKit;

extern vector < XROSOMA > vecDNK;
string SamplName;   //vector < SAMPLE >  vecSAMPL;


#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

int parsFORMATrec( int recNum, char *pB, int *ind_Xro, MUTANT &rMUT );
int parsINFOrec( int recNum, char *pB, int *ind_Xro, MUTANT &rMUT );
//void parsMutRec( char *pB, int *xrNum, int *Xpos, char *chREF, char *chALT,
//                string &samp );

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
    int vcfFORMAT=0;        // if 0: get VAF from INFO (t_alt_count, t_ref_count)
                            // if 1: get VAF from FORMAT (AD)
    MUTANT rdMUT, emptyMUT;
    char cBuffer[4096];
    int ind_Xro=-1, curnt_X=-1;
//    int nMutXro=0;
//    clock_t startP=0, startX=0, stopP=0, stopX=0;
//    double duration;

    long startP = clock();
    
    printf("\nLoading mutations from '%s' ......\n", vcf_Fname);
    sprintf(cBuffer, "%s%s", ArgKit.MUTdir.c_str(), vcf_Fname );
    FILE *f_Mutas=fopen( cBuffer, "rb");
    if ( f_Mutas==NULL ) {
        printf ("\nFailed to open VCF file: '%s'\n", cBuffer );
        return -1;
    }

    if ( ! xtrSamplName ( vcf_Fname, cBuffer) ) {
        fclose(f_Mutas);
        return -1;
    }
    
    if  ( ! SamplName.empty() )
        clearSampl( );
        SamplName = cBuffer;
//    }
    
    int cntAllRec=0;
    int cntIgnor=0;
    while ( 1 ) {
        memset( (void *)cBuffer, '\0', sizeof(cBuffer));
        if ( fgets_ShortRec(cBuffer, sizeof(cBuffer)-1, f_Mutas )==0 )
            break;
        if ( cBuffer[0] == '#' ) {
            if ( (cBuffer[1] != '#') && strstr(cBuffer, "FORMAT" ) )
                vcfFORMAT = 1;
            continue;
        }
        rdMUT = emptyMUT;
        cntAllRec += 1;

        int retC = ( vcfFORMAT )   ? parsFORMATrec( cntAllRec,cBuffer, &ind_Xro, rdMUT )
                                   : parsINFOrec( cntAllRec, cBuffer, &ind_Xro, rdMUT );
        if ( retC < 0 ) {
            fclose(f_Mutas);
            return -1;
        }
        if ( retC == 0 )    {
            cntIgnor++;
            continue;
        }
        
//        vecDNK[ind_Xro].testValidDNK(rdMUT.nucPos, rdMUT.nucREF );   do it at parser
//          cAPO = ( vecDNK[ind_Xro].APOtest(rdMUT.nucPos, rdMUT.nucREF, crdMUT.nucALT) > 0 ) ? 1 : 0;
        if ( vecDNK[ind_Xro].vMutAPO.empty() || vecDNK[ind_Xro].vMutAPO.back().nucPos <= rdMUT.nucPos )
            vecDNK[ind_Xro].vMutAPO.push_back( rdMUT );
        else {
            vector <MUTANT>::iterator Iter =
            lower_bound( vecDNK[ind_Xro].vMutAPO.begin( ), vecDNK[ind_Xro].vMutAPO.end( ),
                        rdMUT,  lesser_MUT);
            if ( Iter->nucPos == rdMUT.nucPos ) {
                printf("\tLine %d IGNORED: DublPos of Mutation: %s %d REF'%c' ALT'%s'\n", cntAllRec,
                       vecDNK[ind_Xro].XroID.c_str(), rdMUT.nucPos, rdMUT.nucREF, rdMUT.nucALT);
                cntIgnor++;
            }
            else
                vecDNK[ind_Xro].vMutAPO.insert( Iter, rdMUT );
        }

    }
    
    fclose(f_Mutas);
    
    for ( int nX=0; nX<vecDNK.size(); nX++ )
        if ( vecDNK[nX].maxMUTsize < vecDNK[nX].vMutAPO.size() )
            vecDNK[nX].maxMUTsize = vecDNK[nX].vMutAPO.size() + 1;
    
    long stopP = clock();
    double duration = (double)(stopP - startP) / CLOCKS_PER_SEC;
    printf("=======endLoading: MUTcnt=%d Ignored=%d -> %d  dT=%5.2f\n\n", cntAllRec, cntIgnor,
           cntAllRec-cntIgnor, duration );
    
    return cntAllRec;
}
/////////////////////////////////////////////////////////////////////////////////

void  clearSampl( )
{
    int n;
    for ( n=0; n<vecDNK.size(); n++ )    {
        if ( vecDNK[n].vMutAPO.empty() )
            continue;
        vecDNK[n].vMutAPO.clear();
        vecDNK[n].vClust.clear();
        vecDNK[n].vAggrC.clear();
        vecDNK[n].PorogMut = 0;
        vecDNK[n].PorogClust = 0;
    }
}
/////////////////////////////////////////////////////////////////////////////////

/*
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
*/
/////////////////////////////////////////////////////////////////////////////////

int parsFORMATrec( int recNum, char *pRec, int *ind_Xro, MUTANT &rMUT )
{
    char cREF[1024]="?",  cALT[1024]="?", cXr[XRO_ID_SIZE]="?";
    char cFORMAT[2048]="?", N845_2[2048]="?";
    
//                       #CHROM POS ID  REF ALT  QUAL FILTER INFO FORMAT N845-2
    int nFld = sscanf(pRec, "%s\t%d\t%*s\t%s\t%s\t%*s\t%*s\t%*s\t%s\t%s\n", cXr, &rMUT.nucPos, cREF, cALT, cFORMAT, N845_2);
    if ( nFld != 6 ) {
        printf("Line %d:: Failed to parse: Unknown format: [X=%s, Pos=%d, Ref=%s, Alt=%s\n\t'%s'\n",
                   recNum+1,  cXr, rMUT.nucPos, cREF, cALT, pRec);
        return -1;
    }

    *ind_Xro = findXroByID(cXr);
    if ( *ind_Xro  < 0 ) {
        printf("\tLine %d:: Unknown chromo_name %s\n\t'%s'\n", recNum+1, cXr, pRec);
        return -1;
    }
    
    int cntNuc=0;
    char *pA=cALT;
    while ( *pA )   {
        if ( *pA==',' )
            cntNuc--;
        else
            cntNuc++;
        pA++;
    }
    if ( cntNuc > 1 || strlen(cREF) > 1  ) {          //    || strlen(cALT) > 1 ) {
//        printf("\tLine %d IGNORED: Insertion: %s %d REF'%s' ALT'%s'\n", recNum+1, cXr, rMUT.nucPos, cREF, cALT);
        return 0;
    }
    
    rMUT.nucREF = cREF[0];
    strncpy(rMUT.nucALT, cALT, 8);
    if ( ! vecDNK[*ind_Xro].testValidDNK(rMUT.nucPos, rMUT.nucREF ) )
        return -1;
    
    char *pAD = cFORMAT;
    int nCol =0;
    while ( *pAD ) {
        if ( *pAD == 'A' && *(pAD+1) == 'D')
            break;
        if ( *pAD == ':' )
            nCol++;
        pAD++;
    }
    if ( ! (*pAD)  )    {
        printf("\tLine %d:: Not found AD_field at FORMAT '%s'\n\t'%s'\n", recNum+1, cFORMAT, pRec);
        return 0;
    }
    
    pAD = N845_2;
    int n=0;
    while (  *pAD ) {
        if ( *pAD == ':' )
            n++;
        pAD++;
        if ( n == nCol )
            break;
    }
    if ( ! (*pAD)  )    {
        printf("\tLine %d:: Not found AD_values\n\tFORMAT '%s'\n\t       '%s'\n", recNum+1, cFORMAT, N845_2);
        return 0;
    }
    
    int counts[4]={0,0,0,0};
    int nAD=0;
    int sumAD =0;
    while ( *pAD && nAD < 4) {
        counts[nAD] = atoi(pAD);
        sumAD += counts[nAD];
        while ( isdigit(*pAD) ) pAD++;
        if ( *pAD == ',' )  {
            pAD++;
            nAD++;
            continue;
        }
        break;
    }
    for ( n=1; n<=nAD; n++ )
        rMUT.calcVAF[n-1] = (float)counts[n] / (float)sumAD;
        
    return 1;
}
/////////////////////////////////////////////////////////////////////////////////

int parsINFOrec( int recNum, char *pRec, int *ind_Xro, MUTANT &rMUT )
{
    char cREF[1024]="?",  cALT[1024]="?", cXr[XRO_ID_SIZE]="?";
    char cINFO[2048]="?";
//    int XroID;

//                 #CHROM POS  ID   REF  ALT QUAL FILTER INFO
    if ( sscanf(pRec, "%s\t%d\t%*s\t%s\t%s\t%*s\t%*s\t%s\n", cXr, &rMUT.nucPos, cREF, cALT, cINFO) != 5 ) {
        printf("Line %d:: Failed to parse: Unknown format: [X=%s, Pos=%d, Ref=%s, Alt=%s\n\t'%s'\n",
                   recNum+1,  cXr, rMUT.nucPos, cREF, cALT, pRec);
//            fclose(f_Mutas);
            return -1;
    }

    *ind_Xro = findXroByID(cXr);
    if ( *ind_Xro  < 0 ) {
        printf("\tLine %d:: Unknown chromo_name %s\n\t'%s'\n", recNum+1, cXr, pRec);
        return -1;
    }
    if ( strlen(cREF) > 1  || strlen(cALT) > 1 ) {
//        printf("\tLine %d IGNORED: Insertion: %s %d REF'%s' ALT'%s'\n", recNum+1, cXr, rMUT.nucPos, cREF, cALT);
        return 0;
    }
    
    rMUT.nucREF = cREF[0];
    strncpy(rMUT.nucALT, cALT, 8);
    
    if ( !vecDNK[*ind_Xro].testValidDNK(rMUT.nucPos, rMUT.nucREF ) )
        return -1;
    
    int alt_count, ref_count;
    char *pCount;
    if ( ! ( pCount= strstr(cINFO, "t_alt_count=") ) )
        return 1;
    while ( *pCount != '=' ) pCount++;
        pCount++;  // =
    alt_count = atoi(pCount);

    if ( ! ( pCount= strstr(cINFO, "t_ref_count=") ) )
        return 1;
    while ( *pCount != '=' ) pCount++;
        pCount++;  // =
    ref_count = atoi(pCount);

    if ( ref_count > 0 )
        rMUT.calcVAF[0] = (float)alt_count / (float)ref_count;
    
    return 1;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

