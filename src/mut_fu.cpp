
//#include <dirent.h>
#include <string.h>
#include <time.h>
#include <zlib.h>

#include "cmain.h"
#include "xrosoma.h"

#define VCF_REC_SIZE 8192

extern PROGARGS ArgKit;
extern vector < XROSOMA > vecDNK;
extern FILE *Ftrace;

int parsFORMATrec( int recNum, char *pB, int *ind_Xro, MUTANT &rMUT );
int parsINFOrec( int recNum, char *pB, int *ind_Xro, MUTANT &rMUT );

////////////////////////////////////////////////////////////////////////

bool lesser_MUT ( const MUTANT &x1, const MUTANT &x2 )
{
    return x1.nucPos < x2.nucPos;
}
/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

class VCF_IN_FILE {
private:
    //    std::ifstream plainFile;
    FILE *vcfFile = NULL;
    gzFile gzFilePtr = nullptr;
    //    char buffer[4096];
    std::string ErrTxt;
    bool useGz = false;
    bool opened = false;
    
    static bool isGzFile(const std::string& filename) {
        return filename.size() >= 3 && filename.compare(filename.size() - 3, 3, ".gz") == 0;
    }
    
public:
    VCF_IN_FILE() = default;  // Default constructor
    
//--------- Open the file called later
    void open(const std::string& filename) {
        useGz = isGzFile(filename);
        if (useGz) {
            gzFilePtr = gzopen(filename.c_str(), "rb");
            if (!gzFilePtr) {
                ErrTxt = "Failed to open : " + filename;
                throw std::runtime_error(ErrTxt);
            }
        } else {
            //            plainFile.open(filename);
            if ( ! (vcfFile = fopen(filename.c_str(), "rb")) ) {
                //            if (!plainFile.is_open()) {
                ErrTxt = "Failed to open : " + filename;
                throw std::runtime_error(ErrTxt);
            }
        }
        opened = true;
    }
//-------------------------
    
    bool getline(char * pBuff) {
        if (!opened) {
            throw std::runtime_error("VCF_file  not opened");
        }
        if (useGz) {
            if (!gzFilePtr)
                return false;
            if ( ! gzgets(gzFilePtr, pBuff, VCF_REC_SIZE) )
                return false;
            return true;
        } else {
            if ( ! vcfFile )
                return false;
            if ( fgets_ShortRec(pBuff, VCF_REC_SIZE-1, vcfFile )==0 )
                return false;
            //            if (!std::getline(plainFile, line)) return false;
            return true;
        }
    }
    //-------------------------
    
    void close() {
        if (gzFilePtr) {
            gzclose(gzFilePtr);
            gzFilePtr = nullptr;
        }
        if ( vcfFile ) {
            fclose(vcfFile);
        }
        opened = false;
        useGz = false;
    }
    
    ~VCF_IN_FILE() {
        close(); // Ensure files are closed on destruction
    }
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void  xtrSamplName ( const char  *vcf_Fpath, string &sName )    //const char  *vcf_Fname,
{
    const char *pFp = vcf_Fpath + strlen(vcf_Fpath) - 1;
    
    while ( *pFp != '/' ) {
        if ( pFp == vcf_Fpath )
            break;
        pFp--;
    }
    if ( *pFp == '/' )
        pFp++;
    
//    const char *pExt = strstr(pFp, ".vcf");
    const char *pExt = vcf_Fpath + strlen(vcf_Fpath) - 1;
    while ( *pExt != '.' ) {
        if ( pExt == vcf_Fpath )
            break;
        pExt--;
    }
    if ( *pExt == '.' )
        pExt--;
    
    sName.clear(); // ArgKit.MutSampl.clear();
    for ( const char *p=pFp; p<=pExt;  p++ )
        sName += *p; //ArgKit.MutSampl += *p;
    
//    clearSampl( );
    
    return;
}
////////////////////////////////////////////////////////////////////////

int loadVCFMutation( const char *vcf_Fpath )
{
    int vcfFORMAT=0;        // if 0: get VAF from INFO (t_alt_count, t_ref_count)
                            // if 1: get VAF from FORMAT (AD)
    MUTANT rdMUT, emptyMUT;
    int ind_Xro=-1;
    long startP = clock();
    
    printf("\nLoading mutations from '%s' ......\n", vcf_Fpath);
    fprintf(Ftrace, "\nLoading mutations from '%s' ......\n", vcf_Fpath);

    VCF_IN_FILE inputFile;
    inputFile.open( vcf_Fpath );
//    FILE *f_Mutas=fopen( vcf_Fpath, "rb");
//    if ( f_Mutas==NULL ) {
//        printf ("\nFailed to open VCF file: '%s'\n", vcf_Fpath );
//        return -1;
//    }

    xtrSamplName ( vcf_Fpath, ArgKit.MutSampl);
    clearSampl( );
        
    int cntAllRec=0;
    int cntIgnor=0;
    int cntWarning=0;
    char cBuffer[VCF_REC_SIZE];
    
    while ( 1 ) {
        memset( (void *)cBuffer, '\0', sizeof(cBuffer));
        if ( ! inputFile.getline(cBuffer) )
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
            inputFile.close();
            return -1;
        }
        if ( rdMUT.isDel() && rdMUT.isIns() )   {
            fprintf(Ftrace, "Error at REF/ALT fields: X.%s nucPos=%d '%s' -> '%s'\n",
                    vecDNK[ind_Xro].XroID.c_str(), rdMUT.nucPos, rdMUT.nucREF.c_str(), rdMUT.nucALT.c_str() );
            cntWarning++;
        }
        if ( retC == 0 )    {
            cntIgnor++;
            continue;
        }
        
//        vecDNK[ind_Xro].testValidDNK(rdMUT.nucPos, rdMUT.nucREF );   do it at parser
//        cAPO = ( vecDNK[ind_Xro].APOtest(rdMUT.nucPos, rdMUT.nucREF, crdMUT.nucALT) > 0 ) ? 1 : 0;
        if ( vecDNK[ind_Xro].vMutAPO.empty() || vecDNK[ind_Xro].vMutAPO.back().nucPos <= rdMUT.nucPos )
            vecDNK[ind_Xro].vMutAPO.push_back( rdMUT );
        else {
            vector <MUTANT>::iterator Iter =
            lower_bound( vecDNK[ind_Xro].vMutAPO.begin( ), vecDNK[ind_Xro].vMutAPO.end( ),
                        rdMUT,  lesser_MUT);
            if ( Iter->nucPos != rdMUT.nucPos ) {
                vecDNK[ind_Xro].vMutAPO.insert( Iter, rdMUT );
                continue;
            }
            printf("\tLine %d IGNORED: DublPos of Mutation: %s %d REF'%s' ALT'%s'\n", cntAllRec,
                       vecDNK[ind_Xro].XroID.c_str(), rdMUT.nucPos, rdMUT.nucREF.c_str(), rdMUT.nucALT.c_str());
            fprintf(Ftrace, "\tLine %d IGNORED: DublPos of Mutation: %s %d REF'%s' ALT'%s'\n", cntAllRec,
                       vecDNK[ind_Xro].XroID.c_str(), rdMUT.nucPos, rdMUT.nucREF.c_str(), rdMUT.nucALT.c_str());
            cntIgnor++;
            cntWarning++;
        }

    }
    inputFile.close();
    
    int cntMut=0;
    for ( int nX=0; nX<vecDNK.size(); nX++ )    {
        cntMut += vecDNK[nX].vMutAPO.size();
        if ( vecDNK[nX].maxMUTsize < vecDNK[nX].vMutAPO.size() )
            vecDNK[nX].maxMUTsize = vecDNK[nX].vMutAPO.size() + 1;
    }
    
    long stopP = clock();
    double duration = (double)(stopP - startP) / CLOCKS_PER_SEC;
    if ( cntWarning > 0 )
        printf("\n!!! CHECK VCF_file. Finded %d errors(Look file 'mutrace.txt')\n\n", cntWarning );
    
    printf("=== Loaded %d mutations; ignored=%d   dT=%5.2f\n", cntMut, cntIgnor, duration );
    fprintf(Ftrace,"=== Loaded %d mutations; ignored=%d   dT=%5.2f\n", cntMut, cntIgnor, duration );
    
    return cntMut;
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

void MUTANT:: identiMuType(  )
{
    int cntNuc=0;
    
    if ( nucREF.size() == 1 ) {
        if ( nucALT.size() == 1 )   {
            setSingl();
            return;
        }
        for ( int n=0; n<nucALT.size(); n++)     {
            if ( nucALT[n]==',' )
                cntNuc--;
            else
                cntNuc++;
        }
        if ( cntNuc > 1 )
            setIns();
        else
            setSingl();
        return;
    }
//----- sREF.size() > 1
    if ( nucALT.size() == 1 )
        setDel();
    else    {
        setDel();
        setIns();
    }
    
    return;
}
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
        printf("\tLine %d:: Unknown CHROM_name %s\n\t'%s'\n", recNum+1, cXr, pRec);
        return -1;
    }
    rMUT.nucREF = cREF;
    rMUT.nucALT = cALT;
    rMUT.identiMuType( );
    if ( rMUT.isSingl() ) {
        if (  ! ArgKit.isArg_SBS() )       // set only param '-id'
            return 0;
    }
    else    {
        if (  ! ArgKit.isArg_ID() )         // unset param '-id'
           return 0;
    }
    
    if ( ! vecDNK[*ind_Xro].testValidDNK(rMUT.nucPos, rMUT.nucREF[0] ) )
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
//        printf("\tLine %d:: Not found AD_field at FORMAT '%s'\n\t'%s'\n", recNum+1, cFORMAT, pRec);
        return 1;   // 0
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
//        printf("\tLine %d:: Not found AD_values\n\tFORMAT '%s'\n\t       '%s'\n", recNum+1, cFORMAT, N845_2);
        return 1;   // 0
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
    rMUT.nucREF = cREF;
    rMUT.nucALT = cALT;
    rMUT.identiMuType( );
//    if ( rMUT.mutType != MUT_SINGL )
//        return 0;
    if ( rMUT.isSingl() ) {
        if (  ! ArgKit.isArg_SBS() )       // set only param '-id'
            return 0;
    }
    else    {
        if (  ! ArgKit.isArg_ID() )         // unset param '-id'
            return 0;
    }
    
    if ( !vecDNK[*ind_Xro].testValidDNK(rMUT.nucPos, rMUT.nucREF[0] ) )
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
        rMUT.calcVAF[0] = (float)alt_count / (float)(ref_count+alt_count);
    
    return 1;
}
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////

