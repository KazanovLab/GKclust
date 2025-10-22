// Based at mClust_patch progect
//      Добавлена возможность обработки вставок в мутациях
//      Добавлена выдача 'mini'
//      Добавлена обработка статистики по геному и выдача 'stat'
//      Добавлена обработка GZ_файлов
//#include <iostream>
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/stat.h>


#include "cmain.h"
#include "xrosoma.h"


long XROSOMA::maxXsize = 0;
long XROSOMA::maxMUTsize = 0;
int XROSOMA::genCLUSTsize=0;
int XROSOMA::chrIDmode =-1; 
extern vector < XROSOMA > vecDNK;

PROGARGS ArgKit;
int procArg( int argc, char* argv[] );

FILE *Ftrace=NULL;


/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    double duration;
    int errCnt=0;
    char Fpath[1024];
    
    DIR *dir_MUT=NULL;
    struct dirent *f_MUT;
    FILE *fMutSet=NULL;
    
    clock_t start = clock();

    if ( ArgKit.procArg( argc, argv) < 0 )
        return -1;
    
    sprintf(Fpath, "%smutrace.txt", ArgKit.OUTdir.c_str() );
    Ftrace = fopen(Fpath, "w");
    
    if ( LoadHuGen( ArgKit.HUGpath.c_str() ) <= 0 )
        return -1;
    
    if ( ArgKit.MUTdir.size() > 0 ) {           //==============  -i option ================
        dir_MUT = opendir(ArgKit.MUTdir.c_str());
        if ( dir_MUT == NULL ) {
            printf ("Cannt OPEN <input_directory> containing VCF files: '%s'\n", ArgKit.MUTdir.c_str());
            return -1;
        }
    }
    else    {                                   //==============  -l option ==================
        fMutSet = fopen(ArgKit.Mutlist.c_str(), "r");
        if ( fMutSet==NULL )    {
            printf ("Cannt OPEN File containing list VCF_files: '%s'\n", ArgKit.Mutlist.c_str() );
            return -1;
        }
    }
    
    while ( getNextMutFile( dir_MUT, fMutSet, Fpath ) ) {
        clock_t startS = clock();
        if ( loadVCFMutation( Fpath ) < 0 )
            continue;
        
        if ( ArgKit.openOutFiles (  ) < 0 ) 
            continue;
        
        srand(SRAND_VALUE);
        memRezerv( );
        
        printf("Processing :   ................\n");
        int clustID = 1;
        int sumMut=0;
        errCnt=0;
        for ( int nX=0; nX<vecDNK.size(); nX++ ) {
            if ( vecDNK[nX].vMutAPO.size() == 0 )
                continue;
            clock_t begT = clock();
            if (  vecDNK[nX].calcMutPorog( ) > 0 )
                vecDNK[nX].findClusters( &clustID );
            
            if (  vecDNK[nX].calcClustPorog( ) > 0 )
                vecDNK[nX].ClustConnector( &clustID);
            
            sumMut +=  vecDNK[nX].vMutAPO.size();
            errCnt += vecDNK[nX].checkMutUnClust( );
            vecDNK[0].genCLUSTsize += (int)vecDNK[nX].vClust.size();
            
            clock_t endT = clock();
            duration = (double)(endT - begT) / CLOCKS_PER_SEC;
            printf("%s: Mutations %ld\t Initial clusters %ld\t Merged clusters %ld\t dT=%5.2f\n", vecDNK[nX].XroID.c_str(),
                   vecDNK[nX].vMutAPO.size(), vecDNK[nX].vClust.size(),vecDNK[nX].vAggrC.size(),
                   duration );
        }
        
        if ( ArgKit.isArg_S() )
            printCluMut( 1, ArgKit.foutClust );
        if ( ArgKit.isArg_F() )
            printCluMut( 0, ArgKit.foutMu_Clu );
        if ( ArgKit.isArg_D() )
            printCluTrace(  );
        if ( ArgKit.isArg_M()  )
            printCluMini(  );
        if ( ArgKit.isArg_STAT()  )
            printCluStat(  );
        
        clock_t stopS = clock();
        duration = (double)(stopS - startS) / CLOCKS_PER_SEC;
        printf("\tProcessed %d mutations in %5.2f seconds.\n", sumMut, duration);
    }
    if ( dir_MUT )
        closedir(dir_MUT);
    if ( fMutSet )
        fclose(fMutSet);
    
    ArgKit.closeOutFiles (  );
    
    if ( errCnt > 0 )
        printf("\nDetecned %d ERRORS. Look  file 'mutrace.txt'\n", errCnt);
    
    clock_t finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("\nTotal execution time: %5.2f seconds.\n", duration);
    fprintf(Ftrace,"\nTotal execution time: %5.2f seconds.\n", duration);
    fclose (Ftrace);
    
    return 0;

}
/////////////////////////////////////////////////////////////////////////

int getNextMutFile( DIR *dir_MUT, FILE *Flist_MUT, char *F_path )
{
    struct dirent *entry;
    if ( dir_MUT ) {
        while ( (entry=readdir(dir_MUT)) != NULL )  {
            if ( ! strstr( entry->d_name, ".vcf") )
                continue;
            strcpy(F_path, ArgKit.MUTdir.c_str() );
            strcat(F_path, entry->d_name);
            if ( ! is_file( F_path ) )
                continue;
            return 1;
        }
        return 0;
    }
    char buff[1024];
    char *pb;
    while ( fgets(buff, sizeof(buff)-1, Flist_MUT ) ) {
        if ( ! strstr(buff, ".vcf") )
            continue;
        pb = &buff[strlen(buff)-1];
        if ( *pb=='\n' )
            *pb = '\0';
        pb = buff;
        while ( *pb && *pb <= ' ' )  pb++;
        if ( ! (*pb) )
            continue;
        if ( *pb == '#' )
            continue;
        if ( ! is_file( pb ) )
            continue;
        
        strcpy(F_path, pb);
        return 1;
    }
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int PROGARGS:: procArg( int argc, char* argv[] )
{
/*
Usage:
      gkclust (-s | -f | -c | -a) (-i <input_directory> | -p <path_list>) -o <output_directory> -g <reference_genome> [-t <p_value>]
    
    At least one of (-s | -f | -c | -a) must be specified.
    At least one of (-i | -p) must be specified.
    
Options:
      -s                     Generate short output (only clustered mutations).
      -f                     Generate full output (all mutations from the input VCF files).
      -c                     Generate only a list of detected clusters.
      -a                     Generate a file with statistics for all processed VCFs.
      -i <input_directory>   Directory containing VCF files.
      -p <path_list>         File containing a list of VCF file paths.
          -b                     Consider only single-base substitutions (default).
          -n                     Consider indels.
          -o <output_directory>  Directory to save output files. (required)
          -g <reference_genome>  Path to the reference genome. (required)
          -t <p_value>           P-value for determining the between-mutation distance threshold (default: 0.01).
            
            Examples:
              Short-format output:
                gkclust -s -i /inputdir/ -o /outputdir/ -g /refgenome/hg19.fa
            
              Full-format output:
                gkclust -f -i /inputdir/ -o /outputdir/ -g /refgenome/hg19.fa
*/
    if ( argc < 7 ) {
        printf ("\nUsage:\n gkclust (-s | -f | -c | -a) (-i <input_directory> | -p <path_list>) -o <output_directory> -g <reference_genome> [-t <p_value>]\n\n");
        printf ("At least one of (-s | -f | -c | -a) must be specified.\n");
        printf ("At least one of (-i | -p) must be specified.\n\n");
        printf ("Options:\n");
        printf ("  -s \tGenerate short output (only clustered mutations).\n");
        printf ("  -f \tGenerate full output (all mutations from the input VCF files).\n");
        printf ("  -c \tGenerate only a list of detected clusters.\n");
        printf ("  -a \tGenerate a file with statistics for all processed VCFs.\n");
        
        printf ("  -i <input_directory>\tDirectory containing VCF files.\n");
        printf ("  -p <path_list>      \tFile containing a list of VCF file paths.\n");
        
        printf ("  -b \tConsider only single-base substitutions (default).\n");
        printf ("  -n \tConsider indels.\n");
        
        printf ("  -o <output_directory>\tDirectory to save output files (required).\n");
        printf ("  -g <reference_genome>\tPath to the reference genome (required).\n");
        printf ("  -t <p-value>         \tP-value for determining the between-mutation distance threshold (default 0.01).\n\n");
        printf ("Examples:\n");
        printf ("  Short format output:\n");
        printf ("    gkclust -s -i /inputdir/ -o /outputdir/ -g /refgenome/hg19.fa\n\n");
        printf ("  Full format output:\n");
        printf ("    gkclust -f -i /inputdir/ -o /outputdir/ -g /refgenome/hg19.fa\n\n");

//        printf ("Error: Incorrect usage. Please check the command syntax and try again.\n\n");
        
        return -1;
    }
    
    char what[16];
//  char parList[] = "-g@0 -o@1 -i@2 -l@3 -s@4 -f@5 -d@6 -m@7 -sbs@8 -id@9 -stat@10 -t@11 ";
    char parList[] = "-g@0 -o@1 -i@2 -p@3 -s@4 -f@5 -d@6 -c@7 -b@8 -n@9 -a@10 -t@11 ";
    int parmN;
    
//    char buffr[4096];
    
    int nP=1;
    while ( nP+1 < argc ) {
        snprintf(what, sizeof(what)-4, "%s@", argv[nP] );
        char *pL = strstr(parList, what);
        if ( ! pL ) {
            parmN = -1;
//            return -1;
        }
        else {
                while ( *pL != '@' ) pL++;
                pL++;
                parmN = atoi(pL);
        }
        switch ( parmN )    {
            case 0:
                if ( argv[nP+1][0] == '-' ) // skipped <reference_genome>
                    break;
                HUGpath = argv[nP+1];
                if ( ! is_file(argv[nP+1]) )    {
                    printf ("<reference_genome> does not exist: '%s'\n", HUGpath.c_str());
                    return -1;
                }
                nP++;
                break;
             case 1:
                if ( argv[nP+1][0] == '-' ) // skipped <output directory>
                    break;
                OUTdir = argv[nP+1];
                if (OUTdir.back() != '/' )
                    OUTdir += '/';
                if ( ! is_dir(argv[nP+1]) )    {
                    printf ("<output directory> does not exist: '%s'\n", OUTdir.c_str());
                    return -1;
                }
                nP++;
                break;
            case 2:
                if ( argv[nP+1][0] == '-' ) // skipped <input_directory>
                    break;
                MUTdir = argv[nP+1];
                if (MUTdir.back() != '/' )
                    MUTdir += '/';
                if ( ! is_dir(argv[nP+1]) )    {
                    printf ("<input_directory> containing VCF files does not exist: '%s'\n", MUTdir.c_str());
                    return -1;
                }
                nP++;
                break;
            case 3:
                if ( argv[nP+1][0] == '-' ) // skipped <File_list>
                    break;
                Mutlist = argv[nP+1];
                if ( ! is_file(argv[nP+1]) )    {
                    printf ("File containing list VCF_files does not exist: '%s'\n", Mutlist.c_str());
                    return -1;
                }
                nP++;
                break;
            case 4:
                argTAG |= _ARG_S;
                break;
            case 5:
                argTAG |= _ARG_F;
                break;
            case 6:
                argTAG |= _ARG_D;
                break;
            case 7:
                argTAG |= _ARG_M;
                break;
            case 8:
                argTAG |= _ARG_SBS;
                break;
            case 9:
                argTAG |= _ARG_ID;
                break;
            case 10:        //stat
                argTAG |= _ARG_STAT;
                break;
            case 11:
//                strncpy(chPart, argv[nP+1], sizeof(chPart)-1);
                snprintf(chPart, sizeof(chPart), "%s", argv[nP+1]);
                float fP;
                sscanf(chPart, "%f", &fP);
                if ( fP <= 0.001 || fP >= 0.999 )   {
                    printf ("\nInvalid <p-value> = '%s'. Expected ba number >=0.001  and  <=0.999 \n", chPart );
                    return -1;
                }
                nP++;
                break;
            default:
                printf ("Invalid option: '%s'\n", argv[nP]);
                break;
        }
        nP++;
    }
    
    if ( isArg_SBS() || isArg_ID()  )
    { }
    else
         argTAG |= _ARG_SBS;  // default set _SBS
    
    if ( isArg_S() || isArg_F() || isArg_D() || isArg_M() || isArg_STAT() )
    { }
    else
        argTAG |= _ARG_S;  // default set _S
    
    if ( HUGpath.empty() )    {
        printf ("Not defined:  -g <reference_genome> \n");
        return -1;
    }
    if ( OUTdir.empty() )    {
        printf ("Not defined:  -o <output_directory> \n");
        return -1;
    }
   if ( MUTdir.empty() && Mutlist.empty() )    {
        printf ("Not defined: -i <input_directory>\n         or: -l <list_file>\n ");
        return -1;
    }
    if ( ! MUTdir.empty() && ! Mutlist.empty() )    {
        printf ("Cannt use options  '-i' and '-l' together\n ");
        return -1;
    }
    return nP;
}
/////////////////////////////////////////////////////////////////////////

int PROGARGS::openOutFiles ( )//const char *vcf_Fname )
{
    char commaP[24];
    string strFpath;
    string strFN = OUTdir + MutSampl + "+";
    
    strcpy(commaP, chPart);
    for ( int n=0; n<strlen(commaP); n++) {
        if ( commaP[n]=='.' ) {
            commaP[n] = ',';
            break;
        }
    }
    strFN += commaP;
    
    if ( isArg_SBS() && isArg_ID() )  {}
    else
        if ( isArg_SBS() )
            strFN +=  "_sbs";
        else
            strFN += "_indel";     //"_dl";

//  ------------------------
    if ( isArg_M() )    {
        if ( foutMini )
            fclose(foutMini);
        strFpath = strFN + "_clusterlist.txt";     //"_Mini.txt";
        if ( ! (foutMini=fopen(strFpath.c_str(), "w") ) ) {
            printf("InvOPENop '%s'\n", strFpath.c_str());
            goto BadExit;
        }
    }
    //  ------------------------
    if ( isArg_S() )    {
        if ( foutClust )
            fclose(foutClust);
        strFpath = strFN + "_short.txt";    //"_Clust.txt";
        if ( ! (foutClust=fopen(strFpath.c_str(), "w") ) ) {
            printf("InvOPENop '%s'\n", strFpath.c_str());
            goto BadExit;
        }
    }
//  ------------------------
    if ( isArg_F() )    {
        if ( foutMu_Clu )
            fclose(foutMu_Clu);
        strFpath = strFN + "_full.txt";   //"_CluMut.txt";
        if ( ! (foutMu_Clu=fopen(strFpath.c_str(), "w") ) ) {
            printf("InvOPENop '%s'\n", strFpath.c_str());
            goto BadExit;
        }
    }
//  ------------------------
    if ( isArg_D() )    {
        if ( foutTrace )
            fclose(foutTrace);
        strFpath = strFN + "_Trace.txt";
        if ( ! (foutTrace=fopen(strFpath.c_str(), "w") ) ) {
            printf("Failed to create output file:  '%s'\n", strFpath.c_str() );
            goto BadExit;
        }
//  -----------
        if ( foutRndCl )
            fclose(foutRndCl);
        strFpath = strFN + "_RndCl.txt";
        if ( ! (foutRndCl=fopen(strFpath.c_str(), "w") ) ) {
            printf("Failed to create output file: '%s'\n", strFpath.c_str() );
            goto BadExit;
        }
//  -----------
        if ( foutRndMu )
            fclose(foutRndMu);
        strFpath = strFN + "_RndMu.txt";
        if ( ! (foutRndMu=fopen(strFpath.c_str(), "w") ) ) {
            printf("Failed to create output file: '%s'\n", strFpath.c_str() );
            goto BadExit;
        }
    }
    //  ------------------------
    if ( isArg_STAT() )    {
        HStatpath = OUTdir + HUGname + "+" + commaP;
        if ( isArg_SBS() && isArg_ID() )  {}
        else
            if ( isArg_SBS() )
                HStatpath +=  "_sbs";
            else
                HStatpath += "_indel";
            HStatpath += "_briefstatt.txt";   // "_Stat.txt";
    }
/*    if ( isArg_STAT() )    {
        if ( foutStat )
            fclose(foutStat);
        strFpath = strFN + "_Stat.txt";
        if ( is_file(strFpath.c_str() ) )
            foutStat=fopen(strFpath.c_str(), "a");      // newStatFile == 0
        else    {
            foutStat=fopen(strFpath.c_str(), "w");
            newStatFile = 1;
        }
        if ( ! foutStat ) {
            printf("InvOPENop '%s'\n", strFpath.c_str());
            goto BadExit;
        }
    }
 */

    if ( false )    {
    BadExit:
        closeOutFiles (  );
        return -1;
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

bool is_dir(const char  *path ) {
    struct stat statbuf;
    if (stat(path, &statbuf) != 0) {
        return false; // Ошибка при получении информации о файле
    }
    return S_ISDIR(statbuf.st_mode);
}
//////////////////////////////////////

bool is_file(const char *path ) {
    struct stat statbuf;
    if (stat(path, &statbuf) != 0) {
        return false; // Ошибка при получении информации о файле
    }
    return S_ISREG(statbuf.st_mode);
}
/////////////////////////////////////////////////////////////////////////
/*
int is_directory(const char *path) {
    struct stat path_stat;
    int rc;
    rc  = stat(path, &path_stat);
    if (rc != 0) {
        return 0; // ошибка доступа
    }
    return S_ISDIR(path_stat.st_mode); // 1 если директория, иначе 0
}
 */
/////////////////////////////////////////////////////////////////////////

void PROGARGS::closeOutFiles (  )
{
    if ( foutClust )
        fclose(foutClust);
    if ( foutMu_Clu )
        fclose(foutMu_Clu);
    if ( foutTrace )
        fclose(foutTrace);
    if ( foutRndCl )
        fclose(foutRndCl);
    if ( foutRndMu )
        fclose(foutRndMu);
    if ( foutMini )
        fclose(foutMini);
//    if ( foutMini )
//        fclose(foutStat);
    foutClust = NULL;
    foutMu_Clu = NULL;
    foutTrace = NULL;
    foutRndCl = NULL;
    foutRndMu = NULL;
    foutMini = NULL;
//    foutStat = NULL;
    
    return;
}
/////////////////////////////////////////////////////////////////////////

int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in )
{
    char buff[4096];
    int lng;
    
    if ( fgets(shortRec, sizeRec, f_in)==NULL )
        return 0;
    lng = (int)strlen(shortRec);
    if ( *(shortRec+lng-1) != '\n' )
        while ( 1 ) {
            if ( fgets(buff, sizeof(buff), f_in)==NULL )
                break;
            if ( buff[strlen(buff)-1] == '\n' )
                break;
        }
    
    return lng;
}

////////////////////////////////////////////////////// //////////////////////////////
////////////////////////////////////////////////////////
/*
 //#include <stdlib.h>
 //#include <dirent.h>
 //#include <errno.h>
 
 int main( int argc, char* argv[] ) {
 //    const char *path = "/Users/gena/gena_old/mut/mut_set_orig";
 DIR *dir;
 struct dirent *entry;
 
 // Open the directory
 dir = opendir(path);
 if (dir == NULL) {
 perror("opendir");
 return EXIT_FAILURE;
 }
 
 printf("Contents of the directory '%s':\n", path);
 
 // Read entries in the directory
 while ((entry = readdir(dir)) != NULL) {
 printf("%s\n", entry->d_name);
 }
 
 // Close the directory
 if (closedir(dir) == -1) {
 perror("closedir");
 return EXIT_FAILURE;
 }
 
 return EXIT_SUCCESS;
 }
 */
/////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////// //////////////////////////////
////////////////////////////////////////////////////// //////////////////////////////

