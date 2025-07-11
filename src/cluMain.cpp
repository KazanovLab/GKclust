// dnk.cpp : Defines the entry point for the console application.
//
//#include <iostream>
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include "cmain.h"
#include "xrosoma.h"


long XROSOMA::maxXsize = 0;
long XROSOMA::maxMUTsize = 0;
int XROSOMA::chrIDmode =-1; 
extern vector < XROSOMA > vecDNK;

PROGARGS ArgKit;
int procArg( int argc, char* argv[] );

#ifdef DEBUG_TRACE_TSLD
FILE *tsld=NULL;
#endif

/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
//  1-2: -g Fpath for Xro_file
//  3-4: -sfd Folder for set Mutation_Files (.vcf)
//        -s short(Clust_only)
//        -f full (Mut&Clust)
//        -d debugg (Track, dist_Mut, dist_Clust)
//  5-6: -o Folder for output data       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  7-8: -t % (0.1 - 100) ] default = 1
    
    double duration;
    int errCnt=0;
    
    if ( procArg( argc, argv) < 0 )
        return -1;
    
#ifdef DEBUG_TRACE_TSLD
    char buffr[1024];
    sprintf(buffr, "%smutrace.txt", ArgKit.OUTdir.c_str() );
    tsld = fopen(buffr, "w");
#endif
    
    
    if ( LoadHuGen( ArgKit.HUGpath.c_str() ) <= 0 )
        return -1;
    //    testXsize( );
    
    clock_t start = clock();
    DIR *dir_MUT = opendir(ArgKit.MUTdir.c_str());
    if (dir_MUT == NULL) {
        printf ("<input_directory> containing VCF files does not exist: '%s'\n", ArgKit.MUTdir.c_str());
        return -1;
    }
    
    struct dirent *f_MUT;
    
    while ( (f_MUT=readdir(dir_MUT)) != NULL ) {
        //        if ( f_MUT->d_type != DT_REG )
        //            continue;
        //        if ( f_MUT->d_name[0]=='.' )
        //            continue;
        //
        if ( ! strstr( f_MUT->d_name, ".vcf") )
            continue;
        
        clock_t startS = clock();
        if ( loadVCFMutation( f_MUT->d_name ) < 0 )
            continue;
        
        if ( ArgKit.openOutFiles ( f_MUT->d_name ) < 0 )
            continue;
        
        srand(SRAND_VALUE);
        memRezerv( );
//        printMutRezerv( );
        
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
            errCnt += vecDNK[nX].testMutUnClust( );
            
            clock_t endT = clock();
            duration = (double)(endT - begT) / CLOCKS_PER_SEC;
            printf("%s: Muts %ld\t Clusts %ld\t overClusts %ld\t dT=%5.2f\n", vecDNK[nX].XroID.c_str(),
                   vecDNK[nX].vMutAPO.size(), vecDNK[nX].vClust.size(),vecDNK[nX].vAggrC.size(),
                   duration );
        }
//        printMutRezerv( );
        
        if ( GET_ARG_S(ArgKit.argTAG)  )
            printCluMut( 1, ArgKit.foutClust );
        if ( GET_ARG_F(ArgKit.argTAG)  )
            printCluMut( 0, ArgKit.foutMu_Clu );
        if ( GET_ARG_D(ArgKit.argTAG)  )
            printCluTrace(  );
        
        clock_t stopS = clock();
        duration = (double)(stopS - startS) / CLOCKS_PER_SEC;
        printf("\tProcessed %d mutations in %5.2f seconds.\n", sumMut, duration);
    }
    closedir(dir_MUT);
//    printMutRezerv( );
    
    if ( errCnt > 0 )
        printf("\nDetecned %d ERRORS. Look  file 'mutrace.txt'\n", errCnt);
    
    clock_t finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("\nTotal execution time: %5.2f seconds.\n", duration);

    
    return 0;

}
/////////////////////////////////////////////////////////////////////////

int procArg( int argc, char* argv[] )
{
//  1-2: -g Fpath for Xro_file
//  3-4: -sfd Folder for set Mutation_Files (.vcf)
//        -s short(Clust_only)
//        -f full (Mut&Clust)
//        -d debugg (Track, dist_Mut, dist_Clust)
//  5-6: -o Folder for output data       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  7-8: -t % (0.1 - 100) ] default = 1
    
    
    if ( argc < 7 ) {
        printf ("\nUsage: sbsclust (-s | -f | -sf) <input_directory> -o <output_directory> -g <reference_genome> [-t <p-value>]\n\n");
        printf ("Options:\n");
        printf ("  -s \tGenerate short output (only clustered mutations)\n");
        printf ("  -f \tGenerate full output (all mutations from the input VCF files)\n");
        printf ("  -sf\tGenerate both short and full outputs simultaneously\n");
        printf ("  <input_directory>\tDirectory containing VCF files (required).\n");
        printf ("  -o <output_directory>\tDirectory to save output files (required).\n");
        printf ("  -g <reference_genome>\tPath to the reference genome (required).\n");
        printf ("  -t <p-value>\tP-value for determining the between-mutation distance threshold (default 0.01).\n\n");
        printf ("Examples:\n");
        printf ("  Short format output:\n");
        printf ("    sbsclust -s /inputdir/ -o /outputdir/ -g /refgenome/hg19.fa\n\n");
        printf ("  Full format output:\n");
        printf ("    sbsclust -f /inputdir/ -o /outputdir/ -g /refgenome/hg19.fa\n\n");
        printf ("  Both short and full output:\n");
        printf ("    sbsclust -sf /inputdir/ -o /outputdir/ -g /refgenome/hg19.fa\n\n");
        printf ("Error: Incorrect usage. Please check the command syntax and try again.\n\n");
        
        return -1;
    }
    
    FILE *fTest;
    DIR  *dTest;
    
    int nP=1;
    while ( nP+1 < argc ) {
        if ( argv[nP][0] != '-' )   {
            printf ("Invalid option: '%s'\n", argv[nP]);
            return -1;
        }
        switch (argv[nP][1]) {
            case 'o':
                ArgKit.OUTdir = argv[nP+1];
                if (ArgKit.OUTdir.back() != '/' )
                    ArgKit.OUTdir += '/';
                if ( (dTest=opendir(ArgKit.OUTdir.c_str()) ) )
                    closedir(dTest);
                else    {
                    printf ("<output directory> does not exist: '%s'\n", ArgKit.OUTdir.c_str());
                    return -1;
                }
                break;
            case 'g':
                ArgKit.HUGpath = argv[nP+1];
                if ( (fTest=fopen(ArgKit.HUGpath.c_str(), "r") ) )
                    fclose(fTest);
                else    {
                    printf ("<reference_genome> does not exist: '%s'\n", ArgKit.HUGpath.c_str());
                    return -1;
                }
                break;
            case 's':
            case 'f':
            case 'd':
                ArgKit.MUTdir = argv[nP+1];
                if (ArgKit.MUTdir.back() != '/' )
                    ArgKit.MUTdir += '/';
                if ( (dTest=opendir(ArgKit.MUTdir.c_str()) ) )
                    closedir(dTest);
                else    {
                    printf ("<input_directory> containing VCF files does not exist: '%s'\n", ArgKit.MUTdir.c_str());
                    return -1;
                }
                for ( int i=1; argv[nP][i]; i++ )   {
                    if ( argv[nP][i]=='s' )
                        SET_ARG_S(ArgKit.argTAG);
                    else
                    if ( argv[nP][i]=='f' )
                        SET_ARG_F(ArgKit.argTAG);
                    else
                    if ( argv[nP][i]=='d' )
                        SET_ARG_D(ArgKit.argTAG);
                    else    {
                        printf ("Invalid option: '%s' : '%c'\n", argv[nP], argv[nP][i]);
                        return -1;
                    }
                }
                break;
            case 't':
                strncpy(ArgKit.chPart, argv[nP+1], sizeof(ArgKit.chPart)-1);
                float fP;
                sscanf(ArgKit.chPart, "%f", &fP);
                if ( fP <= 0.001 || fP >= 0.999 )   {
                    printf ("\nInvalid <p-value> = '%s'. Expected ba number >=0.001  and  <=0.999 \n", ArgKit.chPart );
                    return -1;
                }
                
                break;
            default:
                printf ("Invalid option: '%s'\n", argv[nP]);
                return -1;
        }
        nP += 2;
        
    }
    
    if ( ArgKit.OUTdir.empty() )    {
        printf ("Not defined:  -o <output_directory> \n");
        return -1;
    }
    if ( ArgKit.HUGpath.empty() )    {
        printf ("Not defined:  -g <reference_genome> \n");
        return -1;
    }
    if ( ArgKit.MUTdir.empty() )    {
        printf ("Not defined: -s <input_directory>\n");
        return -1;
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int PROGARGS::openOutFiles ( const char *vcf_Fname )
{
    char buffr[2048];
    char inpFname[1024];
    char commaP[8];
    
    const char *pn = strstr(vcf_Fname, ".vcf");
    memset(inpFname, '\0', sizeof(inpFname) );
    strncpy(inpFname, vcf_Fname, (pn-vcf_Fname));
    
    strcpy(commaP, chPart);
    for ( int n=0; n<strlen(commaP); n++) {
        if ( commaP[n]=='.' ) {
            commaP[n] = ',';
            break;
        }
    }
    if ( GET_ARG_S(argTAG)  )    {
        if ( foutClust )
            fclose(foutClust);
        foutClust = NULL;
        sprintf(buffr, "%s%s._%s_Clust.txt", OUTdir.c_str(), inpFname, commaP);  //vecSAMPL[iSamp].SampName.c_str()
        if ( ! (foutClust=fopen(buffr, "w") ) ) {
            printf("InvOPENop '%s'\n", buffr);
            goto BadExit;
        }
    }
    
    if ( GET_ARG_F(argTAG)  )    {
        if ( foutMu_Clu )
            fclose(foutMu_Clu);
        foutMu_Clu = NULL;
        sprintf(buffr, "%s%s._%s_Clu+MU.txt", OUTdir.c_str(), inpFname, commaP);  //vecSAMPL[iSamp].SampName.c_str()
        if ( ! (foutMu_Clu=fopen(buffr, "w") ) ) {
            printf("InvOPENop '%s'\n", buffr);
            goto BadExit;
        }
    }
    
    if ( GET_ARG_D(argTAG)  )    {
        if ( foutTrace )
            fclose(foutTrace);
        foutTrace = NULL;
        if ( foutRndCl )
            fclose(foutRndCl);
        foutRndCl = NULL;
        if ( foutRndMu )
            fclose(foutRndMu);
        foutRndMu = NULL;
        
        sprintf(buffr, "%s%s._%s_Trace.txt", OUTdir.c_str(), inpFname, commaP);  //vecSAMPL[iSamp].SampName.c_str()
        if ( ! (foutTrace=fopen(buffr, "w") ) ) {
            printf("Failed to create output file:  '%s'\n", buffr);
            goto BadExit;
        }
        sprintf(buffr, "%s%s._%s_RndCl.txt", OUTdir.c_str(), inpFname, commaP);  //vecSAMPL[iSamp].SampName.c_str()
        if ( ! (foutRndCl=fopen(buffr, "w") ) ) {
            printf("Failed to create output file: '%s'\n", buffr);
            goto BadExit;
        }
        sprintf(buffr, "%s%s._%s_RndMu.txt", OUTdir.c_str(), inpFname, commaP);  //vecSAMPL[iSamp].SampName.c_str()
        if ( ! (foutRndMu=fopen(buffr, "w") ) ) {
            printf("Failed to create output file: '%s'\n", buffr);
            goto BadExit;
        }
    }
    
    if ( false )    {
    BadExit:
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
        foutClust = NULL;
        foutMu_Clu = NULL;
        foutTrace = NULL;
        foutRndCl = NULL;
        foutRndMu = NULL;
        
        return -1;
    }
    
    return 0;
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
 dir = opendir(path);//(FOLDER_MUT_DATA);
 if (dir == NULL) {
 perror("opendir");
 return EXIT_FAILURE;
 }
 
 printf("Contents of the directory '%s':\n", path);//FOLDER_MUT_DATA);
 
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

