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

#ifdef DEBUG_TRACE_TSLD
	FILE *tsld=NULL;
#endif

long XROSOMA::maxXsize = 0;
extern vector < XROSOMA > vecDNK;
long SAMPLE::maxMUTsize = 0;
extern vector < SAMPLE >  vecSAMPL;

PROGARGS ArgKit;
int procArg( int argc, char* argv[] );

/*const char CancerID[NO_CANCER_ID][8] = { "BLCA", "BRCA", "CESC", "HNSC", "LUAD", "LUSC"};
int DNKset[NO_HUMAN_XRO][2] =
{   {1, 252811431}, {2, 246673736}, {3, 200851408}, {4, 193885138},
    {5, 183499849}, {6, 173559654}, {7, 161412159}, {8, 148455023},
    {9, 143230852}, {10,137471045}, {11,136935267}, {12,135764152},
    {13,116815249}, {14,108883191}, {15,103996213}, {16, 91645622},
    {17, 82355229}, {18, 79192724}, {19, 59973769}, {20, 63925972},
    {21, 48817551}, {22, 52037576}, {23,157488797}, {24, 60221845}
};
*/

/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
//  1-2: -g Fpath for Xro_file
//  3-4: -sfd Folder for set Mutation_Files (.vcf)
//        -s short(Clust_only)
//        -f full (Mut&Clust)
//        -d debugg (Track, dist_Mut, dist_Clust)
//  5-6: -o Folder for output data       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  7-8: -t % (0.1 - 100) ] default = 1
{
    char buffr[1024];
    double duration;

    
    if ( procArg( argc, argv) < 0 )
        return -1;

#ifdef DEBUG_TRACE_TSLD
    sprintf(buffr, "%smutrace.txt", ArgKit.OUTdir.c_str() );
    tsld = fopen(buffr, "w");
#endif
    
    
    if ( LoadHuGen( ArgKit.HUGpath.c_str() ) <= 0 )
        return -1;
//    testXsize( );
//    return 0;

    clock_t start = clock();
    DIR *dir_MUT = opendir(ArgKit.MUTdir.c_str());
    if (dir_MUT == NULL) {
        printf ("INV.opendir ='%s'\n", ArgKit.MUTdir.c_str());
        return -1;
    }
    
    struct dirent *f_MUT;
    int iSamp=0;
    while ( (f_MUT=readdir(dir_MUT)) != NULL ) {
        if ( f_MUT->d_name[0]=='.' )
            continue;
        
        clock_t startS = clock();
        if ( loadVCFMutation( f_MUT->d_name ) < 0 )
            continue;
        
        if ( ArgKit.openOutFiles ( iSamp ) < 0 )
            continue;

        int clustID = 1;
        int sumMut=0;
        for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
            vecSAMPL[iSamp].PorogMut[nX] = vecSAMPL[iSamp].calcMutPorog( nX ); 
            if (  vecSAMPL[iSamp].PorogMut[nX] > 0 )
                vecSAMPL[iSamp].findClusters( nX,  &clustID, vecSAMPL[iSamp].PorogMut[nX] );
        
            vecSAMPL[iSamp].PorogClust[nX] = vecSAMPL[iSamp].calcClustPorog( nX );
            if (  vecSAMPL[iSamp].PorogClust[nX] > 0 )
                vecSAMPL[iSamp].ClustConnector( nX,  &clustID,  vecSAMPL[iSamp].PorogClust[nX] );
            sumMut +=  vecSAMPL[iSamp].vMutAPO[nX].size();
        
        }
    
        if ( GET_ARG_S(ArgKit.argTAG)  )
            vecSAMPL[iSamp].printCluMut( 1, ArgKit.foutClust );
        if ( GET_ARG_F(ArgKit.argTAG)  )
            vecSAMPL[iSamp].printCluMut( 0, ArgKit.foutMu_Clu );
        if ( GET_ARG_D(ArgKit.argTAG)  )
            vecSAMPL[iSamp].printCluTrace(  );
    
        clock_t stopS = clock();
        duration = (double)(stopS - startS) / CLOCKS_PER_SEC;
        printf("S.%d #Mut=%d dT=%5.2fsec\n", iSamp, sumMut, duration);
    }
    closedir(dir_MUT);
        
    clock_t finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("\nAllTime dT=%5.2f sec\n", duration);
    
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
    
    
    if ( argc < 3 ) {
        printf ("\nUsage: %s  -g 'HUGpath' -o 'OUTdir' -[s][f][d] 'MUTdir' [-t 'Part(0-100)']\n", argv[0]);
        return -1;
    }
    
    FILE *fTest;
    DIR  *dTest;
    
    int nP=1;
    while ( nP+1 < argc ) {
        if ( argv[nP][0] != '-' )   {
            printf ("Inv. argID='%s'\n", argv[nP]);
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
                    printf ("Not found Output folder ='%s'\n", ArgKit.OUTdir.c_str());
                    return -1;
                }
                break;
            case 'g':
                ArgKit.HUGpath = argv[nP+1];
                if ( (fTest=fopen(ArgKit.HUGpath.c_str(), "r") ) )
                    fclose(fTest);
                else    {
                    printf ("FILE not found ='%s'\n", ArgKit.HUGpath.c_str());
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
                    printf ("Not found Mutations folder ='%s'\n", ArgKit.MUTdir.c_str());
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
                        printf ("Inv. argID='%s' : '%c'\n", argv[nP], argv[nP][i]);
                        return -1;
                    }
                }
                break;
            case 't':
                strncpy(ArgKit.chPart, argv[nP+1], sizeof(ArgKit.chPart)-1);
                float fP;
                sscanf(ArgKit.chPart, "%f", &fP);
                if ( fP <= 0.05 || fP >= 100 )   {
                    printf ("\nBAD percentValue = '%s'  MustBe >0.0  &  <100.0 \n", ArgKit.chPart );
                    return -1;
                }
                
                break;
            default:
                printf ("Inv. argID='%s'\n", argv[nP]);
                return -1;
        }
        nP += 2;
        
    }
    
    if ( ArgKit.OUTdir.empty() )    {
        printf ("Not defined Output folder (-o path)\n");
        return -1;
    }
    if ( ArgKit.HUGpath.empty() )    {
        printf ("Not defined GenomeFile (-g filename)\n");
        return -1;
    }
    if ( ArgKit.MUTdir.empty() )    {
        printf ("Not defined Mutations folder (-s path)\n");
        return -1;
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int PROGARGS::openOutFiles ( int iSamp )
{
    char buffr[2048];
    
    if ( GET_ARG_S(argTAG)  )    {
        if ( foutClust )
            fclose(foutClust);
        foutClust = NULL;
        sprintf(buffr, "%s%s._%s_Clust.txt", OUTdir.c_str(), vecSAMPL[iSamp].SampName.c_str(), chPart);
        if ( ! (foutClust=fopen(buffr, "w") ) ) {
            printf("InvOPENop '%s'\n", buffr);
            goto BadExit;
        }
    }
    
    if ( GET_ARG_F(argTAG)  )    {
        if ( foutMu_Clu )
            fclose(foutMu_Clu);
        foutMu_Clu = NULL;
        sprintf(buffr, "%s%s._%s_Clu+MU.txt", OUTdir.c_str(), vecSAMPL[iSamp].SampName.c_str(), chPart);
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
        
        sprintf(buffr, "%s%s._%s_Trace.txt", OUTdir.c_str(), vecSAMPL[iSamp].SampName.c_str(), chPart);
        if ( ! (foutTrace=fopen(buffr, "w") ) ) {
            printf("InvOPENop '%s'\n", buffr);
            goto BadExit;
        }
        sprintf(buffr, "%s%s._%s_RndCl.txt", OUTdir.c_str(), vecSAMPL[iSamp].SampName.c_str(), chPart);
        if ( ! (foutRndCl=fopen(buffr, "w") ) ) {
            printf("InvOPENop '%s'\n", buffr);
            goto BadExit;
        }
        sprintf(buffr, "%s%s._%s_RndMu.txt", OUTdir.c_str(), vecSAMPL[iSamp].SampName.c_str(), chPart);
        if ( ! (foutRndMu=fopen(buffr, "w") ) ) {
            printf("InvOPENop '%s'\n", buffr);
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
    srand(SRAND_VALUE);
    
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

