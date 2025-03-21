
#ifndef _CMAIN_H__
#define _CMAIN_H__

#include "vector"
#include "string"
using namespace std;

//#define DEBUG_TRACE_TSLD

struct PROGARGS {
    string HUGpath;
    string MUTdir;
    string OUTdir;
    char argTAG;
    char chPart[8];
//    float PartOfMut;
    FILE * foutClust;
    FILE * foutMu_Clu;
    FILE * foutTrace;
    FILE * foutRndMu;
    FILE * foutRndCl;
    
    PROGARGS() { argTAG='\0'; strcpy (chPart,"1");      //PartOfMut = 1;
                foutClust=NULL; foutMu_Clu=NULL; foutTrace=NULL; foutRndMu=NULL; foutRndCl=NULL; };
    int openOutFiles ( int iSamp );
};

#define PROG_ARG_ID "sfdgot"
#define _ARG_S 0x01
#define _ARG_F 0x02
#define _ARG_D 0x04

#define SET_ARG_S(_TAG)    _TAG |= _ARG_S
#define GET_ARG_S(_TAG)  (((_TAG & _ARG_S)==0) ? 0 : 1 )

#define SET_ARG_F(_TAG)    _TAG |= _ARG_F
#define GET_ARG_F(_TAG)  (((_TAG & _ARG_F)==0) ? 0 : 1 )

#define SET_ARG_D(_TAG)    _TAG |= _ARG_D
#define GET_ARG_D(_TAG)  (((_TAG & _ARG_D)==0) ? 0 : 1 )


#define NO_HUMAN_XRO 24
#define NO_CANCER_ID 6

#define SRAND_VALUE 199
#define RANDOM_CYCLES 1000

int  xtrSamplName ( const char  *vcf_Fpath, char *SamplName);   //, char *Folder=NULL);
int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in );

void tst_loadMut();
void print_MutSize();
void print_N_zone( );
void testXsize( );


#endif

