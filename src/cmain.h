
#ifndef _CMAIN_H__
#define _CMAIN_H__

#include <dirent.h>
#include "vector"
#include "string"
using namespace std;

//#define DEBUG_TRACE_TSLD

#define _ARG_S 0x01
#define _ARG_F 0x02
#define _ARG_D 0x04
#define _ARG_M 0x08

#define _ARG_SBS 0x10
#define _ARG_ID  0x20
#define _ARG_STAT 0x40

struct PROGARGS {
    string HUGpath;
    string MUTdir;
    string Mutlist;
    string OUTdir;
    string MutSampl;
    string HUGname;
    string HStatpath;
    unsigned char argTAG;        // s c f d m  ::  sbs id
    char chPart[8];     // -t value
    FILE * foutClust;
    FILE * foutMu_Clu;
    FILE * foutTrace;
    FILE * foutRndMu;
    FILE * foutRndCl;
    FILE * foutMini;
    FILE * foutHStat;
    
    PROGARGS() { argTAG='\0'; strcpy (chPart,"0.01");  foutHStat = NULL;
                foutClust=NULL; foutMu_Clu=NULL; foutTrace=NULL; foutRndMu=NULL; foutRndCl=NULL; foutMini=NULL;};
    int procArg( int argc, char* argv[] );
    int openOutFiles ( );//const char *vcf_Fname );
    void closeOutFiles (  );
    bool isArg_STAT() { return ( (argTAG & _ARG_STAT) > 0); };
    bool isArg_M() { return ( (argTAG & _ARG_M) > 0); };
    bool isArg_S() { return ( (argTAG & _ARG_S) > 0); };
    bool isArg_F() { return ( (argTAG & _ARG_F) > 0); };
    bool isArg_D() { return ( (argTAG & _ARG_D) > 0); };
    
    bool isArg_SBS() { return ( (argTAG & _ARG_SBS) > 0); };
    bool isArg_ID() { return ( (argTAG & _ARG_ID) > 0); };
};


#define SRAND_VALUE 199
#define RANDOM_CYCLES 1000

void  xtrSamplName ( const char  *vcf_Fpath, string &sName );
bool is_dir(const char  *path );
bool is_file(const char  *path );
int getNextMutFile( DIR *dir_MUT, FILE *Flist_MUT, char *F_path );
int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in );

void tst_loadMut();
void print_N_zone( );

#endif

