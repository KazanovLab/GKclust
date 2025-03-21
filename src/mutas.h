
#ifndef _MUTAS_H__
#define _MUTAS_H__

#include "vector"
using namespace std;

//#define PRINT_HISTO_SEGMUT
#define DEBUG_TRACE_TSLD


//#define XRO_SET_FILENAME "xromo_set.txt"
//#define FOLDER_HUMAN "/Users/gena/mutas/humanXRO/"


//#define XRO_SET_FILENAME "hg19.fa"
//#define FOLDER_HUMAN "/Users/gena/mutas/humanGen/"
// path= /Users/gena/mutas/humanGen/hg19.fa

//#define FOLDER_MUT_DATA "/Users/gena/mutas/mut_set_orig/"
// path=  /Users/gena/mutas/mut_set_orig/cda1a403-16b6-487c-a82a-c377d1d0f89d.consensus.20160830.somatic.snv_mnv.vcf

#define FOLDER_OUT_DATA "/Users/gena/CluSet/Classik/outd/"
#define FOLDER_TEST_OUT "/Users/gena/CluSet/Classik/outd/"

//#define PROFCLUST_FILE_PATH "/Users/gena/CluSet/Clusters/SigProfilerCLusters_cda1a403-16b6-487c-a82a-c377d1d0f89d.txt" // #17
//#define PROFCLUST_MUT 17
//int  LoadXclust( );
void printCompare( );

#define NO_HUMAN_XRO 24
#define NO_CANCER_ID 6

//#define SEGM_SIZE 500000
#define RANDOM_CYCLES 1000

//void getFileName( const char *Canc_File, char *f_name );
int  xtrSamplName ( const char  *vcf_Fpath, char *SamplName, char *Folder=NULL);
int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in );

//const char * getNameMutaF ( int indCancer );

void tst_loadMut();
void print_MutSize();
void print_N_zone( );
void testXsize( );


#endif

// /Users/gena/mutas/mut_set/snv_mnv_BLCA-US.tsv
