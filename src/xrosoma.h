#ifndef _XROSOMA_H__
#define _XROSOMA_H__


#include "vector"
#include "string"
using namespace std;

#define _C 0
#define _G 1
#define _A 2
#define _T 3

//////////////////////////////////////////////////////////

struct MUSTAT
{
    int Dmin;       // min   distance between mutations
    int Dmax;       // max     --       --       --
    int Hmidl;      // middle  --       --       --    in Histo
    int Hmedi;      // mediana --       --       --       --
    int cluSiz;     // Cluster Size == parametr
    vector < pair <int,int> > Hist_Dist;  // < Distance, cntMut >
    
    MUSTAT() { Dmin=0; Dmax=0; Hmidl=0; Hmedi=0; };
};
//////////////////////////////////////////////////////////

struct MUTANT
{
		long nucPos;        // from 1
		char nucREF;
		char nucALT;
        char APOtag;
        int iClust;
        int iAggrCL;

    MUTANT() {nucPos=-1; nucREF='?'; nucALT='?'; APOtag=0; iClust = -1; iAggrCL = -1; }; //_cID2 = "NONE";
    MUTANT(long P) {nucPos=P; };
    MUTANT(long P, char R, char A, char APO)
        {nucPos=P; nucREF=R; nucALT=A; APOtag=APO; iClust = -1; iAggrCL = -1; }; //_cID2 = "NONE"; 
};
bool lesser_MUT ( const MUTANT &x1, const MUTANT &x2 );

//////////////////////////////////////////////////////////

struct CLUSTER
{
    int iBegMu;          // first index at vMutAPO
    int iEndMu;          // last index at vMutAPO
    int cLng;           //  Cluster Size
    int cID;
    CLUSTER() { iBegMu=-1; iEndMu=-1; cLng=0; cID=0;};
    CLUSTER(int bP, int eP) { iBegMu=bP; iEndMu=eP; cLng=0; cID=0; };
};

//////////////////////////////////////////////////////////

class SAMPLE {
public:
    static long maxMUTsize;
    string SampName;
    vector < MUTANT > vMutAPO[NO_HUMAN_XRO];
    int PorogMut[NO_HUMAN_XRO];   //ClustPorog
    
    vector < CLUSTER > vClust[NO_HUMAN_XRO];
    vector < CLUSTER > vAggrC[NO_HUMAN_XRO];
    int PorogClust[NO_HUMAN_XRO];
    
    SAMPLE(string S) { SampName=S; };
    
    void clearSampl( );
    int calcMutPorog( int Xsiz ); //, float PartOfMut);
    int calcClustPorog( int iXro); //, float PartOfMut);
    void  findClusters( int iXro, int *clustID, int Porog );
    void ClustConnector(  int iXro, int *clustID, int Porog  );
    
    void printCluMut( int clustOnly, FILE * fRezult );
    void printCluTrace(  );
};
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

class  XROSOMA {
private:
//		vector <HUGEN> :: iterator begSearch;
public:
    static long maxXsize;
	string XfName;
	string Xref_NC;
	int  chrNum;				// get from Xref_NC
	string version;
//
//  ===== real body of XRO will be read from files =============================
    long Xsize;
    char *Xbody;
    char *Xtag;   
    
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
/*  ===== Targets (TCx, xGA) for APOBEC mutations  =====================
    long TCxMemSize;    // = sizeof (TCXtag)
    char *TCXtag;       //contains tags from '*Xtag' with TCx_TAG
                        //                      for current segment
    long TCXsize;      // = real no Targets in *TCXtag  (<=TCxMemSize)
*/

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

///////////
    XROSOMA() { chrNum=0; Xsize=0L; Xbody=NULL; };
    XROSOMA(int cN) { chrNum=cN; Xsize=0L; Xbody=NULL; };
    

    int testValidDNK( int Pos, char chREF, char chALT );
    int APOtest( long Pos, char chREF, char chALT );
    int CLtest_N( pair<long,int> Clust );

    int CLtest_Cross( pair<long,int> Clust, vector < pair<long,int> >::iterator Iter);
};
/////////////////////////////////////////////////////////

//int LoadDNKset( const char *pFolder, const char *fSetName );


int LoadHuGen( const char *fPath );
int loadVCFMutation( const char *vcf_Fname );

///////////////////////////////////////////////

int findXroByID( int xID, int say=1 );
int getNucID(const char Nuc);
int selectXid ( char *buff );

#endif
