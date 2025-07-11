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
		int nucPos;        // from 1
		char nucREF;
		char nucALT[8];
        char APOtag;
        float calcVAF[3];
        int iClust;
        int iAggrCL;

    MUTANT() {nucPos=-1; nucREF='?'; nucALT[0]='?'; APOtag=0; calcVAF[0]=-1;calcVAF[1]=-1;calcVAF[2]=-1;
        iClust = -1; iAggrCL = -1; }; //_cID2 = "NONE";
    MUTANT(int P) {nucPos=P; };
    MUTANT(int P, char R, char *pA, char APO)
        {nucPos=P; nucREF=R; strncpy(nucALT, pA, 8); APOtag=APO; iClust = -1; iAggrCL = -1; }; //_cID2 = "NONE";
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

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

#define XRO_ID_SIZE 64

class  XROSOMA {
private:
//		vector <HUGEN> :: iterator begSearch;
public:
    static long maxXsize;
    static long maxMUTsize;
    static int  chrIDmode;
	string XfName;
	string  XroID;
//	int  chrNum;				// get from Xref_NC
	string version;
//
//  ===== real body of XRO will be read from files =============================
    long Xsize;
    char *Xbody;
    char *Xtag;
//
// ======= SAMPLE part ======================
    vector < MUTANT > vMutAPO;
    int PorogMut;
    
    vector < CLUSTER > vClust;
    vector < CLUSTER > vAggrC;
    int PorogClust;
    
    void addDublRndMut( );
    int calcMutPorog( );
    void  findClusters( int *clustID );
    void ClustConnector( int *clustID );
    void addDublRndClust( );
    int calcClustPorog( ); 
    
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
///////////
    XROSOMA() { XroID="?"; Xsize=0L; Xbody=NULL; };
    XROSOMA(char *xName) { XroID=xName; Xsize=0L; Xbody=NULL; };
    

    int testValidDNK( int Pos, char chREF );
    int APOtest( long Pos, char chREF, char chALT );
    int CLtest_N( pair<long,int> Clust );

    int CLtest_Cross( pair<long,int> Clust, vector < pair<long,int> >::iterator Iter);
    int testMutUnClust( );
    int  getGapClust(MUTANT *pMUt, int next) ;  // +1 = next; -1 = prev
};
/////////////////////////////////////////////////////////

int LoadHuGen( const char *fPath );
int loadVCFMutation( const char *vcf_Fname );
void clearSampl( );
void memRezerv( );
void printMutRezerv( );

void printCluMut( int clustOnly, FILE * fRezult );
void printCluTrace(  );
///////////////////////////////////////////////

int findXroByID( char *xID, int say=1 );
int getNucID(const char Nuc);
int selectXid ( char *buff );       // see it at 'mClust' prog

#endif
