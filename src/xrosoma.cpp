


#include <string.h>
#include <time.h>
#include <ctype.h>
#include <algorithm>
#include <functional>      // For greater<int>( )

#include "cmain.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
	extern FILE *tsld;
#endif

vector < XROSOMA > vecDNK;

char complment[4][2] = { {'C','G'},  {'G','C'}, {'A','T'}, {'T','A'} };
char NucID[5] = "CGAT";
int  cmpInd[4] = { _G, _C, _T, _A };

///////////////////////////////////////////////////////

int selectXid ( char *buff )
{
    int xID;
    
    if ( strstr(buff,">chr") != buff )
        return 0;
    
    if ( *(buff+4) == 'X')
        xID = 23;
    else
        if ( *(buff+4) == 'Y')
            xID = 24;
        else
            xID = atoi((buff+4));
    
    return xID;
}
///////////////////////////////////////////////////////

int LoadHuGen( const char *fPath )
{
    FILE *fHuGen=NULL;
    char fBuff[1024];
    char *pb;
    int Xro_ID, n;
    string sGen;
    
    printf("\nLoad HuGen  from '%s'.....\n", fPath);
    
    if ( ( fHuGen=fopen(fPath, "r"))==NULL )    {
        printf("File '%s' : ERR_OPENop\n",fPath);
        return -1;
    }
    int cntXro = 0;
    XROSOMA *pCurXro = NULL;
    while ( 1 )    {
        if( fgets(fBuff, sizeof(fBuff)-1, fHuGen) == NULL )
            break;
        if ( fBuff[0] != '>' )    {
            if ( (pb = strchr(fBuff, '\r')) )   //for Windows
                *pb = '\0';
            else
                if ( (pb = strchr(fBuff, '\n')) )
                    *pb = '\0';
            pb = fBuff;
            while ( *pb )   {
                *pb = toupper(*pb);
                pb++;
            }
            sGen += fBuff;
            continue;
        }
// new XRO
        if ( pCurXro )  {
            pCurXro->Xsize = sGen.size() ;
            pCurXro->Xbody =  new char[pCurXro->Xsize+1];//+3];
            if ( pCurXro->Xsize > pCurXro->maxXsize )
                pCurXro->maxXsize = pCurXro->Xsize;
            strcpy(pCurXro->Xbody, sGen.c_str() );
            cntXro++;
            printf("X.%02d  ", pCurXro->chrNum);
        }
        Xro_ID = selectXid ( fBuff );
        if ( ! Xro_ID )
            break;
        if ( findXroByID(Xro_ID, 0) >= 0 )     // patch
            break;
        if ( (vecDNK.size() == 0) || ( vecDNK.back().chrNum < Xro_ID )  ) {
            vecDNK.push_back(XROSOMA(Xro_ID));
            n = (int)vecDNK.size() - 1;
        } else {
            for ( n=0; n<vecDNK.size(); n++ )
                if ( vecDNK[n].chrNum > Xro_ID )
                    break;
            vecDNK.insert(vecDNK.begin()+n, XROSOMA(Xro_ID));
        }
        pCurXro = &vecDNK[n];
        pCurXro->XfName.assign(fPath);
        sGen.clear();
    }
    fclose (fHuGen );
    return cntXro;
}
////////////////////////////////////////////////////////////////////////////////////

int XROSOMA::APOtest( long Pos, char chREF, char chALT )
{
    // APOBEC:
    //    a) TCx ---> 'C' > ['G' | 'T' ]
    //    b) xGA ---> 'G' > ['C' | 'A' ]
    // Returns: 1=APOBEC; 2=APOBEC in ZONE; 0= any other
    //
    
    char *pB = Xbody+Pos-1;
    int retC=0;
    
    if (  *pB != chREF ) {
        printf( "\n !!! NUCLEO_mismatch:: CHR_%d[%ld]='%c';  MUTATION: '%c' > '%c'\n",
               chrNum, Pos, *pB, chREF, chALT);
        return 0;
    }
    
    switch ( chREF ) {
        case 'C':
            if ( ! (chALT == 'G' || chALT == 'T') )
                return 0;
            if ( *(pB-1) != 'T' )
                return 0;
            retC = 1;
            break;
        case 'G':
            if ( ! (chALT == 'C' || chALT == 'A') )
                return 0;
            if ( *(pB+1) != 'A' )
                return 0;
            retC = 1;
            break;
        default:
            return 0;
    }
    
    return retC;
}
////////////////////////////////////////////////////////////////////////////

int XROSOMA::testValidDNK( int Pos, char chREF, char chALT )
{
    //    char *pB = Xbody+Pos-1;
    
    if (  *(Xbody+Pos-1) != chREF ) {
        printf( "\n !!! NUCLEO_mismatch:: CHR_%d[%d]='%c';  MUTATION: '%c' > '%c'\n",
               chrNum,Pos,*(Xbody+Pos-1), chREF, chALT);
        return 0;
    }
    return 1;
}
////////////////////////////////////////////////////////////////////////////////////

int  findXroByID( int xID, int say )
{
    int indX;
    
    for ( indX=0; indX<vecDNK.size(); indX++ )
        if ( vecDNK[indX].chrNum == xID )
            return indX;
    if (say)
        printf("ERR: XroID=%d NOT FOUND\n", xID);
    return -1;
}
 
////////////////////////////////////////////////////// //////////////////////////////

int getNucID ( const char Nuc )
{
    //    char *pp = strchr(NucID, Nuc);
    //    return ( ( !pp ) ? -1 : (int)(pp-NucID) );
    switch (Nuc) {
        case 'C':
            return _C;
        case 'G':
            return _G;
        case 'A':
            return _A;
        case 'T':
            return _T;
        default:
            break;
    }
    return -1;
}
////////////////////////////////////////////////////////////////////////////////////

void print_N_zone( )
{
    char *pBeg, *pEnd;
    char *pBody;
    
#ifdef DEBUG_TRACE_TSLD
    for ( int nX=0; nX < NO_HUMAN_XRO; nX++)    {
        pBody = vecDNK[nX].Xbody;
        fprintf(tsld, "X.%d", nX );
        pBeg = vecDNK[nX].Xbody;
        while ( *pBeg ) {
            while ( *pBeg != 'N' && *pBeg )
                pBeg++;
            if ( !(*pBeg) )
                break;
            pEnd = pBeg;
            while ( *pEnd == 'N' )
                pEnd++;
            fprintf(tsld, "\t%ld : %ld", (pBeg-pBody), (pEnd-pBody) );
            pBeg = pEnd;
        }
        fprintf(tsld, "\n");

    }
#endif
}
////////////////////////////////////////////////////////////////////////////////////

void testXsize( )
{
#ifdef DEBUG_TRACE_TSLD
    fprintf (tsld, "#X\tXsiz\tsLen\n");
    for ( int nX=0; nX < NO_HUMAN_XRO; nX++)
        fprintf (tsld, "%d\t%ld\t%ld\n", vecDNK[nX].chrNum, vecDNK[nX].Xsize, strlen(vecDNK[nX].Xbody) );
#endif
        return;
}
////////////////////////////////////////////////////////////////////////////////////
 
