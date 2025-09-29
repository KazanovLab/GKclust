


#include <string.h>
#include <time.h>
#include <ctype.h>
#include <algorithm>
#include <functional>      // For greater<int>( )

#include "cmain.h"
#include "xrosoma.h"

vector < XROSOMA > vecDNK;

char complment[4][2] = { {'C','G'},  {'G','C'}, {'A','T'}, {'T','A'} };
char NucID[5] = "CGAT";
int  cmpInd[4] = { _G, _C, _T, _A };

#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

///////////////////////////////////////////////////////

void xtrctXID( char *pBuff)
{
    char xID[XRO_ID_SIZE];
    char *pX = xID ;
    char *pB = pBuff ;
    
    pB++;   // >
    while ( *pB && *pB==' ') pB++;
    strncpy(xID, pB, XRO_ID_SIZE-1);
    
    while ( *pX && *pX > ' ') pX++;
    *pX = '\0';
    strcpy(pBuff,  xID);
    
    return;
}
///////////////////////////////////////////////////////

int LoadHuGen( const char *fPath )
{
    FILE *fHuGen=NULL;
    char fBuff[1024];
    char *pb;
//    int n;  //, Xro_ID,
    string sGen;
    
    printf("\nLoading genome from '%s'.....\n", fPath);
    
    if ( ( fHuGen=fopen(fPath, "r"))==NULL )    {
        printf("Failed to open reference genome: '%s'\n",fPath);
        return -1;
    }
    int cntXro = 0;
    XROSOMA *pCurXro = NULL;
    while ( 1 )    {
        if( fgets(fBuff, sizeof(fBuff)-1, fHuGen) == NULL )
            break;
        if ( (pb = strchr(fBuff, '\r')) )   //for Windows
            *pb = '\0';
        else
            if ( (pb = strchr(fBuff, '\n')) )
                *pb = '\0';
        if ( fBuff[0] != '>' )    {
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
            printf("%s  ", pCurXro->XroID.c_str() );
        }
        
        xtrctXID( fBuff );
        if ( findXroByID(fBuff, 0) >= 0 ) {
            printf("Redefinition chroID=%s", fBuff );
            fclose (fHuGen );
            return -1;
        }
        
//        if ( (vecDNK.size() == 0) || ( vecDNK.back().chrNum < Xro_ID )  ) {
//            n = (int)vecDNK.size() - 1;
//        } else {
//            for ( n=0; n<vecDNK.size(); n++ )
//                if ( vecDNK[n].chrNum > Xro_ID )
//                    break;
//            vecDNK.insert(vecDNK.begin()+n, XROSOMA(Xro_ID));
//        }
        
        vecDNK.push_back(XROSOMA(fBuff));
        pCurXro = &vecDNK[vecDNK.size()-1];
        pCurXro->XfName.assign(fPath);
        if ( pCurXro->chrIDmode < 0 )      // first
            pCurXro->chrIDmode =  (strstr(fBuff, "chr") == fBuff) ? 1 : 0;

        sGen.clear();
    }
//
    pCurXro->Xsize = sGen.size() ;
    pCurXro->Xbody =  new char[pCurXro->Xsize+1];//+3];
    if ( pCurXro->Xsize > pCurXro->maxXsize )
        pCurXro->maxXsize = pCurXro->Xsize;
    strcpy(pCurXro->Xbody, sGen.c_str() );
    cntXro++;
    printf("%s  ", pCurXro->XroID.c_str() );
//    pCurXro->xroSAMPL.
    printf("\n=== Loaded %d chro\n", cntXro);
    
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
        printf( "\n !!! NUCLEO_mismatch:: CHR_%s[%ld]='%c';  MUTATION: '%c' > '%c'\n",
               XroID.c_str(), Pos, *pB, chREF, chALT);
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

int XROSOMA::testValidDNK( int Pos, char chREF ) 
{
    //    char *pB = Xbody+Pos-1;
    
    if (  *(Xbody+Pos-1) != chREF ) {
        printf( "\nReference base mismatch at %s:%d. Reference genome has '%c', but VCF REF = '%c'.\n",
               XroID.c_str(), Pos,*(Xbody+Pos-1), chREF);
        return 0;
    }
    return 1;
}
////////////////////////////////////////////////////////////////////////////////////

int findXro( char *xID )    {
    for ( int indX=0; indX<vecDNK.size(); indX++ )  {
        if ( vecDNK[indX].XroID == xID )
            return indX;
    }
    return -1;
}
////////////////////////////////////////////////////////////////////////////////////

int  findXroByID( char *xID, int say )
{
    int indX;
    int stat;
    char my_xID[XRO_ID_SIZE];
    
    if  ( strstr(xID,"chr")==xID )
        stat = ( vecDNK[0].chrIDmode ) ? 1 : 2;
    else
        stat = ( vecDNK[0].chrIDmode ) ? 3 : 4;
//     vcf     genome
// 1=(chrID) : (chrID)
// 2=(chrID) : (   ID)
// 3=(   ID) : (chrID)
// 1=(   ID) : (   ID)
    
    strncpy(my_xID, xID, XRO_ID_SIZE-1);
    switch (stat) {
            
        case 2:         // (chrID) : (   ID)
            if ( (indX = findXro((xID+3)) ) >= 0 )    // (chrID->ID) : ( ID )
                break;
            
            if ( strcmp(xID, "chrM")==0 )           // (chrM->MT) : ( ID )
                strcpy(my_xID, "MT");
            if ( (indX = findXro(my_xID)) >= 0 )
                break;
            
            if ( (indX = findXro(xID)) >= 0 )      // (chrID) : ( ID )   for finding in patches
                break;
            
            indX =-1;
            break;
            
        case 3:         // (   ID) : (chrID)
            strcpy(my_xID,"chr");
            strcat(my_xID, xID);
            if ( (indX = findXro(my_xID)) >= 0 )      // (ID->chrID) : (chrID)
                break;
            
            if ( strcmp(xID, "MT")==0 )
                strcpy(my_xID, "chrM");
            if ( (indX = findXro(my_xID)) >= 0 )      // (MT->chrM) : ( chrID )
                break;
            
            if ( (indX = findXro(xID)) >= 0 )      // (ID) : (chrID)   for finding in patches
                break;
            
            indX =-1;
            break;

        default:        // case1:: (chrID) : (chrID)
                        // case4:: (   ID) : (   ID)
            indX = findXro(xID);
            break;

    }
    
    if (say && indX < 0 )
        printf("ERR: XroID='%s' NOT FOUND\n", xID);
    
    return indX;
}
////////////////////////////////////////////////////////////////////////////////////

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
    for ( int nX=0; nX < vecDNK.size()-1; nX++)    {
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
 
