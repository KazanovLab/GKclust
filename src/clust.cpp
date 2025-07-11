//
//  clust.cpp
//  APOclust
//
//  Created by Gena on 2/26/24.
//  Copyright © 2024 Gena. All rights reserved.
//

#include <stdlib.h>
#include "cmain.h"
#include "xrosoma.h"

extern PROGARGS ArgKit;

extern vector < XROSOMA > vecDNK;
extern string SamplName;

vector < long > vRndMut;               // posMUT for Random simulator
vector < pair <long, int> > vRndClu;   // <begPos, cLng>

vector < int > MutSpace;  // Space between Mutations

bool lesser_pair ( const pair<long,int> &x1, const pair<long,int> &x2 ) {
    return x1.first < x2.first;
}

#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

//////////////////////////////////////////////////////////////////////////////////

void memRezerv( )
{
    
    if ( vRndMut.capacity() < XROSOMA::maxMUTsize )
        vRndMut.reserve(XROSOMA::maxMUTsize);
    if ( MutSpace.capacity() < RANDOM_CYCLES * XROSOMA::maxMUTsize )
        MutSpace.reserve(RANDOM_CYCLES * XROSOMA::maxMUTsize );
    
    return;
}
//////////////////////////////////////////////////////////////////////////////////

void printMutRezerv( )
{
    printf ( "\nvRndMut.capacity=%ld\tMutSpace.reserve=%ld\nmaxMUTsize=%ld\tnmaxMUTsize*nC=%ld\n",
            vRndMut.capacity(), MutSpace.capacity(), XROSOMA::maxMUTsize, RANDOM_CYCLES * XROSOMA::maxMUTsize);
    return;
}
//////////////////////////////////////////////////////////////////////////////////

void XROSOMA:: findClusters( int *clustID )      //iXro , int Porog
{
    int cntM;
    int MaxP = (int)(vMutAPO.size() -1);
    CLUSTER Clust;
    
    int nM=0;
    while  ( nM <= MaxP ) {
        cntM = 0;
        for ( int n=nM+1; n <= MaxP; n++,cntM++ ) {
            if ( (vMutAPO[n].nucPos -vMutAPO[n-1].nucPos) > PorogMut )
                break;
        }
        if ( cntM > 0 ) {
            Clust.iBegMu = nM;
            Clust.iEndMu = nM + cntM;
            Clust.cLng = (int)(vMutAPO[Clust.iEndMu].nucPos - vMutAPO[Clust.iBegMu].nucPos);
            Clust.cID   = *clustID;
            vClust.push_back(Clust);
            for ( int n=0; n <= cntM; n++ ) {
                vMutAPO[nM+n].iClust = (int)vClust.size() - 1 ;
            }
        }

        if ( cntM > 0 )
            (*clustID)++ ;
        nM += (cntM+1);
    }
    return;
    
}
///////////////////////////////////////////////////////////////////////////////////

void XROSOMA:: ClustConnector( int *clustID )    // int iXro, int Porog
{
    int cnt;
    int MaxCL = (int)(vClust.size() -1);
    CLUSTER CombinedCL;
    int headMut, tailMut_1;
    
    int nCL=0;
    while  ( nCL <= MaxCL ) {
        cnt = 0;
        for ( int n=nCL+1; n <= MaxCL; n++,cnt++ ) {
            headMut = vClust[n].iBegMu;
            tailMut_1 = vClust[n-1].iEndMu;
            if ( (vMutAPO[headMut].nucPos - vMutAPO[tailMut_1].nucPos) > PorogClust )
                break;
        }
        if ( cnt > 0 ) {
            CombinedCL.iBegMu = vClust[nCL].iBegMu;
            CombinedCL.iEndMu = vClust[nCL+cnt].iEndMu;
            CombinedCL.cLng = (int)(vMutAPO[CombinedCL.iEndMu].nucPos -
                                    vMutAPO[CombinedCL.iBegMu].nucPos);
            CombinedCL.cID   = *clustID;
            vAggrC.push_back(CombinedCL);
            for ( int n=CombinedCL.iBegMu; n <= CombinedCL.iEndMu; n++ )    {
                vMutAPO[n].iAggrCL = (int)vAggrC.size() - 1 ;
            }
        }

        if ( cnt > 0 )
            (*clustID)++ ;
        
        nCL += (cnt+1);
    }
    return;

}
///////////////////////////////////////////////////////////////////////////////////


void XROSOMA:: addDublRndMut( )
{
    int nM;
    int Msiz = (int)vMutAPO.size();
    unsigned long Rnd_pos;
    vector <long>::iterator Iter;
    
    nM = 1;
    while (nM < vRndMut.size() ) {
        if ( vRndMut[nM] == vRndMut[nM-1] )   {
            vRndMut.erase(vRndMut.begin()+nM);
            continue;
        }
        nM++;
    }
    if ( vRndMut.size() == vMutAPO.size() )
        return;
    
    while ( vRndMut.size() < vMutAPO.size() )   {
        Rnd_pos = ( (double)rand() / (double)RAND_MAX ) * Xsize;
        if ( *(Xbody+Rnd_pos)=='N' )
            continue;
        if ( vRndMut.back() < Rnd_pos )   {      // ? <=  | < ?
            vRndMut.push_back( Rnd_pos );
            continue;
        }
        Iter = lower_bound( vRndMut.begin( ), vRndMut.end( ), Rnd_pos );   // MUTANT(NucPos),  lesser_MUT);
        if ( *Iter == Rnd_pos )
            continue;
        vRndMut.insert( Iter, Rnd_pos );
    }
    
    return;
}
//////////////////////////////////////////////////////////////////////////////////

int XROSOMA:: calcMutPorog(  )
{
//    int Xsiz = (int)vecDNK[iXro].Xsize;
    int Msiz = (int)vMutAPO.size();

    unsigned int nM, L;
    unsigned long Rnd_pos;
    vector <long>::iterator Iter;
    
    if ( Msiz < 3 )
        return 0;
 
    clock_t startP = clock();
    
    MutSpace.clear();
    
//    int cntN=0, cntDubl=0;
    for ( int nCycl=0; nCycl<RANDOM_CYCLES; nCycl++ )    {
        nM = 0;
        vRndMut.clear();
        while (nM < Msiz) {
            Rnd_pos = ( (double)rand() / (double)RAND_MAX ) * Xsize;
            if ( *(Xbody+Rnd_pos)=='N' )
                continue;

            vRndMut.push_back(Rnd_pos);
            nM++;
        }

        sort(vRndMut.begin(), vRndMut.end());
        
        addDublRndMut( );
        
        for ( nM=1; nM<Msiz; nM++ ) {
            L = (int)(vRndMut[nM] - vRndMut[nM-1]);
            MutSpace.push_back(L);
        }
        
    }    // of limon cycles
    
    sort(MutSpace.begin(), MutSpace.end());
    
    float fP;
    sscanf(ArgKit.chPart, "%f", &fP);
    
    int indPor = fP * RANDOM_CYCLES  * (Msiz-1);  //(fP/100)
    PorogMut  = MutSpace[indPor];
    
    if ( GET_ARG_D(ArgKit.argTAG) ) {
        char chP[16];
        sprintf(chP, "%d", PorogMut);
        int step = (int)(MutSpace.size() / 1000);
        nM = 0;
        fprintf(ArgKit.foutRndMu, "Chr\tChrSize\tnMut\tSpace\tGap\n");
        while ( nM < MutSpace.size() ) {
            fprintf(ArgKit.foutRndMu,"%s\t%ld\t%d\t%d\t%s\n", XroID.c_str(), Xsize,
                    Msiz, MutSpace[nM], ( ( abs((int)nM-indPor) <= step ) ? chP : " " ) );
            nM += step;
        }

    }
    clock_t stopP = clock();
    double duration = (double)(stopP - startP) / CLOCKS_PER_SEC;
//    printf("%s.calcMutPorog [%d] dT=%5.2f\n", XroID.c_str(), Msiz, duration );
    
    return PorogMut;
}
//////////////////////////////////////////////////////////////////////////////////

void XROSOMA:: addDublRndClust( )
{
    int nCL;
    unsigned long Rnd_pos;
    pair < long, int > Rnd_CL;
    vector < pair<long,int> >::iterator Iter;
    
    vector < pair<long,int> > vRw;
    vRw = vRndClu;
    vRndClu.clear();
    
    vRndClu.push_back(vRw[0]);
    nCL= 1;
    int vCLu_siz = (int)vClust.size();
    
    while (nCL < vCLu_siz ) {
        if ( vRw[nCL].first == vRndClu.back().first )   {
            vRw[nCL].first = -1;
            nCL++;
            continue;
        }
        if ( CLtest_Cross( vRw[nCL],  vRndClu.end() ) > 0 )    {
            vRw[nCL].first = -1;
            nCL++;
            continue;
        }
        vRndClu.push_back(vRw[nCL]);
        nCL++;
    }
    
//    if ( vRndClu.size() == vClust.size() )
//        return;
    
    nCL= 0;
    while ( vRndClu.size() < vClust.size() ) {
        if ( vRw[nCL].first >= 0 )  {
            nCL++;
            continue;
        }
        Rnd_pos = ( (double)rand() / (double)RAND_MAX ) * Xsize;
        if ( ( Rnd_pos + vRw[nCL].second ) >= Xsize )
            continue;   // try again at nCL
        Rnd_CL.first  = Rnd_pos;
        Rnd_CL.second = vRw[nCL].second;
        if ( CLtest_N( Rnd_CL ) > 0 )
            continue;    // try again at nCL
        
        Iter = lower_bound( vRndClu.begin( ), vRndClu.end( ), Rnd_CL, lesser_pair);
        if ( CLtest_Cross(Rnd_CL,  Iter) > 0 )
            continue;       // try again at nCL
        
        vRndClu.insert( Iter, Rnd_CL );
        nCL++;
    }
    
    return;
}
/////////////////////////////////////////////////////////////////////////

int XROSOMA:: calcClustPorog(  ) //int iXro )
{
//    int Xsiz = (int)vecDNK[iXro].Xsize;
    int nCL;
    unsigned long Rnd_pos;
    pair < long, int > Rnd_CL;
//    vector < pair<long,int> >::iterator Iter;

    CLUSTER *pCLU;
    
    int vCLu_siz = (int)vClust.size();
    
    if ( vCLu_siz < 3 )
        return 0;
    
    clock_t startP = clock();
    
    MutSpace.clear();
    
//    int cntN=0, cntDubl=0;
    for ( int nCycl=0; nCycl<RANDOM_CYCLES; nCycl++ )    {
        nCL = 0;
        vRndClu.clear();
        while (nCL < vCLu_siz) {
            Rnd_pos = ( (double)rand() / (double)RAND_MAX ) * Xsize;
            pCLU = &vClust[nCL];
            if ( (Rnd_pos + pCLU->cLng) >= Xsize )
                continue;
            Rnd_CL.first  = Rnd_pos;
            Rnd_CL.second = pCLU->cLng;
            if ( CLtest_N( Rnd_CL ) > 0 )
                continue;
            
            vRndClu.push_back( Rnd_CL );
            nCL++;
        }
        
        sort(vRndClu.begin(), vRndClu.end(), lesser_pair );
        
        addDublRndClust( );
        
        int L;
        for ( nCL=1; nCL<vCLu_siz; nCL++ ) {
            L = (int)(vRndClu[nCL].first - (vRndClu[nCL-1].first + vRndClu[nCL-1].second) );
            if ( L <= 0 )   {
                printf("ERR Space=%d: RND_CL_#%d[%ld, %d] :: CL_#%d[%ld, %d]\n", L,
                       nCL, vRndClu[nCL].first, vRndClu[nCL].second,
                       nCL-1, vRndClu[nCL-1].first, vRndClu[nCL-1].second);
            }
            MutSpace.push_back(L);
        }
        
    }    // of limon cycles
            
    sort( MutSpace.begin(), MutSpace.end() );

    float fP;
    sscanf(ArgKit.chPart, "%f", &fP);
    
    int indPor = fP * RANDOM_CYCLES  * (vCLu_siz-1); //(fP/100)
    PorogClust  = MutSpace[indPor];
    
    if ( GET_ARG_D(ArgKit.argTAG) ) {
        char chP[16];
        sprintf(chP, "%d", PorogClust);
        int step = (int)(MutSpace.size() / 1000);
        nCL = 0;
        fprintf(ArgKit.foutRndCl, "Chr\tChrSize\tnMut\tSpace\tGap\n");
        while ( nCL < MutSpace.size() ) {
            fprintf(ArgKit.foutRndCl,"%s\t%ld\t%d\t%d\t%s\n", XroID.c_str(), Xsize,
                    vCLu_siz, MutSpace[nCL], ( ( abs((int)nCL-indPor) <= step ) ? chP : " " ) );
            nCL += step;
        }

    }

    clock_t stopP = clock();
    double duration = (double)(stopP - startP) / CLOCKS_PER_SEC;
//    printf("%s.calcClustPorog [%d] dT=%5.2f\n\n", XroID.c_str(), vCLu_siz, duration );
    
    return PorogClust;
    
}

/////////////////////////////////////////////////////////////////////////
/*
int XROSOMA:: calcClustPorog(  ) //  insert variant
{
//    int Xsiz = (int)vecDNK[iXro].Xsize;
    int nCL;
    unsigned long Rnd_pos;
    pair < long, int > Rnd_CL;
    vector < pair<long,int> >::iterator Iter;
    CLUSTER *pCLU;
    
    int vCLu_siz = (int)vClust.size();
    
    if ( vCLu_siz < 3 )
        return 0;
    
    clock_t startP = clock();
    
    MutSpace.clear();
    
//    int cntN=0, cntDubl=0;
    for ( int nCycl=0; nCycl<RANDOM_CYCLES; nCycl++ )    {
        nCL = 0;
        vRndClu.clear();
        while (nCL < vCLu_siz ) {
            Rnd_pos = ( (double)rand() / (double)RAND_MAX ) * Xsize;
            pCLU = &vClust[nCL];
            if ( (Rnd_pos + pCLU->cLng) >= Xsize )
                continue;
            Rnd_CL.first  = Rnd_pos;
            Rnd_CL.second = pCLU->cLng;
            if ( CLtest_N( Rnd_CL ) > 0 )
                continue;
            
            if ( vRndClu.empty() )  {       //}|| vRndMut.back() < Rnd_pos )   {
                vRndClu.push_back( Rnd_CL );
                nCL++;
                continue;
            }
            
            Iter = lower_bound( vRndClu.begin( ), vRndClu.end( ), Rnd_CL, lesser_pair);
            if ( CLtest_Cross(Rnd_CL,  Iter) > 0 )
                continue;
            
            vRndClu.insert( Iter, Rnd_CL );
            nCL++;
        }
        
        
        int L;
        for ( nCL=1; nCL<vCLu_siz; nCL++ ) {
            L = (int)(vRndClu[nCL].first - (vRndClu[nCL-1].first + vRndClu[nCL-1].second) );
            if ( L <= 0 )   {
                printf("ERR Space=%d: RND_CL_#%d[%ld, %d] :: CL_#%d[%ld, %d]\n", L,
                       nCL, vRndClu[nCL].first, vRndClu[nCL].second,
                       nCL-1, vRndClu[nCL-1].first, vRndClu[nCL-1].second);
            }
            MutSpace.push_back(L);
        }
        
    }    // of limon cycles
    
    sort( MutSpace.begin(), MutSpace.end() );
    
    float fP;
    sscanf(ArgKit.chPart, "%f", &fP);
    
    int indPor = fP * RANDOM_CYCLES  * (vCLu_siz-1); //(fP/100)
    PorogClust  = MutSpace[indPor];
    
    if ( GET_ARG_D(ArgKit.argTAG) ) {
        char chP[16];
        sprintf(chP, "%d", PorogClust);
        int step = (int)(MutSpace.size() / 1000);
        nCL = 0;
        fprintf(ArgKit.foutRndCl, "Chr\tChrSize\tnMut\tSpace\tGap\n");
        while ( nCL < MutSpace.size() ) {
            fprintf(ArgKit.foutRndCl,"%s\t%ld\t%d\t%d\t%s\n", XroID.c_str(), Xsize,
                    vCLu_siz, MutSpace[nCL], ( ( abs((int)nCL-indPor) <= step ) ? chP : " " ) );
            nCL += step;
        }
        
    }
    
    clock_t stopP = clock();
    double duration = (double)(stopP - startP) / CLOCKS_PER_SEC;
    printf("%s.calcClustPorog [%d] dT=%5.2f\n\n", XroID.c_str(), vCLu_siz, duration );
    
    return PorogClust;
    
}
*/
//////////////////////////////////////////////////////////////////////////////////

int XROSOMA :: CLtest_N( pair<long,int> Clust) //int xPos, int Lng)
{
    for ( int n=0; n<=Clust.second ; n++) {
        if ( *(Xbody+Clust.first+n)=='N' )
            return 1;
    }
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////

int XROSOMA :: CLtest_Cross(pair<long,int> newClust, vector <pair <long,int> >::iterator ItCLset )
{

    int stat=0;

    if ( ItCLset == vRndClu.begin( ) )
        stat = 1;
    if ( ItCLset == vRndClu.end( ) )
        stat = 2;
    
    switch (stat) {
        case 1:     // insert to head
            if ( (newClust.first + newClust.second) >= ItCLset->first )     //cross with next
                return 1;
            break;
        case 2: // insert to tail
            if ( ( (ItCLset-1)->first + (ItCLset-1)->second) >= newClust.first)     //cross with prev
                return 1;
            break;
        default:
            if ( (newClust.first + newClust.second) >= ItCLset->first )
                return 1;
            if ( ( (ItCLset-1)->first + (ItCLset-1)->second ) >= newClust.first)
                return 1;
            break;
    }
  
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////
//===============================================================================

void printCluTrace(  )
{
    char cName_1[16], cName_2[16];
    char cMutN_1[8], cMutN_2[8];
    char cPorog[8], cDprio[8], cDnext[8];
    int stat_1=0, stat_2=0;   // =1 if NOT Clust;  =2 if itis Clust & 1_st Mut ;  =3 Clust & next Mut
    CLUSTER *pCLU, *pAggr;
    XROSOMA *pXro;

    fprintf(ArgKit.foutTrace ,"Chr\tPos\tRef\tAlt\t"); //APO\t");
    fprintf(ArgKit.foutTrace,"ClustID_1\tMutCnt_1\tGap_Mut\tSpace_up\tSpace_down\t");
    fprintf(ArgKit.foutTrace,"ClustID_2\tMutCnt_2\tGap_Clust\tSpace_up\tSpace_down\n");
    
    for ( int nX=0; nX < vecDNK.size(); nX++ )    {
        pXro = &vecDNK[nX];
        if ( pXro->vMutAPO.size()==0 )
            continue;
        for ( int nMu=0; nMu < pXro->vMutAPO.size(); nMu++ )   {
            int indCl = pXro->vMutAPO[nMu].iClust;
            pCLU =  (indCl < 0) ? NULL : &(pXro->vClust[indCl]);
            indCl = pXro->vMutAPO[nMu].iAggrCL;
            pAggr = (indCl < 0) ? NULL : &(pXro->vAggrC[indCl]);
            strcpy(cMutN_1, " "); strcpy(cMutN_2, " ");
            strcpy(cName_1, " "); strcpy(cName_2, " ");
            //---------
            if ( pXro->vMutAPO[nMu].iClust < 0 )
                stat_1 =1;
            else
                stat_1 = ( nMu==0 || (pXro->vMutAPO[nMu].iClust != pXro->vMutAPO[nMu-1].iClust) ) ? 2 : 3;
            
            if ( pXro->vMutAPO[nMu].iAggrCL < 0 )
                stat_2 =1;
            else
                stat_2 = ( nMu==0 || (pXro->vMutAPO[nMu].iAggrCL != pXro->vMutAPO[nMu-1].iAggrCL) ) ? 2 : 3;
            
            //fprintf(fMutClu,"Xr\tPOS\tREF\tALT\t  //APO
            fprintf(ArgKit.foutTrace,"%s\t%d\t%c\t%s\t",    //%s\t",
                    pXro->XroID.c_str(), (int)pXro->vMutAPO[nMu].nucPos,
                    pXro->vMutAPO[nMu].nucREF, pXro->vMutAPO[nMu].nucALT );
            //                    ( (vMutAPO[nX][nMu].APOtag == 0 ) ? " " : "APO" ) );
            
//----- my_Clust
            if ( stat_1 != 1 )
                sprintf(cName_1, "cID_%d", pCLU->cID );
            if ( stat_1 == 2 )
                sprintf(cMutN_1, "%d", (pCLU->iEndMu - pCLU->iBegMu) + 1 );
            
            if ( stat_1 == 2  )  //
                sprintf(cPorog, "%d", pXro->PorogMut);
            else
                strcpy(cPorog, " ");
            
            if ( stat_1 != 1 )  { //}|| statX != 1 ) {
                sprintf(cDprio, "%d", (  (nMu==0)
                                        ? pXro->vMutAPO[nMu].nucPos
                                        : pXro->vMutAPO[nMu].nucPos - pXro->vMutAPO[nMu-1].nucPos) );
                sprintf(cDnext, "%d", ( (nMu == pXro->vMutAPO.size()-1)
                                        ? (int)pXro->Xsize - pXro->vMutAPO[nMu].nucPos
                                        : pXro->vMutAPO[nMu+1].nucPos - pXro->vMutAPO[nMu].nucPos) );
            } else {
                strcpy(cDprio, " ");
                strcpy(cDnext, " ");
            }
            //                     "cID_1\tnMu_1\tPorog_1\tGAP_up\tGAP_down\t");
            fprintf(ArgKit.foutTrace,"%s\t%s\t%s\t%s\t%s\t", cName_1, cMutN_1, cPorog, cDprio, cDnext );
            
//----- compound_Clust
            if ( stat_2 != 1 )
                sprintf(cName_2, "cID_%d", pAggr->cID );
            if ( stat_2 == 2 )
                sprintf(cMutN_2, "%d", (pAggr->iEndMu - pAggr->iBegMu) + 1 );
            
            if ( stat_2 == 2 )
                sprintf(cPorog, "%d", pXro->PorogClust);
            else
                strcpy(cPorog, " ");
            
            if ( stat_2 != 1 ) {
                sprintf(cDprio, "%d", (  (nMu==0)
                                        ? pXro->vMutAPO[nMu].nucPos
                                        : pXro->vMutAPO[nMu].nucPos - pXro->vMutAPO[nMu-1].nucPos) );
                sprintf(cDnext, "%d", ( (nMu == pXro->vMutAPO.size()-1)
                                        ? (int)pXro->Xsize - pXro->vMutAPO[nMu].nucPos
                                        : pXro->vMutAPO[nMu+1].nucPos - pXro->vMutAPO[nMu].nucPos) );
            } else {
                strcpy(cDprio, " ");
                strcpy(cDnext, " ");
            }
            
            //                       "cID_2\tnMu_2\tPorog_2\tGAP_up\tGAP_down\n");
            fprintf(ArgKit.foutTrace,"%s\t%s\t%s\t%s\t%s\n", cName_2, cMutN_2, cPorog, cDprio, cDnext );
            
        }
    }

    return;
    
}
//===============================================================================

void printCluMut( int clustOnly, FILE * fRezult )
{

    char clustID[16];
    char cntMut[8];
    string valVAF;
    char cValV[128];
    int nV;
    int stat_1=0, stat_2=0;   // =1 if NOT Clust;  =2 if itis Clust & 1_st Mut ;  =3 Clust & next Mut
    CLUSTER *pCLU, *pAggr;
    XROSOMA *pXro;

    fprintf(fRezult, "Chr\tPos\tRef\tAlt\tVAF\t"); //APO\t");
    fprintf(fRezult, "ClusterID\tMutCount\n");
    
    
    for ( int nX=0; nX < vecDNK.size(); nX++ )    {
        pXro = &vecDNK[nX];
        if ( pXro->vMutAPO.size()==0 )
            continue;
        for ( int nMu=0; nMu < pXro->vMutAPO.size(); nMu++ )   {
            if ( pXro->vMutAPO[nMu].iClust < 0 )
                stat_1 = 1;
            else
                stat_1 = ( nMu==0 || (pXro->vMutAPO[nMu].iClust != pXro->vMutAPO[nMu-1].iClust) ) ? 2 : 3;
            
            if ( pXro->vMutAPO[nMu].iAggrCL < 0 )
                stat_2 = 1;
            else
                stat_2 = ( nMu==0 || (pXro->vMutAPO[nMu].iAggrCL != pXro->vMutAPO[nMu-1].iAggrCL) ) ? 2 : 3;
            
            valVAF.clear();
            nV=0;
            while ( 1 )   {
                sprintf(cValV, "%.2f", pXro->vMutAPO[nMu].calcVAF[nV]);
                valVAF += cValV;
                if ( nV==2 || pXro->vMutAPO[nMu].calcVAF[nV+1] < 0 )
                    break;
                valVAF += ", ";
                nV++;
            }
    
            if (stat_1==1 && stat_2==1) {
                if ( ! clustOnly )
                    fprintf(fRezult,"%s\t%d\t%c\t%s\t%s\n", // %s\t",
                            pXro->XroID.c_str(), (int)pXro->vMutAPO[nMu].nucPos,
                            pXro->vMutAPO[nMu].nucREF, pXro->vMutAPO[nMu].nucALT, valVAF.c_str() );
                continue;
            }
            int indCL = pXro->vMutAPO[nMu].iClust;
            pCLU =  (indCL < 0)  ? NULL : &(pXro->vClust[indCL]);
            indCL = pXro->vMutAPO[nMu].iAggrCL;
            pAggr = (indCL < 0) ? NULL : &(pXro->vAggrC[indCL]);
//---------
            //fprintf(fMutClu,"#S\tXr\tPOS\tREF\tALT\t  //APO
            fprintf(fRezult,"%s\t%d\t%c\t%s\t%s\t", // %s\t",
                    pXro->XroID.c_str(), (int)pXro->vMutAPO[nMu].nucPos,
                    pXro->vMutAPO[nMu].nucREF, pXro->vMutAPO[nMu].nucALT, valVAF.c_str() );
            //                    ( (vMutAPO[nX][nMu].APOtag == 0 ) ? " " : "APO" ) );
            
            strcpy(cntMut,  " ");
            strcpy(clustID, " ");
            if ( stat_2 == 1 )  {   // it's not compound cluster
                sprintf(clustID, "cID_%d", pCLU->cID );
                //              if ( stat_1 == 2 )      // print nMut only for 1_st Mut
                sprintf(cntMut, "%d", (pCLU->iEndMu - pCLU->iBegMu) + 1 );
            }
            else {                  // Compound cluster
                sprintf(clustID, "cID_%d", pAggr->cID );
                //                if ( stat_2 == 2 )      // print nMut only for 1_st Mut
                sprintf(cntMut, "%d", (pAggr->iEndMu - pAggr->iBegMu) + 1 );
            }
            //fprintf(fMutClu,"cluID\tnMut\tsampl\n");
            fprintf(fRezult,"%s\t%s\n", clustID, cntMut );
            
        }
    }
    return;
}

//===============================================================================
/*
curCL   nextCL  curAGR  nextAGR     stat
 -      -           -       -       1       singlMUT
 -      -           -       +     x               НЕ кластерная агрегация
 -      -           +       -     x
 -      -           +       +     x
 
 -      +           -       -       2       Claster start .......
 -      +           -       +       3       Clust & Aggr start .......... dL > PorogMut
 -      +           +       -     x
 -      +           +       +     x
 
 +      -           -       -       4       ..... Cluster STOP
 +      -           -       +     x
 +      -           +       -       5       ..... Clust & Aggr STOP
 +      -           +       +     x
 
 +      +           -       -       6       if ( == ) {in sameCluster}  else {test_PorogMut}
 +      +           -       +       7       if ( == ) {ERROR}  else Aggr start {dL > PorogMut;}
 +      +           +       -       8       if ( == ) {ERROR}  else Aggr STOP  {dL > PorogMut;}
 +      +           +       +       9

*/
//////////////////////////////////////////////////////////////////////////////////

int XROSOMA:: testMutUnClust( )
{
    int nM, stat;
    char desiz[8];
    int dLnuc, dLclu;
    MUTANT *pcurnMut, *pnextMut;
    struct  {   //Tdecision
        char CASE[8];
        int state;
    } td[16] = {    {"0000",1}, {"0001",-11}, {"0010",-12}, {"0011",13},
                    {"0100",2}, {"0101",3}, {"0110",-22}, {"0111",23},
                    {"1000",4}, {"1001",-31}, {"1010",5}, {"1011",33},
                    {"1100",6}, {"1101",7}, {"1110",8}, {"1111",9}
                };
    int errCnt=0;
    
    if ( vMutAPO.size() == 0 )
        return 0;
    if ( (this - &vecDNK[0]) == 0 )
        fprintf (tsld, "\n ClusterTEST for '%s'\n", SamplName.c_str() );
    fprintf (tsld, "%s PorogMUT=%d PorogCLust=%d\n", XroID.c_str(), PorogMut, PorogClust);
    memset (desiz, '\0', 8);
    for ( nM=1; nM<vMutAPO.size(); nM++ )   {
        pcurnMut = &vMutAPO[nM-1];
        pnextMut = &vMutAPO[nM];
        dLnuc = pnextMut->nucPos - pcurnMut->nucPos;
        
        desiz[0] = ( pcurnMut->iClust >= 0) ? '1' : '0';
        desiz[1] = ( pnextMut->iClust >= 0) ? '1' : '0';
        desiz[2] = ( pcurnMut->iAggrCL >= 0) ? '1' : '0';
        desiz[3] = ( pnextMut->iAggrCL >= 0) ? '1' : '0';
        stat=0;
        for ( int n=0; n<16; n++ )
            if ( strcmp(desiz, td[n].CASE)==0 ) {
                stat = td[n].state;
                break;
            }
        switch (stat) {
            case 1:     // singl MUT
                if ( dLnuc > PorogMut )
                    continue;
                fprintf (tsld, "ER.%d '%s': Д/Б в кластере [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                break;
                
            case 2:     // Singl ->Claster. Claster start .......
                if ( dLnuc > PorogMut )
                    continue;
                fprintf (tsld, "ER.%d '%s': Не включен в кластер с головы [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                break;
                
            case 3:     // Singl ->(Claster & Aggr). Clust & Aggr start ..........
                if ( dLnuc > PorogMut )
                    continue;
                fprintf (tsld, "ER.%d '%s': Не включен в класт+агр с головы [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                break;
                
            case 4:     // Claster -> Singl ..... Cluster STOP
                if ( dLnuc > PorogMut )
                    continue;
                fprintf (tsld, "ER.%d '%s': Не включен в кластер с хвоста [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                break;

            case 5:     // ..... Clust & Aggr STOP
                if ( dLnuc > PorogMut )
                    continue;
                fprintf (tsld, "ER.%d '%s': Не включен в класт+агр с хвоста [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                break;
                
            case 6:     // "1100" in Clust not in AGGR
                if ( pcurnMut->iClust == pnextMut->iClust )  {   // same Clust
                    if ( dLnuc <= PorogMut)
                        continue;
                    fprintf (tsld, "ER.%d '%s': В кластере далекие мутации [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                }
                if ( pcurnMut->iClust != pnextMut->iClust )  {   // different Clust
                    if ( dLnuc > PorogMut)
                        continue;
                    fprintf (tsld, "ER.%d '%s': Близкие мутации в разных кластерах [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                }
                break;
                    
            case 7:     // if ( == ) {ERROR}  else Aggr start {dL > PorogMut;}
                if ( pcurnMut->iClust != pnextMut->iClust )  {   // different Clust: new AGGR start ...
                    if ( dLnuc > PorogClust )
                        continue;
                    fprintf (tsld, "ER.%d '%s': Не включен в агрегат с головы [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                    break;
                }
                if ( pcurnMut->iClust == pnextMut->iClust )     // same Clust
                    fprintf (tsld, "ER.%d '%s': Агрегат начинается внутри кластера [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                break;
                
            case 8:     // if ( == ) {ERROR}  else Aggr STOP  {dL > PorogMut;}
                if ( pcurnMut->iClust != pnextMut->iClust )  {   // different Clust: ....AGGR stop
                    if ( dLnuc > PorogClust )
                        continue;
                    fprintf (tsld, "ER.%d '%s': Не включен в агрегат с хвоста [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                    break;
                }
                if ( pcurnMut->iClust == pnextMut->iClust )     // same Clust
                    fprintf (tsld, "ER.%d '%s': Агрегат заканчивается внутри кластера [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                break;
                
            case 9:     // in Clust && in AGGR
                if ( pcurnMut->iAggrCL == pnextMut->iAggrCL )
                    continue;
                if ( pcurnMut->iClust != pnextMut->iClust )
                    continue;
                fprintf (tsld, "ER.%d '%s': Кластер в разных агрегатах [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                break;
                
            case 13:     // "0011"
                if ( pcurnMut->iAggrCL == pnextMut->iAggrCL )
                    continue;
                fprintf (tsld, "ER.%d '%s': Некластерная агрегация [M%d]->[A%d]\n", stat, desiz,
                         nM-1, pcurnMut->iAggrCL);
                fprintf (tsld, "ER.%d '%s': Некластерная агрегация [M%d]->[A%d]\n", stat, desiz,
                         nM, pnextMut->iAggrCL);
                break;
                
            case 23:     // "0111" Clust -> Mut  at AGGR
                if ( pcurnMut->iAggrCL != pnextMut->iAggrCL )   {
                    fprintf (tsld, "ER.%d '%s': Некластерная агрегация [M%d]->[A%d]\n", stat, desiz,
                             nM-1, pcurnMut->iAggrCL);
                    break;
                }
                dLclu = getGapClust(pnextMut, -1);  //Prev
                if ( dLclu <= PorogClust )
                    continue;
                fprintf (tsld, "ER.%d '%s': Некластерная агрегация [M%d]->[A%d]\n", stat, desiz,
                         nM-1, pcurnMut->iAggrCL);
                fprintf (tsld, "\t\tОшибка агрегации кластеров [prev]+[M%d]->[A%d] dL=%d\n",
                         nM-1, pcurnMut->iAggrCL, dLclu);
                break;
                
            case 33:     // "1011" Clust -> Mut  at AGGR
                if ( pcurnMut->iAggrCL != pnextMut->iAggrCL )   {
                    fprintf (tsld, "ER.%d '%s': Некластерная агрегация [M%d]->[A%d]\n", stat, desiz,
                             nM, pnextMut->iAggrCL);
                    break;
                }
                dLclu = getGapClust(pcurnMut, +1);  //Next
                if ( dLclu <= PorogClust )
                    continue;
                fprintf (tsld, "ER.%d '%s': Некластерная агрегация [M%d]->[A%d]\n", stat, desiz,
                             nM, pnextMut->iAggrCL);
                fprintf (tsld, "\t\tОшибка агрегации кластеров [С%d]+[next]->[A%d] dL=%d\n",
                             pcurnMut->iClust, pcurnMut->iAggrCL, dLclu);
                break;
                
                
                
            default:
                fprintf (tsld, "ER.%d '%s': Некластерная агрегация [%d - %d] = %d\n", stat, desiz, nM-1, nM, dLnuc);
                break;
        }
        fprintf (tsld, "\t%d\tC%d\tA%d\n\t%d\tC%d\tA%d\n",
                 pcurnMut->nucPos, pcurnMut->iClust, pcurnMut->iAggrCL,
                 pnextMut->nucPos, pnextMut->iClust, pnextMut->iAggrCL);
        errCnt++;

    }
    
    return errCnt;
}
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////

int XROSOMA:: getGapClust(MUTANT *pMUt, int next)   // +1 = next; -1 = prev
{
    int neighbour = ( next>0 ) ? pMUt->iClust+1 : pMUt->iClust-1;
    
    if ( neighbour >= vClust.size() || neighbour < 0 )
        return (int)Xsize;
  
    CLUSTER *pClu = &vClust[neighbour];
    
    int dL = ( next>0 ) ? vMutAPO[pClu->iBegMu].nucPos - pMUt->nucPos
                        : pMUt->nucPos - vMutAPO[pClu->iEndMu].nucPos;
    
    return dL;
}
//////////////////////////////////////////////////////////////////////////////////

