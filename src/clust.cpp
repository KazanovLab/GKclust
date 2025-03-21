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

#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

extern PROGARGS ArgKit;

extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

vector < long > vRndMut;               // posMUT for Random simulator
vector < pair <long, int> > vRndClu;   // <begPos, cLng>

vector < int > MutSpace;  // Space between Mutations

bool lesser_pair ( const pair<long,int> &x1, const pair<long,int> &x2 ) {
    return x1.first < x2.first;
}
//////////////////////////////////////////////////////////////////////////////////

void SAMPLE:: findClusters( int iXro, int *clustID, int Porog )
{
    int cntM;
    int MaxP = (int)(vMutAPO[iXro].size() -1);
    CLUSTER Clust;
    
    int nM=0;
    while  ( nM <= MaxP ) {
        cntM = 0;
        for ( int n=nM+1; n < MaxP; n++,cntM++ ) {
            if ( (vMutAPO[iXro][n].nucPos -vMutAPO[iXro][n-1].nucPos) > Porog )
                break;
        }
        if ( cntM > 0 ) {
            Clust.iBegMu = nM;
            Clust.iEndMu = nM + cntM;
            Clust.cLng = (int)(vMutAPO[iXro][Clust.iEndMu].nucPos - vMutAPO[iXro][Clust.iBegMu].nucPos);
            Clust.cID   = *clustID;
            vClust[iXro].push_back(Clust);
            for ( int n=0; n <= cntM; n++ ) {
                vMutAPO[iXro][nM+n].iClust = (int)vClust[iXro].size() - 1 ;
            }
        }

        if ( cntM > 0 )
            (*clustID)++ ;
        nM += (cntM+1);
    }
    return;
    
}
///////////////////////////////////////////////////////////////////////////////////

void SAMPLE:: ClustConnector(  int iXro, int *clustID, int Porog  ) 
{
    int cnt;
    int MaxCL = (int)(vClust[iXro].size() -1);
    CLUSTER CombinedCL;
    int headMut, tailMut_1;
    
    int nCL=0;
    while  ( nCL <= MaxCL ) {
        cnt = 0;
        for ( int n=nCL+1; n < MaxCL; n++,cnt++ ) {
            headMut = vClust[iXro][n].iBegMu;
            tailMut_1 = vClust[iXro][n-1].iEndMu;
            if ( (vMutAPO[iXro][headMut].nucPos -vMutAPO[iXro][tailMut_1].nucPos) > Porog )
                break;
        }
        if ( cnt > 0 ) {
            CombinedCL.iBegMu = vClust[iXro][nCL].iBegMu;
            CombinedCL.iEndMu = vClust[iXro][nCL+cnt].iEndMu;
            CombinedCL.cLng = (int)(vMutAPO[iXro][CombinedCL.iEndMu].nucPos -
                                    vMutAPO[iXro][CombinedCL.iBegMu].nucPos);
            CombinedCL.cID   = *clustID;
            vAggrC[iXro].push_back(CombinedCL);
            for ( int n=CombinedCL.iBegMu; n <= CombinedCL.iEndMu; n++ )    {
                vMutAPO[iXro][n].iAggrCL = (int)vAggrC[iXro].size() - 1 ;
            }
        }

        if ( cnt > 0 )
            (*clustID)++ ;
        
        nCL += (cnt+1);
    }
    return;


}
//////////////////////////////////////////////////////////////////////////////////


int SAMPLE:: calcMutPorog( int iXro ) //, float PartOfMut)
{
    int Xsiz = (int)vecDNK[iXro].Xsize;
    int Msiz = (int)vMutAPO[iXro].size();
    unsigned int nM, L;
    unsigned long Rnd_pos;
    vector <long>::iterator Iter;
    
    if ( Msiz < 3 )
        return 0;
    
    if ( vRndMut.capacity() < Msiz )
        vRndMut.reserve(maxMUTsize+1);
    if ( MutSpace.capacity() < maxMUTsize )
        MutSpace.reserve(RANDOM_CYCLES * maxMUTsize + 1);

    MutSpace.clear();
    
    int cntN=0, cntDubl=0;
    for ( int nCycl=0; nCycl<RANDOM_CYCLES; nCycl++ )    {
        nM = 0;
        vRndMut.clear();
        while (nM < Msiz) {
            Rnd_pos = ( (double)rand() / (double)RAND_MAX ) * Xsiz;
            if ( *(vecDNK[iXro].Xbody+Rnd_pos)=='N' )   {
                cntN++;
                continue;
            }
            if ( vRndMut.empty() || vRndMut.back() < Rnd_pos )   {      // ? <=  | < ?
                vRndMut.push_back( Rnd_pos );
                nM++;
                continue;
            }
            Iter = lower_bound( vRndMut.begin( ), vRndMut.end( ), Rnd_pos );   // MUTANT(NucPos),  lesser_MUT);
            if ( *Iter == Rnd_pos ) {
                cntDubl++;
                continue;
            }
            vRndMut.insert( Iter, Rnd_pos );
            nM++;
        }


        for ( nM=1; nM<Msiz; nM++ ) {
            L = (int)(vRndMut[nM] - vRndMut[nM-1]);
            MutSpace.push_back(L);
        }
        
    }    // of limon cycles
    
    sort(MutSpace.begin(), MutSpace.end());
    
    float fP;
    sscanf(ArgKit.chPart, "%f", &fP);
    
    int indPor = (fP/100) * RANDOM_CYCLES  * (Msiz-1);
    int Porog  = MutSpace[indPor];
    
    if ( GET_ARG_D(ArgKit.argTAG) ) {
        char chP[16];
        sprintf(chP, "%d", Porog);
        int step = (int)(MutSpace.size() / 1000);
        nM = 0;
        fprintf(ArgKit.foutRndMu, "Chr\tChrSize\tnMut\tSpace\tGap\n");
        while ( nM < MutSpace.size() ) {
            fprintf(ArgKit.foutRndMu,"%d\t%d\t%d\t%d\t%s\n", vecDNK[iXro].chrNum, Xsiz,
                    Msiz, MutSpace[nM], ( ( abs((int)nM-indPor) <= step ) ? chP : " " ) );
                //, SampName.c_str());
            nM += step;
        }

    }
    
    return Porog;

}
/////////////////////////////////////////////////////////////////////////

int SAMPLE:: calcClustPorog( int iXro ) //, float PartOfMut)
{
    int Xsiz = (int)vecDNK[iXro].Xsize;
    int nCL;
    unsigned long Rnd_pos;
    pair < long, int > Rnd_CL;
    vector < pair<long,int> >::iterator Iter;

    CLUSTER *pCLU;
    
    int CLsiz = (int)vClust[iXro].size();  //vMutAPO[iXro].size();
    
    if ( CLsiz < 3 )
        return 0;
    
    MutSpace.clear();
    
    int cntN=0, cntDubl=0;
    for ( int nCycl=0; nCycl<RANDOM_CYCLES; nCycl++ )    {
        nCL = 0;
        vRndClu.clear();
        while (nCL < CLsiz) {
            Rnd_pos = ( (double)rand() / (double)RAND_MAX ) * Xsiz;
            pCLU = &vClust[iXro][nCL];
            if ( (Rnd_pos + pCLU->cLng) >= Xsiz )
                continue;
            Rnd_CL.first  = Rnd_pos;
            Rnd_CL.second = pCLU->cLng;
            if ( vecDNK[iXro].CLtest_N( Rnd_CL ) > 0 )  {
                cntN++;
                continue;
            }
            if ( vRndClu.empty() )  {       //}|| vRndMut.back() < Rnd_pos )   {
                vRndClu.push_back( Rnd_CL );
                nCL++;
                continue;
            }
            
            Iter = lower_bound( vRndClu.begin( ), vRndClu.end( ), Rnd_CL, lesser_pair);
            if ( vecDNK[iXro].CLtest_Cross(Rnd_CL,  Iter) > 0 )  {
                cntDubl++;
                continue;
            }
            vRndClu.insert( Iter, Rnd_CL );
            nCL++;
        }
        
        
        int L;
        for ( nCL=1; nCL<CLsiz; nCL++ ) {
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
    
    int indPor = (fP/100) * RANDOM_CYCLES  * (CLsiz-1);
    int Porog  = MutSpace[indPor];
    
    if ( GET_ARG_D(ArgKit.argTAG) ) {
        char chP[16];
        sprintf(chP, "%d", Porog);
        int step = (int)(MutSpace.size() / 1000);
        nCL = 0;
        fprintf(ArgKit.foutRndCl, "Chr\tChrSize\tnMut\tSpace\tGap\n");
        while ( nCL < MutSpace.size() ) {
            fprintf(ArgKit.foutRndCl,"%d\t%d\t%d\t%d\t%s\n", vecDNK[iXro].chrNum, Xsiz,
                    CLsiz, MutSpace[nCL], ( ( abs((int)nCL-indPor) <= step ) ? chP : " " ) );
                //vecSAMPL[iSamp].SampName.c_str());
            nCL += step;
        }

    }


    return Porog;
    
}
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

void SAMPLE::printCluTrace(  )
{
    char cName_1[16], cName_2[16];
    char cMutN_1[8], cMutN_2[8];
    char cPorog[8], cDprio[8], cDnext[8];
    int stat_1=0, stat_2=0;   // =1 if NOT Clust;  =2 if itis Clust & 1_st Mut ;  =3 Clust & next Mut
    CLUSTER *pCLU, *pAggr;

    fprintf(ArgKit.foutTrace ,"Chr\tPos\tRef\tAlt\t"); //APO\t");
    fprintf(ArgKit.foutTrace,"ClustID_1\tMutCnt_1\tGap_Mut\tSpace_up\tSpace_down\t");
    fprintf(ArgKit.foutTrace,"ClustID_2\tMutCnt_2\tGap_Clust\tSpace_up\tSpace_down\n");
    
    for ( int nX=0; nX < NO_HUMAN_XRO; nX++ )    {
        for ( int nMu=0; nMu < vMutAPO[nX].size(); nMu++ )   {
            pCLU =  (vMutAPO[nX][nMu].iClust < 0) ? NULL : &vClust[nX][vMutAPO[nX][nMu].iClust];
            pAggr = (vMutAPO[nX][nMu].iAggrCL < 0) ? NULL : &vAggrC[nX][vMutAPO[nX][nMu].iAggrCL];
            strcpy(cMutN_1, " "); strcpy(cMutN_2, " ");
            strcpy(cName_1, " "); strcpy(cName_2, " ");
            //---------
            if ( vMutAPO[nX][nMu].iClust < 0 )
                stat_1 =1;
            else
                stat_1 = ( nMu==0 || (vMutAPO[nX][nMu].iClust != vMutAPO[nX][nMu-1].iClust) ) ? 2 : 3;
            
            if ( vMutAPO[nX][nMu].iAggrCL < 0 )
                stat_2 =1;
            else
                stat_2 = ( nMu==0 || (vMutAPO[nX][nMu].iAggrCL != vMutAPO[nX][nMu-1].iAggrCL) ) ? 2 : 3;
            
            //fprintf(fMutClu,"Xr\tPOS\tREF\tALT\t  //APO
            fprintf(ArgKit.foutTrace,"%d\t%d\t%c\t%c\t",    //%s\t",
                    nX+1, (int)vMutAPO[nX][nMu].nucPos,
                    vMutAPO[nX][nMu].nucREF, vMutAPO[nX][nMu].nucALT );
            //                    ( (vMutAPO[nX][nMu].APOtag == 0 ) ? " " : "APO" ) );
            
            //----- my_Clust
            if ( stat_1 != 1 )
                sprintf(cName_1, "cID_%d", pCLU->cID );
            if ( stat_1 == 2 )
                sprintf(cMutN_1, "%d", (pCLU->iEndMu - pCLU->iBegMu) + 1 );
            
            if ( stat_1 == 2  )  //
                sprintf(cPorog, "%d", PorogMut[nX]);
            else
                strcpy(cPorog, " ");
            
            if ( stat_1 != 1 )  { //}|| statX != 1 ) {
                sprintf(cDprio, "%ld", (  (nMu==0)
                                        ? vMutAPO[nX][nMu].nucPos
                                        : vMutAPO[nX][nMu].nucPos - vMutAPO[nX][nMu-1].nucPos) );
                sprintf(cDnext, "%ld", ( (nMu == vMutAPO[nX].size()-1)
                                        ? vecDNK[nX].Xsize - vMutAPO[nX][nMu].nucPos
                                        : vMutAPO[nX][nMu+1].nucPos - vMutAPO[nX][nMu].nucPos) );
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
                sprintf(cPorog, "%d", PorogClust[nX]);
            else
                strcpy(cPorog, " ");
            
            if ( stat_2 != 1 ) {
                sprintf(cDprio, "%ld", (  (nMu==0)
                                        ? vMutAPO[nX][nMu].nucPos
                                        : vMutAPO[nX][nMu].nucPos - vMutAPO[nX][nMu-1].nucPos) );
                sprintf(cDnext, "%ld", ( (nMu == vMutAPO[nX].size()-1)
                                        ? vecDNK[nX].Xsize - vMutAPO[nX][nMu].nucPos
                                        : vMutAPO[nX][nMu+1].nucPos - vMutAPO[nX][nMu].nucPos) );
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

void SAMPLE::printCluMut( int clustOnly, FILE * fRezult )
{
    //    FILE * fRezult;
    char clustID[16];
    char cntMut[8];
    int stat_1=0, stat_2=0;   // =1 if NOT Clust;  =2 if itis Clust & 1_st Mut ;  =3 Clust & next Mut
    CLUSTER *pCLU, *pAggr;
    /*
     fRezult = fopen( fName.c_str(), "w");
     if ( fRezult==NULL ) {
     printf ("\nINV_CreateFile =""%s"" not found\n", fName.c_str() );
     return ;
     }
     */
    fprintf(fRezult, "Chr\tPos\tRef\tAlt\t"); //APO\t");
    fprintf(fRezult, "ClusterID\tMutCount\n");
    
    
    for ( int nX=0; nX < NO_HUMAN_XRO; nX++ )    {
        for ( int nMu=0; nMu < vMutAPO[nX].size(); nMu++ )   {
            if ( vMutAPO[nX][nMu].iClust < 0 )
                stat_1 = 1;
            else
                stat_1 = ( nMu==0 || (vMutAPO[nX][nMu].iClust != vMutAPO[nX][nMu-1].iClust) ) ? 2 : 3;
            
            if ( vMutAPO[nX][nMu].iAggrCL < 0 )
                stat_2 = 1;
            else
                stat_2 = ( nMu==0 || (vMutAPO[nX][nMu].iAggrCL != vMutAPO[nX][nMu-1].iAggrCL) ) ? 2 : 3;
            
            if (stat_1==1 && stat_2==1) {
                if ( ! clustOnly )
                    fprintf(fRezult,"%d\t%d\t%c\t%c\n", // %s\t",
                            nX+1, (int)vMutAPO[nX][nMu].nucPos,
                            vMutAPO[nX][nMu].nucREF, vMutAPO[nX][nMu].nucALT );
                continue;
            }
            
            pCLU =  (vMutAPO[nX][nMu].iClust < 0)  ? NULL : &vClust[nX][vMutAPO[nX][nMu].iClust];
            pAggr = (vMutAPO[nX][nMu].iAggrCL < 0) ? NULL : &vAggrC[nX][vMutAPO[nX][nMu].iAggrCL];
            //---------
            //fprintf(fMutClu,"#S\tXr\tPOS\tREF\tALT\t  //APO
            fprintf(fRezult,"%d\t%d\t%c\t%c\t", // %s\t",
                    nX+1, (int)vMutAPO[nX][nMu].nucPos,
                    vMutAPO[nX][nMu].nucREF, vMutAPO[nX][nMu].nucALT );
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
            fprintf(fRezult,"%s\t%s\n", clustID, cntMut ); //SampName.c_str()
            
        }
    }
    return;
}

//===============================================================================


//////////////////////////////////////////////////////////////////////////////////
