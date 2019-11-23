/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "libmycu/cudp/CuBatchDP_init_btck.cuh"
#include "libmycu/cudp/CuBatchDP_final.cuh"
#include "libmycu/cuss/CuBatchSS_com.h"
#include "CuBatchDP_aln_final.cuh"

// #define CUDP_ALN_FINAL_TESTPRINT 27 //449 //1286 //-1 //0//3024//4899

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for finalizing MAP dynamic programming;
// NOTE: thread block is supposed to consist of one warp;
// NOTE: keep #registers below 64!
// printsss, flag of whether to print secondary structure states along 
// alignments;
// logevthld, log e-value threshold;
// ndb1pros, number of profiles in the first profile data buffer db1;
// querprosOmtd, number of query profiles up to this query profile;
// ndb1prosOmtd, number of profiles missed up to the first one in db1;
// ndbCprosOmtd, number of profiles missed up to the first one in dbC;
// nqyposs, number of query positions;
// ndb1poss, number of cached db profile positions (phase 1);
// ndbCposs, number of new db profile positions (phase 1);
// querposoffset, offset from the origin of the device buffers allocated for 
// queries;
// bdb1posoffset, offset from the origin of the device buffers allocated for 
// cached db profile data;
// bdbCposoffset, offset from the origin of the device buffers allocated for 
// new (read) db profile data;
// dblen2, total number of positions of profiles passed to phase 2 PLUS 
// padding for this data;
// dbxpad2, number of padded positions for memory alignment;
// dbalnlen2, total size of alignments over all profiles beiong processed in 
// phase 2;
//

// -------------------------------------------------------------------------
// FinalizeDP_ALN: device code for finalizing dynamic programming in 
// parallel over profiles;
// NOTE: memory pointers should be aligned!
// scores, scores used as input for delineating positive substitution scores;
// tmpdpdiagbuffers, temporary buffers of calculated diagonal scores;
// maxscoordsbuf, coordinates of maximum alignment scores;
// btckdata, backtracking information data;
// dp2alndatbuffers, profile-specific statistics from (previous) phase 1;
// outalns, buffer for alignments themselves;
// 
__global__ void FinalizeDP_ALN(
    bool printsss,
    float logevthld,
    uint ndb1pros,
    uint /*querprosOmtd*/, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint /*bdbCposoffset*/,
    uint dblen2,
    uint dbxpad2,
    uint dbalnlen2,
    const CUBSM_TYPE* __restrict__ scores, 
    const float* __restrict__ tmpdpdiagbuffers,
    const uint* __restrict__ maxscoordsbuf,
    const char* __restrict__ btckdata,
    float* __restrict__ dp2alndatbuffers,
    char* __restrict__ outalns )
{
    LNTYPE dbprodstCache;//distance in positions to db profile blockIdx.x
    INTYPE dbprolenCache;//length of profile blockIdx.x
    //NOTE: keep SMEM below 2KB to consistently ensure high occupancy;
    //reuse maxscCache
    //__shared__ CUBSM_TYPE
    //        scoreCache[CUDP_2DCACHE_DIM_D+1];
    //cache of maximum scores from the last processed diagonals
    // (add 1 to dimensions to avoid bank conflicts):
    __shared__ CUBDP_TYPE maxscCache[TIMES2(CUDP_2DCACHE_DIM_D)+1];
    __shared__ int maxscndxCache[TIMES2(CUDP_2DCACHE_DIM_D)];//indices
    __shared__ float outPstsCoordsCache[3];//#posistives and beg-end coordinates
    __shared__ float outIdnsGaps[2];//#identities and gaps in the alignment
    __shared__ char outAlnCache[nTDP2OutputAlignmentSSS][CUDP_2DCACHE_DIM_D];
    //
    float maxsc;
    int ndx;
    //
    __shared__ uint proattr[2];//original profile NUMBER and new DISTANCE from the beginning
    float eval;//log e-value
    uint dbfldsndx;
    uint pronr;


    // blockIdx.x is the profile serial number in phase 2;
    //{{get the original profile number and e-value too
    if( threadIdx.x < 2 ) {
        proattr[threadIdx.x] = 
            *(uint*)(dp2alndatbuffers + nTDP2OutputAlnData*blockIdx.x+dp2oadOrgProNo+threadIdx.x);
        if( threadIdx.x == 0 )
            eval = dp2alndatbuffers[nTDP2OutputAlnData*blockIdx.x+dp2oadEvalue];
    }
    __syncwarp();
    eval = __shfl_sync(0xffffffff, eval, 0);
    //}}


    //all threads exit if calculated e-value exceeds the threshold
    if( logevthld < eval ) {
        //overwrite the dp2oadEpA field with the zero alignment length
        if( threadIdx.x == 0 )
            *(int*)(dp2alndatbuffers + nTDP2OutputAlnData*blockIdx.x+dp2oadEpA) = 0;
        return;
    }


    //NOTE: protection against overflow ensured on the host side
    if( proattr[0] < ndb1pros ) { pronr = proattr[0] + ndb1prosOmtd;
                dbfldsndx = pmv2DTotFlds;
    } else {    pronr = proattr[0] - ndb1pros + ndbCprosOmtd;//jump to section ndbCposs
                //bdb1posoffset = bdbCposoffset;
                dbfldsndx = TIMES2(pmv2DTotFlds);
    }

    if( threadIdx.x == 0 ) {
        dbprodstCache = ((LNTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DDist]))[pronr];
        dbprolenCache = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DLen]))[pronr];
    }

    //{{use registers efficiently
    __syncwarp();
    dbprodstCache = __shfl_sync(0xffffffff, dbprodstCache, 0);
    dbprolenCache = __shfl_sync(0xffffffff, dbprolenCache, 0);
    //}}


    if(dblen2 < proattr[1]+dbprolenCache+dbxpad2)
        return;

    const uint dblen = ndb1poss + ndbCposs + dbxpad;
    //dbpos is the beginning x position in the score matrix;
    const uint dbpos = (proattr[0] < ndb1pros)? dbprodstCache - bdb1posoffset: dbprodstCache + ndb1poss;
    //db profle position calculated from the end, including the offset determined by tid:
    int dbposoff = dbpos + dbprolenCache-1-threadIdx.x;
    //dbpos2 is the starting position of the profile in phase 2
    int dbpos2 = proattr[1];
    //db profle position (phase 2) calculated from the end, including the offset determined by tid:
    int dbpos2off = dbpos2 + dbprolenCache-1-threadIdx.x;


    bool qrysssinuse, trgsssinuse;

    if( threadIdx.x == 0 ) {
        //check whether SSS ionformation is present
        qrysssinuse = ((LNTYPE*)(dc_pm2dvfields_[pmv2DAddrSS]))[querposoffset] != (LNTYPE)-1;
        trgsssinuse = ((LNTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrSS]))[dbprodstCache] != (LNTYPE)-1;
    }
    //(qry/trg)sssinuse are used only by thread 0; do not sync
    //qrysssinuse = __shfl_sync(0xffffffff, qrysssinuse, 0);
    //trgsssinuse = __shfl_sync(0xffffffff, trgsssinuse, 0);


    int qpos = dpdssDiagM * nTDPDiagScoreSubsections * dblen;
// // //     int qpos = dpdssDiagM * nTDPDiagScoreSubsections * dblen2;
    maxscCache[threadIdx.x] = CUBDP_Q(0);
    maxscCache[threadIdx.x+CUDP_2DCACHE_DIM_D] = CUBDP_Q(0);
    //read maximum alignment scores at the end
    if( threadIdx.x < dbprolenCache)
        maxscCache[threadIdx.x] = tmpdpdiagbuffers[qpos + dbposoff];
// // //         maxscCache[threadIdx.x] = tmpdpdiagbuffers[qpos + dbpos2off];
    if( threadIdx.x+CUDP_2DCACHE_DIM_D < dbprolenCache )
        maxscCache[threadIdx.x+CUDP_2DCACHE_DIM_D] = 
            tmpdpdiagbuffers[qpos + dbposoff-CUDP_2DCACHE_DIM_D];
// // //             tmpdpdiagbuffers[qpos + dbpos2off-CUMAPDP_2DCACHE_DIM_D];

    __syncthreads();

// if((int)pronr==CUDP_ALN_FINAL_TESTPRINT){
//     printf(" ALNTEST: bid= %u tid= %u: pronr= %u len= %d addr= %u score= %.6f xx= %x\n",
//         blockIdx.x,threadIdx.x,pronr,
//         dbprolenCache,dbprodstCache,maxscCache[threadIdx.x],maxscoordsbuf[dbposoff]
//     );
//     if( threadIdx.x+CUDP_2DCACHE_DIM_D < dbprolenCache )
//     printf(" ALNTEST: bid= %u tid= %u: pronr= %u len= %d addr= %u score= %.6f xx= %x\n",
//         blockIdx.x,threadIdx.x+CUDP_2DCACHE_DIM_D,pronr,
//         dbprolenCache,dbprodstCache,maxscCache[threadIdx.x+CUDP_2DCACHE_DIM_D],
//         maxscoordsbuf[dbposoff-CUDP_2DCACHE_DIM_D]
//     );
//     for(size_t _k=0;_k<10000000000UL;_k++)clock();
// }

    //find the maximum score and its index (by index propagation);
    //all pairs are independent to avoid race conditions;
    dpfinmaxndxinit( maxscCache, maxscndxCache, threadIdx.x, threadIdx.x+CUDP_2DCACHE_DIM_D );
    __syncthreads();
    for(int i = CUDP_2DCACHE_DIM_D >> 1; i >= 1; i >>= 1 ) {
        dpfinmaxndx( maxscCache, maxscndxCache, threadIdx.x, threadIdx.x+i );
        __syncthreads();
    }

    maxsc = maxscCache[0];
    ndx = maxscndxCache[0];

    int xx = 0, x, y;
    //# identities, # positive substitution scores, # gaps:
    int idts = 0, psts = 0, gaps = 0;
    char btck = dpbtckDIAG;
    //read the coordinates of the maximum score
    if( 0.0f < maxsc && threadIdx.x == 0 )
        xx = maxscoordsbuf[dbposoff-ndx];
// // //         xx = maxscoordsbuf[dbpos2off-ndx];
    //coordinates xx, x, and y are used only by tid 0
    //xx = __shfl_sync(0xffffffff, xx, 0);
    y = GetCoordY(xx);
    x = GetCoordX(xx);


#ifdef CUDP_ALN_FINAL_TESTPRINT
    if((int)pronr==CUDP_ALN_FINAL_TESTPRINT){
        printf(" ALN Fin: bid= %u tid= %u: pronr= %u len= %d addr= %u maxsc= %.6f ndx= %d  "
            "(%x: y= %d x= %d)\n",
            blockIdx.x,threadIdx.x,pronr,
            dbprolenCache,dbprodstCache,maxsc,ndx,xx,y,x
        );
    }
#endif


    int alnlen = 0;

    //set the termination symbol for character strings;
    //NOTE: use the number of db positions calculated for the phase-2 output
    if( threadIdx.x == 0) {
        int i;
        qpos = dbpos2off + (blockIdx.x+1) * nqyposs;
        #pragma unroll
        for( i = 0; i < nTDP2OutputAlignment; i++, qpos += dbalnlen2 )
            outalns[qpos] = 0;
        if( printsss ) {
            #pragma unroll
            for( ; i < nTDP2OutputAlignmentSSS; i++, qpos += dbalnlen2 )
                outalns[qpos] = 0;
        }
    }

    //backtrace over the alignment to save it;
    while( btck != dpbtckSTOP ) {
        if( threadIdx.x == 0 ) {
            //only thread 0 caches data
            for( ndx = 0; ndx < blockDim.x; ndx++ ) {
                if( x < 0 || y < 0 ) {
                    btck = dpbtckSTOP;
                    break;
                }
                qpos = y * dblen + dbpos + x;
// // //                 qpos = y * dblen2 + dbpos2 + x;
                btck = btckdata[qpos];//READ
                if( btck == dpbtckSTOP )
                    break;
                outAlnCache[dp2oaQuerySSS][ndx] = ' ';
                outAlnCache[dp2oaTargetSSS][ndx] = ' ';
                if( btck == dpbtckUP ) { 
                    gaps++; 
                    outAlnCache[dp2oaQuery][ndx] = (char)//READ
                        ((CHTYPE*)(dc_pm2dvfields_[pmv2Daa]))[querposoffset+y];
                    outAlnCache[dp2oaTarget][ndx] = '-';
                    outAlnCache[dp2oaMiddle][ndx] = ' ';
                    if( printsss && qrysssinuse )
                        outAlnCache[dp2oaQuerySSS][ndx] = (char)//READ
                            ((CHTYPE*)(dc_pm2dvfields_[pmv2DSSstate]))[querposoffset+y];
                    y--; 
                    continue; 
                }
                else if( btck == dpbtckLEFT ) { 
                    gaps++; 
                    outAlnCache[dp2oaQuery][ndx] = '-';
                    outAlnCache[dp2oaTarget][ndx] = (char)//READ
                        ((CHTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2Daa]))[dbprodstCache+x];
                    outAlnCache[dp2oaMiddle][ndx] = ' ';
                    if( printsss && trgsssinuse )
                        outAlnCache[dp2oaTargetSSS][ndx] = (char)//READ
                            ((CHTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DSSstate]))[dbprodstCache+x];
                    x--; 
                    continue; 
                }
                //(btck == dpbtckDIAG)
                outAlnCache[dp2oaQuery][ndx] = (char)//READ
                    ((CHTYPE*)(dc_pm2dvfields_[pmv2Daa]))[querposoffset+y];
                outAlnCache[dp2oaTarget][ndx] = (char)//READ
                    ((CHTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2Daa]))[dbprodstCache+x];
                if( printsss ) {
                    if( qrysssinuse )
                        outAlnCache[dp2oaQuerySSS][ndx] = (char)//READ
                            ((CHTYPE*)(dc_pm2dvfields_[pmv2DSSstate]))[querposoffset+y];
                    if( trgsssinuse )
                        outAlnCache[dp2oaTargetSSS][ndx] = (char)//READ
                            ((CHTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DSSstate]))[dbprodstCache+x];
                }
                qpos = y * dblen + dbpos + x;
                if( outAlnCache[dp2oaQuery][ndx] == outAlnCache[dp2oaTarget][ndx] ) {
                    outAlnCache[dp2oaMiddle][ndx] = outAlnCache[dp2oaQuery][ndx];
                    idts++;
                } else if( scores[qpos] > 0.0f ) {//READ
                    outAlnCache[dp2oaMiddle][ndx] = '+';
                    psts++;
                } else
                    outAlnCache[dp2oaMiddle][ndx] = ' ';
                x--; y--;
            }
        }
        __syncwarp();
        btck = __shfl_sync(0xffffffff, btck, 0);
        ndx = __shfl_sync(0xffffffff, ndx, 0);
        //next, write part of the alignment to GMEM
        if( threadIdx.x < ndx ) {
            //allocated alignment length for a pair is l_query + l_target;
            //write the alignment in the reverse order starting at the end of the 
            // section dedicated to this db profile;
            // dbpos2off points to the end minus thread index;
            //subtract the number of positions that has been written already:
            qpos = dbpos2off + (blockIdx.x+1) * nqyposs  - alnlen-1;//-1 for termination 0
            outalns[qpos] = outAlnCache[dp2oaQuery][threadIdx.x]; qpos += dbalnlen2;
            outalns[qpos] = outAlnCache[dp2oaMiddle][threadIdx.x]; qpos += dbalnlen2;
            outalns[qpos] = outAlnCache[dp2oaTarget][threadIdx.x]; qpos += dbalnlen2;
            if( printsss ) {
                outalns[qpos] = outAlnCache[dp2oaQuerySSS][threadIdx.x]; qpos += dbalnlen2;
                outalns[qpos] = outAlnCache[dp2oaTargetSSS][threadIdx.x]; qpos += dbalnlen2;
            }
        }
        //alignment length has increased by ndx
        alnlen += ndx;
    }


    if( threadIdx.x == 0 ) {
        outIdnsGaps[0] = idts;
        outIdnsGaps[1] = gaps;
        outPstsCoordsCache[0] = psts;
        //beginning coordinates: beg state is always a match, add 1
        *(uint*)(outPstsCoordsCache+1) = CombineCoords(x+1,y+1);
        *(uint*)(outPstsCoordsCache+2) = xx;//alignment end coordinates
    }

    __syncwarp();

    //overwrite alignment information obtained in phase 1
    if( threadIdx.x < 3 ) {
        dp2alndatbuffers[nTDP2OutputAlnData*blockIdx.x+dp2oadPstvs+threadIdx.x] = 
            outPstsCoordsCache[threadIdx.x];
        if( threadIdx.x < 2 )
            dp2alndatbuffers[nTDP2OutputAlnData*blockIdx.x+dp2oadIdnts+threadIdx.x] = 
                outIdnsGaps[threadIdx.x];
        //overwrite the dp2oadEpA field with the alignment length
        if( threadIdx.x == 0 )
            *(int*)(dp2alndatbuffers + nTDP2OutputAlnData*blockIdx.x+dp2oadEpA) = alnlen;
    }

#ifdef CUDPMAPDP_FINAL_ALNS_AT_BEG
    //rewrite the saved alignment so that it starts at the beginning of the db 
    // profile section;
    //alnlen+1 to include the termination character;
    for( ndx = 0; ndx < alnlen+1; ndx += blockDim.x) {
        int i;
        //first, read part of the alignment
        if( threadIdx.x < alnlen+1-ndx ) {
            qpos = dbpos2 + (blockIdx.x+1) * nqyposs  + dbprolenCache-1 - alnlen + ndx;
            #pragma unroll
            for( i = 0; i < nTDP2OutputAlignment; i++, qpos += dbalnlen2 )
                outAlnCache[i][threadIdx.x] = outalns[qpos+threadIdx.x];
            if( printsss ) {
                #pragma unroll
                for( ; i < nTDP2OutputAlignmentSSS; i++, qpos += dbalnlen2 )
                    outAlnCache[i][threadIdx.x] = outalns[qpos+threadIdx.x];
            }
        }
        __syncwarp();
        if( threadIdx.x < alnlen+1-ndx ) {
            //next, shift it towards the beginning and write it 
            qpos = dbpos2 + blockIdx.x * nqyposs  + ndx;
            #pragma unroll
            for( i = 0; i < nTDP2OutputAlignment; i++, qpos += dbalnlen2 )
                outalns[qpos+threadIdx.x] = outAlnCache[i][threadIdx.x];
            if( printsss ) {
                #pragma unroll
                for( ; i < nTDP2OutputAlignmentSSS; i++, qpos += dbalnlen2 )
                    outalns[qpos+threadIdx.x] = outAlnCache[i][threadIdx.x];
            }
        }
        //no need to sync: next reading is beyond the indices used in this iteration
    }
#endif


__syncwarp();


#ifdef CUDP_ALN_FINAL_TESTPRINT
    __syncwarp();
    if( threadIdx.x == 0 ) {
        if((int)pronr==CUDP_ALN_FINAL_TESTPRINT){
            printf(" ALN Fin0 (pronr=%u): maxsc= %.6f y= %d x= %d alnlen= %d(%d) "
                "dbpos2off= %d blockIdx.x= %u nqyposs= %u "
                "idts= %d psts= %d gaps= %d  "
                "beg= (%u,%u) end= (%x: %u,%u)\n",
                pronr,maxsc,y,x,alnlen,
                *(int*)(dp2alndatbuffers + nTDP2OutputAlnData*blockIdx.x+dp2oadEpA),
                dbpos2off, blockIdx.x, nqyposs,
                (int)dp2alndatbuffers[nTDP2OutputAlnData*blockIdx.x+dp2oadIdnts],
                (int)dp2alndatbuffers[nTDP2OutputAlnData*blockIdx.x+dp2oadPstvs],
                (int)dp2alndatbuffers[nTDP2OutputAlnData*blockIdx.x+dp2oadNGaps],
                GetCoordY(*(uint*)(dp2alndatbuffers + nTDP2OutputAlnData*blockIdx.x+dp2oadBegCoords)),
                GetCoordX(*(uint*)(dp2alndatbuffers + nTDP2OutputAlnData*blockIdx.x+dp2oadBegCoords)),
                *(uint*)(outPstsCoordsCache+2),
                GetCoordY(*(uint*)(dp2alndatbuffers + nTDP2OutputAlnData*blockIdx.x+dp2oadEndCoords)),
                GetCoordX(*(uint*)(dp2alndatbuffers + nTDP2OutputAlnData*blockIdx.x+dp2oadEndCoords))
            );
            int i;
            qpos = dbpos2off + (blockIdx.x+1) * nqyposs - alnlen;
#   ifdef CUDPMAPDP_FINAL_ALNS_AT_BEG
            qpos = dbpos2 + blockIdx.x * nqyposs;
#   endif
            for( i = 0; i < nTDP2OutputAlignment; i++, qpos += dbalnlen2 )
                printf(" ALN Fin0 (pronr=%u): %s\n",pronr,outalns+qpos);
            if( printsss ) {
                #pragma unroll
                for( ; i < nTDP2OutputAlignmentSSS; i++, qpos += dbalnlen2 )
                    printf(" ALN Fin0 (pronr=%u): %s\n",pronr,outalns+qpos);
            }
        }
    }
#endif
}
