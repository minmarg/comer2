/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "liblib/msg.h"
#include "liblib/root.h"
// #include "libpro/srcpro/DistributionMatrix.h"
#include "AbstractScoreMatrix.h"
#include "ScoresAttr.h"

// -------------------------------------------------------------------------
// constructor
//
template<typename TScore>
ScoresAttr<TScore>::ScoresAttr(
    AbstractScoreMatrix<TScore>* prnt,
    TProbabilityFunction probfun,
    float reflambda,
    float refH,
    float refK,
    const TScore** scores,
    int scaled_by_factor )
:
    parent_( prnt ),
    parent_probfunction_( probfun ),
    parent_scores_( scores ),
    scaled_by_factor_( scaled_by_factor ),
    calculate_HK_ns_( true ),

    probabilities_( NULL ),
    min_score_( 0 ),
    max_score_( 0 ),
    score_gcd_( 1 ),

    priv_prob_vector_( NULL ),
    prob_vector_size_( 0 ),

    probabilities_ns_( NULL ),
    min_score_ns_( 0 ),
    max_score_ns_( 0 ),
    score_gcd_ns_( 1 ),

    priv_prob_vector_ns_( NULL ),
    prob_vector_size_ns_( 0 ),

    expscore_( 0.0f ),
    allnegatives_( false ),
    lambda_( -1.0f ),
    entropy_( 0.0f ),
    parameterK_( -1.0f ),

    queryLen_( 0 ),
    subjectLen_( 0 )
{
    SetRefLambda( reflambda );
    SetRefH( refH );
    SetRefK( refK );
    NewPrivateProbVector( PARAMETER_K_MAXIT, MAX_RANGE_OF_SCORES );
    NewPrivateProbVectorNS( PARAMETER_K_MAXIT, MAX_RANGE_OF_SCORES >> 3 );
}

// -------------------------------------------------------------------------
// default construction
//
template<typename TScore>
ScoresAttr<TScore>::ScoresAttr()
:
    parent_( NULL ),
    parent_probfunction_( NULL ),
    parent_scores_( NULL ),
    scaled_by_factor_( 1 ),
    calculate_HK_ns_( true ),

    probabilities_( NULL ),
    min_score_( 0 ),
    max_score_( 0 ),
    score_gcd_( 1 ),

    priv_prob_vector_( NULL ),
    prob_vector_size_( 0 ),

    probabilities_ns_( NULL ),
    min_score_ns_( 0 ),
    max_score_ns_( 0 ),
    score_gcd_ns_( 1 ),

    priv_prob_vector_ns_( NULL ),
    prob_vector_size_ns_( 0 ),

    expscore_( 0.0f ),
    allnegatives_( false ),
    lambda_( -1.0f ),
    entropy_( 0.0f ),
    parameterK_( -1.0f ),

    queryLen_( 0 ),
    subjectLen_( 0 )
{
    throw MYRUNTIME_ERROR("ScoresAttr::ScoresAttr: Default initialization is prohibited.");
}

// -------------------------------------------------------------------------
// destructor
//
template<typename TScore>
ScoresAttr<TScore>::~ScoresAttr()
{
    DestroyProbabilities();
    DestroyProbabilitiesNS();
    DestroyPrivateProbVector();
    DestroyPrivateProbVectorNS();
}

// -------------------------------------------------------------------------
// Init: initializes class's private data given query and subject sizes
// IMPORTANT: query and subject sizes are supposed to be initialized
//
template<typename TScore>
void ScoresAttr<TScore>::Init( int querylen, int sbjctlen )
{
    SetQuerySize( querylen );
    SetSubjectSize( sbjctlen );
}

// -------------------------------------------------------------------------
// SearchForHSPs: perform search of multiple high-scoring pairs (HSPs)
//     in the same diagonal of the score system;
// minhspscore, minimum HSP score;
// hsplen, HSP length;
// nohsps, minimum number of HSPs required;
// maxdist, max distance between HSPs;
// possbjct, posquery, subject and query positions of HSPs found in the
//     same diagonal (to be returned);
// return true if such HSPs have been found
//
template<typename TScore>
bool ScoresAttr<TScore>::SearchForHSPs(
    TScore minhspscore, int hsplen, int nohsps, int maxdist,
    int* possbjct, int* posquery )
{
    const mystring  preamb = "ScoresAttr::SearchForHSPs: ";
    bool            bret = false;
    const TScore    hmin = minhspscore;
    TScore          hmax = 0;
    const int       diaglen = nohsps * 3;//nohsps*2 for indices
    TScore**        diagonals = NULL;//diagonals keeping best-scoring hsps
    TScore**        hsps = NULL;
    int**           lens = NULL;
    TScore          loc_score, tvalue;
    myruntime_error mre;
    int             m, n, ind, sum, k, l;
    int             msb, nqu/*, mind*/;

    if( hmin <= 0 ) {
        throw MYRUNTIME_ERROR( preamb + "Non-positive HSP score threshold." );
    }
    if( hsplen <= 0 ) {
        throw MYRUNTIME_ERROR( preamb + "Non-positive HSP length." );
    }
    if( nohsps <= 0 ) {
        throw MYRUNTIME_ERROR( preamb + "Non-positive number of HSPs." );
    }
    if( maxdist <= 0 ) {
        throw MYRUNTIME_ERROR( preamb + "Non-positive distance between HSPs." );
    }
    if( GetSubjectSize() < 1 || GetQuerySize() < 1 )
        return false;

    try {
        diagonals = ( TScore** )malloc( sizeof(TScore*) * ( GetSubjectSize() + GetQuerySize() - 1 ));
        hsps = ( TScore** )malloc( sizeof(TScore*) * GetSubjectSize());
        lens = ( int** )malloc( sizeof(int*) * GetSubjectSize());

        if( !diagonals || !hsps || !lens )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory." );

        for( m = 0; m < GetSubjectSize() + GetQuerySize() - 1; m++ ) {
            diagonals[m] = ( TScore* )malloc( sizeof(TScore) * diaglen );
            if( !diagonals[m])
                throw MYRUNTIME_ERROR( preamb + "Not enough memory." );
            memset( diagonals[m], 0, sizeof(TScore) * diaglen );
        }
        for( m = 0; m < GetSubjectSize(); m++ ) {
            hsps[m] = ( TScore* )malloc( sizeof(TScore) * GetQuerySize());
            lens[m] = ( int* )malloc( sizeof(int) * GetQuerySize());

            if( !hsps[m] || !lens[m] )
                throw MYRUNTIME_ERROR( preamb + "Not enough memory." );

            memset( hsps[m], 0, sizeof(TScore) * GetQuerySize());
            memset( lens[m], 0, sizeof(int) * GetQuerySize());
        }

        //heuristics of several high-scoring pairs in the diagonal
        for( m = 0; m < GetSubjectSize() && !bret; m++ ) {
            for( n = 0; n < GetQuerySize() && !bret; n++ )
            {
                loc_score = GetParentScore( m, n );

                if( loc_score <= SCORE_MIN ) {
                    hsps[m][n] = 0;
                    continue;
                }

                hsps[m][n] = loc_score;
                lens[m][n] = 0;

                if( 0 < loc_score ) {
                    //begin or continue hsp
                    hsps[m][n] = loc_score;
                    lens[m][n] = 1;
                }

                if( 0 < m && 0 < n )
                    //boundary hsp positions should be positive
//                   (( ( lens[m-1][n-1] <= 1 || lens[m-1][n-1] == hsplen - 1 ) && 0 < hsps[m-1][n-1] ) ||
//                   ( 1 < lens[m-1][n-1] && lens[m-1][n-1] < hsplen - 1 ) ))
                {
                    tvalue = hsps[m][n] + hsps[m-1][n-1];
                    if( hsps[m][n] < tvalue && 0 < hsps[m][n]) {
                        if( lens[m-1][n-1] < hsplen ) {
                            hsps[m][n] = tvalue;
                            lens[m][n] = lens[m-1][n-1] + 1;
                        }
                    }
                }

                //record hsp score in the diagonal
                if( hsplen <= lens[m][n] && hmin <= hsps[m][n] )
                {
                    ind = GetQuerySize() - 1 - n + m;

                    sum = 0;
                    for( k = 0; k < nohsps; k++ )
                    {
                        if( diagonals[ind][k+nohsps] &&
                            maxdist < m - diagonals[ind][k+nohsps] - hsplen ) {
                            //distance between two hsps is unsatisf.: omit the hsp;
                            //all hsps are considered in turn
                            sum = -1;//to break
                            break;
                        }
                        if( diagonals[ind][k] < hsps[m][n] )
                        {
                            for( l = nohsps - 1; k < l && diagonals[ind][l-1] <= 0; l-- );
                            for( ; k < l; l-- ) {
                                if( maxdist < m - diagonals[ind][l-1+nohsps] - hsplen ) {
                                    //distance between two hsps is unsatisf.: omit the hsp;
                                    sum = -1;//to break
                                    break;
                                }
                                diagonals[ind][l] = diagonals[ind][l-1];//save hsp score
                                diagonals[ind][l+nohsps] = diagonals[ind][l-1+nohsps];//save coordinates..
                                diagonals[ind][l+nohsps+nohsps] = diagonals[ind][l-1+nohsps+nohsps];
                                sum += diagonals[ind][l];
                            }
                            if( sum < 0 )
                                break;
                            diagonals[ind][k] = hsps[m][n];
                            diagonals[ind][k+nohsps] = m;
                            diagonals[ind][k+nohsps+nohsps] = n;
                            sum += diagonals[ind][k];
                            break;
                        }
                        else
                            sum += diagonals[ind][k];
                    }
                    if( diagonals[ind][nohsps-1] && hmax < sum ) {
                        bret = true;//at least one group of hsps is found
                        hmax = sum;
                        //mind = ind;//NOTE: required for testing block
                        msb = diagonals[ind][nohsps];//positions of max of hsps
                        nqu = diagonals[ind][nohsps+nohsps];
                        if( possbjct )  *possbjct = msb;
                        if( posquery )  *posquery = nqu;
                    }
                }
            }
        }
    } catch( myexception const& ex ) {
        mre = ex;
    }

//{{TESTING BLOCK
// GetParent()->PrintScoringMatrix( stderr );
// fprintf( stderr, "\n\n");
// for( n = 0; n < GetQuerySize(); n++ ) {
//     if( !n ) {
//         fprintf( stderr, "****");
//         for( m = 0; m < GetSubjectSize(); m++ )
//             fprintf( stderr, " %5d", m );
//         fprintf( stderr, "\n");
//         for( m = 0; m < GetSubjectSize(); m++ )
//             fprintf( stderr, "------" );
//         fprintf( stderr, "\n");
//     }
//     fprintf( stderr, "%3d:", n );
//     for( m = 0; m < GetSubjectSize(); m++ ) {
//         fprintf( stderr, " %5.1f", ( float )hsps[m][n] / GetScoresFactor());
//     }
//     fprintf( stderr, "\n");
// }
// fprintf( stderr, "\n%d %5.1f: %d %d\n\n", bret, ( float )hmax/* / GetScoresFactor()*/, nqu, msb );
// if( bret ) {
//     for( k = 0; k < nohsps; k++ ) {
//         fprintf( stderr, "  score=%-.1f at (%d,%d)\n",
//             ( float )diagonals[mind][k]/*/GetScoresFactor()*/,
//             diagonals[mind][k+nohsps+nohsps], diagonals[mind][k+nohsps]);
//     }
//     fprintf( stderr, "--------\n");
//     for( ind = 0; ind < GetSubjectSize() + GetQuerySize() - 1; ind++ ) {
//         if( diagonals[ind][nohsps-1] <= 0 )
//             continue;
//         for( k = 0; k < nohsps; k++ ) {
//             fprintf( stderr, "  score=%-.1f at (%d,%d)\n",
//                 ( float )diagonals[ind][k]/*/GetScoresFactor()*/,
//                 diagonals[ind][k+nohsps+nohsps], diagonals[ind][k+nohsps]);
//         }
//         fprintf( stderr, "\n");
//     }
// }
//}}

    for( m = GetSubjectSize(); m < GetSubjectSize() + GetQuerySize() - 1; m++ ) {
        if( diagonals && diagonals[m]) { free( diagonals[m]); diagonals[m] = NULL; }
    }
    for( m = 0; m < GetSubjectSize(); m++ ) {
        if( diagonals && diagonals[m]) { free( diagonals[m]); diagonals[m] = NULL; }
        if( hsps && hsps[m]) { free( hsps[m]); hsps[m] = NULL; }
        if( lens && lens[m]) { free( lens[m]); lens[m] = NULL; }
    }
    if( diagonals ) { free( diagonals ); diagonals = NULL; }
    if( hsps ) { free( hsps ); hsps = NULL; }
    if( lens ) { free( lens ); lens = NULL; }

    if( mre.isset())
        throw mre;

    return bret;
}



// -------------------------------------------------------------------------
// ComputeStatisticalParameters: compute all statistical parameters
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeStatisticalParameters( bool computelambda, bool wrn )
{
    if( computelambda )
        ComputeProbabilitiesAndLambda( wrn );
    ComputeKHGivenLambda();

// GetParent()->PrintScoringMatrix( stderr );
// fprintf( stderr, "  >  lambda=%f, entropy=%f, exp.score=%f\n", GetPrivLambda(), GetEntropy(), GetPrivExpectedScore());

}

// -------------------------------------------------------------------------
// ComputeProbabilitiesAndLambda: compute score probabilities and
//     statistical parameter lambda
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeProbabilitiesAndLambda( bool wrn )
{
    ComputeScoreProbabilities();
    if( GetExpectedScore() < 0.0 )
        ComputeLambda( wrn );
}

// -------------------------------------------------------------------------
// ComputeKHGivenLambda: compute entropy H and Karlin's parameter K
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeKHGivenLambda()
{
    try {
        if( GetExpectedScore() < 0.0f && 0.0f < GetLambda()) {
            ComputeEntropyGivenLambda();
            ComputeKarlinsK();
        }
    } catch( myexception const& ex ) {
        //DestroyProbabilities();
        throw ex;
    }
    //this method is to be last-called in the procedure of computation of staistical parameters
    //we don't need probabilities any more
    //DestroyProbabilities();
}



// -------------------------------------------------------------------------
// ComputeScoreProbabilities: compute probabilities to observe scores
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeScoreProbabilities()
{
    ComputeScoreProbabilitiesHelper();
    if( GetUsingNSscoresforHK())
        ComputeScoreProbabilitiesNS();
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilitiesHelper: compute probabilities to observe scores
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeScoreProbabilitiesHelper()
{
    const mystring  preamb = "ScoresAttr::ComputeScoreProbabilitiesHelper: ";
    const float     accuracy = 1.e-4f;
    char            strbuf[BUF_MAX];
    float           prob, consv;
    int m, n, sc;

#ifdef __DEBUG__
    if( GetParent() == NULL || GetParentProbFunction() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null parent or its function." );
#endif

    if( GetQuerySize() < 1 || GetSubjectSize() < 1 )
        return;

    if( !GetParentScores())
        throw MYRUNTIME_ERROR( preamb + "Null parent scores." );

    TScore  l_minscore = 0,
            l_maxscore = 0;
    TScore  loc_score  = 0;

    for( m = 0; m < GetSubjectSize(); m++ ) {
        for( n = 0; n < GetQuerySize(); n++ ) {

            loc_score = GetParentScore( m, n );

            if( loc_score <= SCORE_MIN )
                continue;

            if( loc_score < l_minscore ) l_minscore = loc_score;
            if( l_maxscore < loc_score ) l_maxscore = loc_score;
        }
    }
    SetMinMaxScores( l_minscore, l_maxscore );

    ( GetParent()->*GetParentProbFunction())( this );

    //verify probabilities
    consv = 0.0f;
    for( sc = GetMinScore(); sc <= GetMaxScore(); sc ++ ) {
        prob = GetProbabilityOf( sc );
        if( !prob )
            continue;
        consv += prob;
    }
    if( consv < 1.0f - accuracy || consv > 1.0f + accuracy ) {
        sprintf( strbuf, "Invalid score probabilities: sum = %f", consv );
        throw MYRUNTIME_ERROR( preamb + strbuf );
    }
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilitiesNS: compute probabilities to observe non-scaled 
// scores; NOTE: the probabilities are calculated based on the probabilities
// obtained for scaled scores
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeScoreProbabilitiesNS()
{
    const mystring  preamb = "ScoresAttr::ComputeScoreProbabilitiesNS: ";
    const int fct = GetScoresFactor();
    float prob;
    int sc, scns;

    if( fct < 2 )
        return;

    TScore minsc = GetMinScore();
    TScore maxsc = GetMaxScore();
    TScore minscns = (TScore)rintf((float)minsc / (float)fct);
    TScore maxscns = (TScore)rintf((float)maxsc / (float)fct);
    
    SetMinMaxScoresNS( minscns, maxscns );

    for( sc = minsc; sc <= maxsc; sc ++ ) {
        prob = GetProbabilityOf( sc );
        if( !prob )
            continue;
        scns = (TScore)rintf((float)sc / (float)fct);
        IncProbabilityOfNS( scns, prob );
    }

#ifdef __DEBUG__
    const float accuracy = 1.e-4f;
    char strbuf[BUF_MAX];
    float consv = 0.0f;
    //probabilities have been already verified for scaled scores;
    //NOTE: this verification is redundant
    for( scns = minscns; scns <= maxscns; scns ++ ) {
        prob = GetProbabilityOfNS( scns );
        if( !prob )
            continue;
        consv += prob;
    }
    if( consv < 1.0f - accuracy || consv > 1.0f + accuracy ) {
        sprintf( strbuf, "Invalid score probabilities: sum = %f", consv );
        throw MYRUNTIME_ERROR( preamb + strbuf );
    }
#endif
}


// -------------------------------------------------------------------------
// ComputeGCD: compute the greatest common divisor of scores
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeGCD()
{
    int gcd = 1;

    gcd = -GetMinScore();
    for( int sc = GetMinScore() + 1; sc <= (int)GetMaxScore() && gcd > 1; sc++ )
        if( 0.0f < GetProbabilityOf( sc ))
            gcd = mygcd( gcd, sc - GetMinScore());

    if( gcd == 0 )
        return;

    SetGCD( gcd );
}

// -------------------------------------------------------------------------
// ComputeGCDNS: compute the greatest common divisor of non-scaled scores
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeGCDNS()
{
    int gcd = 1;

    gcd = -GetMinScoreNS();
    for( int sc = GetMinScoreNS() + 1; sc <= (int)GetMaxScoreNS() && gcd > 1; sc++ )
        if( 0.0f < GetProbabilityOfNS( sc ))
            gcd = mygcd( gcd, sc - GetMinScoreNS());

    if( gcd == 0 )
        return;

    SetGCDNS( gcd );
}

// -------------------------------------------------------------------------
// ComputeLambda: compute statistical parameter lambda
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeLambda( bool wrn )
{
    const mystring preamb = "ScoresAttr::ComputeLambda: ";
    mystring msgbuf = preamb;
    float newlambda = -1.0f;

    SetLambda( newlambda );

    if( 0.0f <= GetExpectedScore()) {
        msgbuf += "Expected score per position is non-negative.";
        warning( msgbuf.c_str());
        return;
    }

    //rare will be cases when interval between adjacent scores is >1;
    //nevertheless, compute the greatest common divisor of all scores;
    //change of variables then will be made via substitution y=exp(-lambda gcd)
    ComputeGCD();
    if( GetUsingNSscoresforHK())
        ComputeGCDNS();

    newlambda = findLambdaRoot(
        GetGCD(), 0.0f + EPSILON, 10.0f, LAMBDA_ACCURACY, MAXIT
    );

    if( newlambda <= 0.0f + EPSILON ) {
        //newlambda = STABLE.StatisParam( Ungapped, Lambda );//to leave score values unchanged
        if( wrn ) {
            msgbuf += "Computation of lambda failed";
            if( GetParent() && GetParent()->GetName())
                    msgbuf += mystring(": ") + GetParent()->GetName();
            else    msgbuf += ".";
            warning( msgbuf.c_str());
        }
        //return;
    }

    SetLambda( newlambda );
}

// -------------------------------------------------------------------------
// findLambdaRoot: helper method to find the root of lambda
//
template<typename TScore>
float ScoresAttr<TScore>::findLambdaRoot( int gcd,
    float x1, float x2, float xacc, int maxit )
{
    float rts = -1.0f;
    const char* retmsg = NULL;
    //{{params
    char    params[BUF_MAX];
    char*   p = params;
    *(void**)p = (void*)this; p += sizeof(void*);
    *(int*)p = gcd;
    //}}

    x1 = expf( -x1 );
    x2 = expf( -x2 );

    retmsg = root_by_NR_and_bisection(
        &lmbd_FdFfunc,
        x1,
        x2,
        xacc,
        maxit,
        (void*)params,
        &rts 
    );

    if( retmsg || rts <= 0.0f )
        return rts;

    return -logf(rts) / (float)gcd;
}

// -------------------------------------------------------------------------
// lmbd_FdFfunc: function describing probability conservation equation
//     which is used to find the root of lambda;
// When fp number is used as a type of scores, then the conservation 
// function is computed as
//
//     1     __  lambda s(m,n)
//   ------ \   e             - 1 = 0
//   l1 l2  /__
//          m,n
//
// Such computation is very slow and perhaps is not accurate.
// When integer scores are used, the function is
//
//    _   lambda s(k)
//   \  e            p(s(k)) - 1
//   /_
//    k
//
// where p(s(k)) are probabilities of scores s(k). To find a solution, a
// procedure as in psi-blast is used: the equation is transformed into
// the polynomial form,
//     exp(max s(k) lambda) * poly(exp(-lambda))
// where poly is a polynomial of exp(-lambda) and has exactly two zeros:
// 1 and the other in interval (0,1)
//
template<typename TScore>
void ScoresAttr<TScore>::lmbd_FdFfunc( float x, float* f, float* df, void* params )
{
    float y = x;//initially x is supposed to be exp(-lambda)
    float ldf = 0.0f;
    float lf = 0.0f;
    //{{params
    char* p = (char*)params;
    ScoresAttr<TScore>* pthis = *(ScoresAttr<TScore>**)p; p += sizeof(void*);
    int gcd = *(int*)p;
    //}}

    for( int sc = pthis->GetMinScore(); sc <= (int)pthis->GetMaxScore(); sc += gcd ) {
        ldf = ldf * y + lf;
        lf  = lf  * y + pthis->GetProbabilityOf(sc);
        //zero will always be reached
        if( sc == 0 )
            lf -= 1.0f;
    }
    *df = ldf;
    *f  = lf;
}



// -------------------------------------------------------------------------
// ComputeEntropyGivenLambda: compute relative entropy given lambda; 
//  We compute entropy before adjusting the scores.
//  Mathematically the result would be obtained the same if computing
//  entropy after the score adjustment, since the probabilities of scores
//  remain the same. The only difference is that the reference lambda
//  instead of computed lambda should be used.
// Relative entropy or information per aligned pair of positions is computed
// as follows:
//
//           _                   lambda s(k)
//   lambda \  s(k) * p(s(k)) * e
//          /_
//          k
//
// where p(s(k)) is a probability of the score s(k)
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeEntropyGivenLambda()
{
    float lmbd = GetLambda();

    if( lmbd < 0.0f ) {
        warning( "ScoresAttr::ComputeEntropyGivenLambda: Unable to compute relative entropy." );
        return;
    }

    TScore (ScoresAttr<TScore>::*getminfunc)() const = &ScoresAttr<TScore>::GetMinScore;
    TScore (ScoresAttr<TScore>::*getmaxfunc)() const = &ScoresAttr<TScore>::GetMaxScore;
    float (ScoresAttr<TScore>::*getprobfunc)(TScore) const = &ScoresAttr<TScore>::GetProbabilityOf;

    if( GetUsingNSscoresforHK()) {
        getminfunc = &ScoresAttr<TScore>::GetMinScoreNS;
        getmaxfunc = &ScoresAttr<TScore>::GetMaxScoreNS;
        getprobfunc = &ScoresAttr<TScore>::GetProbabilityOfNS;
        lmbd *= (float)GetScoresFactor();
    }

    float   y = expf( -lmbd );
    float   m = powf( y, (float)(this->*getmaxfunc)()); // e to -lambda * max_score
    float   H = 0.0f;

    for( int sc = (this->*getminfunc)(); sc <= (int)(this->*getmaxfunc)(); sc ++ )
        H = H * y + sc * (this->*getprobfunc)( sc );

    if( 0.0f < m )
        H /= m;
    else
        if( 0.0f < H )
            H = expf( lmbd * (this->*getmaxfunc)() + logf(H));

    SetH( lmbd * H );
}

// -------------------------------------------------------------------------
// ComputeKarlinsK: computes parameter K analitically as given by
//     Karlin et al. PNAS 87 (1990)
// K is an extreme value distribution parameter related to its location
// parameter mu by relation
//   mu = ln (Kmn) / lambda
// K is computed by the formula
//
//       gcd lambda exp( -2 sigma )
//   K = --------------------------
//       H (1 - exp( -gcd lambda ))
//
// where gcd is the greatest common divisor of the scores, H is relative
// entropy, and sigma is roughly a probability to obtain an arbitrary
// alignment of any length >=1 given scoring scheme
//
//            _   1 (  _           i lambda     _          )
//   sigma = \   -- ( \   P[j](i) e          + \   P[j](i) )
//           /_   j ( /_                       /_          )
//           j>0      i<0                     i>=0
//
// P[j](i) is a probability to obtain an alignment of score i and length j
//              _
//   P[j](i) = \  P[1](k) P[j-1](i-k)
//             /_
//              k
//
// There are three simplifications: if high score = gcd and low score = -gcd,
//
//   K = squared( P[1](gcd) - P[1](-gcd) )  / p[1](-gcd)
//
// if high score = gcd and low score <> -gcd,
//
//   K = H (1 - exp(-gcd lambda )) / (gcd lambda)
//
// if high score <> gcd and low score = -gcd,
//                                                          _
//   K = lambda (1 - exp(-gcd lambda )) / (gcd H) squared( \  i P[1](i) )
//                                                         /_
//                                                        i=l:u
// -------------------------------------------------------------------------
//
template<typename TScore>
void ScoresAttr<TScore>::ComputeKarlinsK()
{
    const mystring  preamb = "ScoresAttr::ComputeKarlinsK: ";
    float   loc_lambda = GetLambda();
    float   loc_expected_score = GetExpectedScore();

    if( 0.0f <= loc_expected_score ) {
        warning( mystring( preamb + "Non-negative expected score.").c_str());
        return;
    }
    if( loc_lambda < 0.0f || GetH() <= 0.0f ) {
        warning( mystring( preamb + "Unable to compute K.").c_str());
        return;
    }

    float K = -1.0f;

    const int       maxrange    = MAX_RANGE_OF_SCORES;
    const int       maxit       = PARAMETER_K_MAXIT;
    const float     precision   = PARAMETER_K_ACCURACY;

    TScore (ScoresAttr<TScore>::*getminfunc)() const = &ScoresAttr<TScore>::GetMinScore;
    TScore (ScoresAttr<TScore>::*getmaxfunc)() const = &ScoresAttr<TScore>::GetMaxScore;
    TScore (ScoresAttr<TScore>::*getgcdfunc)() const = &ScoresAttr<TScore>::GetGCD;
    float (ScoresAttr<TScore>::*getprobfunc)(TScore) const = &ScoresAttr<TScore>::GetProbabilityOf;
    float* (ScoresAttr<TScore>::*getvectfunc)() const = &ScoresAttr<TScore>::GetPrivateProbVector;
    void (ScoresAttr<TScore>::*resetvectfunc)(size_t) = &ScoresAttr<TScore>::ResetPrivateProbVector;

    if( GetUsingNSscoresforHK()) {
        getminfunc = &ScoresAttr<TScore>::GetMinScoreNS;
        getmaxfunc = &ScoresAttr<TScore>::GetMaxScoreNS;
        getgcdfunc = &ScoresAttr<TScore>::GetGCDNS;
        getprobfunc = &ScoresAttr<TScore>::GetProbabilityOfNS;
        getvectfunc = &ScoresAttr<TScore>::GetPrivateProbVectorNS;
        resetvectfunc = &ScoresAttr<TScore>::ResetPrivateProbVectorNS;
        loc_lambda *= (float)GetScoresFactor();
    }



    TScore  gcd = (this->*getgcdfunc)();
    TScore  minscore = (this->*getminfunc)();
    TScore  maxscore = (this->*getmaxfunc)();
    float probminscore = (this->*getprobfunc)(minscore);
    float probmaxscore = (this->*getprobfunc)(maxscore);

    TScore  range = maxscore - minscore;
    float   y = expf( -gcd * loc_lambda );

    if( maxrange < range )
        throw MYRUNTIME_ERROR2( preamb + "Too large range of scores.", SCALING );

    if( gcd <= 0 || range <= 0 )
        throw MYRUNTIME_ERROR2( preamb + "Invalid range of scores.", SCALING );

    if( minscore == -gcd && maxscore == gcd ) {
        if( probminscore )
            K = SQUARE( probminscore - probmaxscore ) / probminscore;
        SetK( K );
        return;
    }

    if( minscore == -gcd || maxscore == gcd ) {
        if( maxscore != gcd )
            K = loc_lambda / ( gcd * GetH()) * ( 1.0f - y ) * SQUARE( loc_expected_score );
        else
            K = GetH() / ( gcd * loc_lambda ) * ( 1.0f - y );
        SetK( K );
        return;
    }

    float   sigma = 0.0f;   //probability of alignment of any length
    float   Pj  = 1.0f;     //probability to obtain arbitrary alignment of length j
    float   Pji = 0.0f;     //probability to obtain alignment of length j with score i
    float*  Pjivector =                 //vector of alignment score probabilities P[j](i),
        (this->*getvectfunc)();         // where j is an alignment length, and i is a score
    TScore  min_j_score = 0;    //minimum score of alignment of length j
    TScore  max_j_score = 0;    //maximum score of alignment of length j
    TScore  lower = 0;          //lower bound of scores
    TScore  upper = 0;          //upper bound of scores

    size_t  maxPjisize = maxit * range + 1;//plus 1 to include zero

    if( !Pjivector ) {
        error( mystring( preamb + "Space for probability vector is not allocated.").c_str());
        return;
    }

    (this->*resetvectfunc)( maxPjisize );
    Pjivector[0] = 1.0f;    //in order to compute the lowest score probability of alignment with length=1

    //compute sigma, probability of alignment of any length
    // Pjivector is reused in each iteration to restore probabilities of
    // alignment of length one smaller
    for( int j = 1; j <= maxit && precision < Pj; j++ )
    {
        TScore  i;
        min_j_score += minscore;
        max_j_score += maxscore;
        lower = maxscore;
        upper = maxscore;

        //compute Pj for all possible scores i
        for( i = max_j_score; min_j_score <= i; i -= gcd )
        {
            Pji = 0.0f;
            //P[j](i) = SUM P[j-1](i-k) P[1](k)
            for( int k = lower; k <= upper; k += gcd ) {
                Pji += Pjivector[(i-min_j_score)-(k-minscore)] * (this->*getprobfunc)( k );
            }

            if( minscore < lower )
                lower -= gcd;
            if( i - min_j_score <= range )
                upper -= gcd;

            Pjivector[i-min_j_score] = Pji;
        }

        Pj = 0.0f;
        //compute SUM (i<0) ( P[j](i) exp(i lambda) );
        //actually we need to compute polynomial of exp(lambda) with
        // coefficients P[j](i)
        for( i = min_j_score; i < 0; i += gcd )
            Pj = Pj * y + Pjivector[i-min_j_score];
        Pj *= y;

        //compute SUM (i>=0) ( P[j](i) );
        //just sum previously computed probabilities P[j](i)
        for( ; i <= max_j_score; i += gcd )
            Pj += Pjivector[i-min_j_score];

        sigma += Pj / j;
    }

    //we have all variables found to compute K
    K = ( gcd * loc_lambda * expf( -2.0f * sigma )) / ( GetH() * ( 1.0f - y ));
    SetK( K );
}

// -------------------------------------------------------------------------
// PrintProbabilities: print score probabilities if available
//
template<typename TScore>
void ScoresAttr<TScore>::PrintProbabilities( FILE* fp )
{
    if( fp == NULL )
        return;

    if( !probabilities_ )
        return;

    fprintf( fp,"%s%5c Score probabilities%s", NL, 32, NL );

    fprintf( fp, "%9c", 32 );

    for( int sc = (int)GetMinScore(); sc <= (int)GetMaxScore(); sc++ ) {
        fprintf( fp, "%s%5d %6.4g", NL, (int)sc, GetProbabilityOf(sc));
    }
    fprintf( fp, "%s", NL );
}

// =========================================================================

template class ScoresAttr<int>;
