/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <stdlib.h>

#include "liblib/msg.h"
#include "liblib/alpha.h"
#include "segdata.h"
#include "SEGAbstract.h"

namespace SEG {

//STATICS:
//default print width
const size_t SEGAbstract::scszDefSegSeqPrintWidth = 60;
//default number of segments
const size_t SEGAbstract::scszDefNoSegments = 20;

// -------------------------------------------------------------------------
// STATIC FUNCTIONS
//
// SimpleIntComparer: comparison function for state vector values
//
static int SimpleIntComparer( const void* key1, const void* key2 )
{
    return (int)( (ssize_t)key1 - (ssize_t)key2 );
}

// /////////////////////////////////////////////////////////////////////////
// CLASS LinearSearchStructure
//
// constructor
//
LinearSearchStructure::LinearSearchStructure( 
    TComparator compfunc, TEqComparer eqfunc, size_t size )
:   BinarySearchStructure( compfunc, size, false/*do not keep duplicates*/),
    eqcomparer( eqfunc )
{
    if( eqcomparer == NULL )
        throw MYRUNTIME_ERROR( "LinearSearchStructure::LinearSearchStructure: "
                               "Null comparison function." );
}

// destructor
//
LinearSearchStructure::~LinearSearchStructure()
{
}

// -------------------------------------------------------------------------
// Find: find a key in the structure and set the location of the found key;
// if the key is not found, then set the location to point to the last 
// element
//
bool LinearSearchStructure::Find( const void* key, int* loc ) const
{
    int n = 0;

    for( n = 0; n < (int)GetSize(); n++ )
        if(( *eqcomparer )( GetValueAt(n), key )) {
            if( loc )
                *loc = n;
            return true;
        }

    if( loc )
        *loc = n;

    return false;
}


// /////////////////////////////////////////////////////////////////////////
// CLASS BinaryModifier
//
// constructor
//
BinaryModifier::BinaryModifier( TComparator compfunc, size_t size )
:   BinarySearchStructure( compfunc, size, true /*keep_duplicates*/)
{
}

// destructor
//
BinaryModifier::~BinaryModifier()
{
}

// -------------------------------------------------------------------------
// IncValue: increase count value; if several count elements with the same 
// value exist, modify the last one but less than the next element
//
bool BinaryModifier::IncValue( tcount value )
{
    size_t  loc;
    tcount  locvalue = 0;
    int     location = -1;

    if( !Find((void*)value, &location ) || location < 0 )
        return false;

    //value found, proceed
    for( loc = location; loc < GetSize() && ((tcount)GetValueAt(loc)) == value; loc++ );

    loc--;

    if( GetSize() <= loc )
        return false;

    //modification keeps the structure ordered
    locvalue = (tcount)GetValueAt(loc) + 1;
    SetValueAt( loc, (void*)locvalue );

    return true;
}

// -------------------------------------------------------------------------
// DecValue: reduce count value; if several count elements with the same 
// value exist, modify the first one but greater than the preceeding element
//
bool BinaryModifier::DecValue( tcount value )
{
    size_t  loc;
    tcount  locvalue = 0;
    int     location = -1;

    if( value == 0 )
        return false;

    if( !Find((void*)value, &location ) || location < 0 )
        return false;

    //value found, proceed
    for( loc = location; (locvalue=(tcount)GetValueAt(loc)) == value && loc; loc-- );

    if( locvalue != value )
        loc++;

    //modification keeps the structure ordered
    locvalue = (tcount)GetValueAt(loc) - 1;
    SetValueAt( loc, (void*)locvalue );

    return true;
}


// /////////////////////////////////////////////////////////////////////////
// CLASS Window
//
// constructor:
//
// NOTE: address is supposed to point to a vector of pointers void*
//
Window::Window(
    TVerifier verffunc,
    TComparator comparator,
    TEqComparer eqcomparer,
    const void* address,
    size_t len,
    size_t szalpha )
:
    verifier( verffunc ),
    window( address ),
    length( len ),
    szalphabet( szalpha ),

    entropy( -1.0f ),
    composition( NULL ),
    compaddresses( NULL ),
    state( NULL )
{
    if( verifier == NULL )
        throw MYRUNTIME_ERROR( "Window::Window: Null verifier function." );

    if( length < 1 )
        throw MYRUNTIME_ERROR( "Window::Window: Window length is 0." );

    size_t reservation = ( szalphabet == 0 )? length: szalphabet;

    composition = new SimpleVector( reservation );

    if( eqcomparer )
        compaddresses = new LinearSearchStructure( comparator, eqcomparer, reservation );
    else
        compaddresses = new BinarySearchStructure( comparator, reservation, false/*do not keep duplicates*/);

    state = new BinaryModifier( SimpleIntComparer, reservation );

    if( composition == NULL || compaddresses == NULL || state == NULL )
        throw MYRUNTIME_ERROR( "Window::Window: Not enough memory." );
}

// destructor
//
Window::~Window()
{
    if( composition )
        delete composition;

    if( compaddresses )
        delete compaddresses;

    if( state )
        delete state;
}

// -------------------------------------------------------------------------
// Initialize: compute composition, arrange state vector, and compute
// entropy
//
void Window::Initialize()
{
    CompositionOn();
    StateOn();
    ComputeEntropy();
}

// -------------------------------------------------------------------------
// CompositionOn: compute the composition of window
//
void Window::CompositionOn()
{
    for( size_t pos = 0; pos < GetLength(); pos++ )
        RefreshComposition( GetWindowAt( pos ));
}

// -------------------------------------------------------------------------
// StateOn: arrange ordered state vector
//
void Window::StateOn()
{
    for( size_t n = 0; n < GetCompositionSize(); n++ )
        PushStateValue( GetCompositionAt( n ));
}

// -------------------------------------------------------------------------
// ComputeEntropy: compute the entropy of window according to 
// information encoded in state vector;
// actually need to compute the sum of P = n[i]/L (log2 n[i]/L),
// where n[i] is a value from state vector; we have precomputed values of
//      p(n[i]) = n[i] log2 n[i];
// so the needed P = 1/L (p(n[i]) - n[i]/L p(L)).
// The sum of P is
//      1/L (SUM p(n[i]) - p(L)), since SUM n[i] == L
//
void Window::ComputeEntropy()
{
    float   partial = 0.0f;
    float   winentr = 0.0f;
    tcount  stateval = 0;
    tcount  L = 0;

    for( size_t n = 0; n < GetStateSize(); n++ ) {
        L += (stateval = GetStateAt(n));
        partial += PRT_ENTROPIES.GetValueOf((size_t)stateval );
    }

    if( 0 < L )
        winentr = ( partial-PRT_ENTROPIES.GetValueOf((size_t)L) ) / (float)L;
    SetEntropy( winentr );
}

// -------------------------------------------------------------------------
// LogProbability: compute the logarithm of probability of this window;
// the probability is equal to the product of multinomial term and the 
// number of compositions divided by the total number of available 
// sequences:
//     P = M * F / N^L,
// where M is the multinomial term, F is the compositional term, N is the 
// size of alphabet, and L is the length of window
//
float Window::LogProbability() const
{
    size_t  N = GetAlphabetSize();
    float   logP = 0.0f;

    if( N == 0 )
        //get window length then
        N = GetLength();

    logP = LogMultinomial() + LogNoCompositions(N) - GetLength() * LOGARITHMS.GetValueOf(N);
    return logP;
}

// -------------------------------------------------------------------------
// LogMultinomial: compute the logarithm of the multinomial term 
// included in the final probability of window composition; the result 
// gives the number of distinct sequences given composition vector
//
float Window::LogMultinomial() const
{
    float   denominator = 0.0f;
    tcount  stateval = 0;
    tcount  L = 0;

    for( size_t n = 0; n < GetStateSize(); n++ ) {
        L += (stateval = GetStateAt(n));
        denominator += LOG_FACT.GetValueOf((size_t)stateval );
    }

//     return LOG_FACT.GetValueOf(( size_t )L ) - denominator;
    return LOG_FACT.GetValueOf(GetLength()) - denominator;
}

// -------------------------------------------------------------------------
// LogNoCompositions: compute the logarithm of the number of distinct 
// compositions the state vector can have; the number is computed as
//     N! / (PRODUCT r[i] (N-K)!),
// N is the size of alphabet, PRODUCT is, up to C, the number of distinct 
// counts in state vector, r[i] is the repetition number of count i in state 
// vector;
//     SUM r[i] == K;
// K is the size of state vector (without the 0 elements).
//
float Window::LogNoCompositions( size_t N ) const
{
    if( N == 0 )
        //get the window length then
        N = GetLength();

    float result = LOG_FACT.GetValueOf( N );

    if( GetStateSize() == 0 )
        return result;

    float   denominator = 0.0f;
    size_t  repnum = 1;
    size_t  K = 0;
    tcount  stateval = 0;
    tcount  statevalp1 = GetStateAt(0);
    tcount  L = 0;

    for( size_t n = 0; n < GetStateSize(); n++ ) {
        L += (stateval = statevalp1);

        if( n + 1 < GetStateSize())
            statevalp1 = GetStateAt(n+1);

        if( stateval == 0 )
            continue;

        K++;

        if( n + 1 < GetStateSize() && stateval == statevalp1 ) {
            repnum++;
            continue;
        }

        denominator += LOG_FACT.GetValueOf( repnum );
        repnum = 1;
    }

#ifdef __DEBUG__
    if( N < K )
        throw MYRUNTIME_ERROR( "Window::LogNoCompositions: Wrong alphabet size." );
#endif

    if( 0 < K ) {
        denominator += LOG_FACT.GetValueOf( N - K );
        result -= denominator;
    }

    return result;
}

// -------------------------------------------------------------------------
// SafeShift: safely shift window by 1 position
//
void Window::SafeShift()
{
    if( Uninitialized()) {
        Initialize();
        return;
    }

    Shift();
}

// -------------------------------------------------------------------------
// Shift: shift window by 1 position to the right; refresh composition and 
// state vectors, recompute entropy
//
void Window::Shift()
{
    if( GetLength() < 1 )
        return;

    size_t      address = 0;
    tcount      begvalue;
    tcount      newvalue;
    const void* value = NULL;
    bool        changed = false;

    value = GetWindowAt( 0 );

    if(( *GetVerifier())( value ))
    {
        //0th element must be found
        if( !FindCompAddress( value, &address ))
            throw MYRUNTIME_ERROR( "Window::Shift: Composition value not found." );

        begvalue = GetCompositionAt( address ); //composition value at the beginning of window

        if( begvalue < 1 )
            throw MYRUNTIME_ERROR( "Window::Shift: Illegal composition value." );

        DecCompositionAt( address );    //decrease composition value since sliding window will change its beginning
        DecStateValue( begvalue );      //decrease state value at the beginning
        changed = true;
    }

    //reinitialize window at the next position
    window = GetWindowAt( 1 );

    value = GetWindowAt( GetLength() - 1 );

    if(( *GetVerifier())( value ))
    {
        if( PushCompAddress( value, &address )) {
            //the first occurence of the element
            InsertCompositionAt( address, 1 );
            PushStateValue( 1 );
        } else {
            //element found, get its composition value
            newvalue = GetCompositionAt( address ); //composition value at the boundary of window

            //newvalue can be 0 if it was met before
            IncCompositionAt( address );
            IncStateValue( newvalue );  //increase state of the newly met value
        }
        changed = true;
    }

    if( changed )
        ComputeEntropy();
}

// =========================================================================
// testing:
// Print: print window composition, state, and entropy
//
void Window::Print( FILE* fp )
{
    if( fp == NULL )
        return;

    fprintf( fp, "%sWindow%s Composition:%s ", NL, NL, NL );
    for( size_t n = 0; n < GetCompositionSize(); n++ )
        fprintf( fp, " %3zu", GetCompositionAt( n ));

    fprintf( fp, "%s State:%s ", NL, NL );
    for( size_t n = 0; n < GetStateSize(); n++ )
        fprintf( fp, " %3zu", GetStateAt( n ));

    fprintf( fp, "%s Entropy:%s  %.3f", NL, NL, GetEntropy());
    fprintf( fp, "%s Log-multinomial term:%s  %.4f", NL, NL, LogMultinomial());
    fprintf( fp, "%s Log-No-compositions term:%s  %.4f", NL, NL, LogNoCompositions());
    fprintf( fp, "%s Log-probability:%s  %.4f%s%s", NL, NL, LogProbability(), NL, NL );
}


// /////////////////////////////////////////////////////////////////////////
// CLASS SEGAbstract
//
// constructor:
//
// NOTE: address is supposed to point to a vector of pointers void*
//
SEGAbstract::SEGAbstract(
    TValidator  valdfunc,
    TVerifier   verffunc,
    TComparator compfunc,
    TEqComparer eqcomfunc,
    void*   address,
    size_t  runlen,
    size_t  winlen,
    float   lowent,
    float   highent,
    size_t  maxdiff,
    size_t  szalpha )
:
    validator( valdfunc ),
    verifier( verffunc ),
    comparator( compfunc ),
    eqcomparer( eqcomfunc ),
    runaddress( NULL ),
    runlength( 0 ),
    winlength( winlen ),
    szalphabet( szalpha ),
    lowentropy( lowent ),
    highentropy( highent ),
    maxdiffextent( maxdiff ),

    inverted( false ),

    entropies( NULL ),
    segments( scszDefNoSegments )
{
    if( validator == NULL || verifier == NULL || comparator == NULL )
        throw MYRUNTIME_ERROR( "SEGAbstract::SEGAbstract: Null function addresses." );

    SetRunAddress( address, runlen );
}

// constructor: alternative
//
SEGAbstract::SEGAbstract(
    TValidator  valdfunc,
    TVerifier   verffunc,
    TComparator compfunc,
    TEqComparer eqcomfunc,
    size_t  winlen,
    float   lowent,
    float   highent,
    size_t  maxdiff,
    size_t  szalpha )
:
    validator( valdfunc ),
    verifier( verffunc ),
    comparator( compfunc ),
    eqcomparer( eqcomfunc ),
    runaddress( NULL ),
    runlength( 0 ),
    winlength( winlen ),
    szalphabet( szalpha ),
    lowentropy( lowent ),
    highentropy( highent ),
    maxdiffextent( maxdiff ),

    inverted( false ),

    entropies( NULL ),
    segments( scszDefNoSegments )
{
    if( validator == NULL || verifier == NULL || comparator == NULL )
        throw MYRUNTIME_ERROR( "SEGAbstract::SEGAbstract: Null function addresses." );
}

// destructor
//
SEGAbstract::~SEGAbstract()
{
    Destroy();
}

// -------------------------------------------------------------------------
// SetRunAddress: reset run at the given address
//
void SEGAbstract::SetRunAddress( void* address, size_t runlen )
{
    Allocate( runlen );

    runaddress = address;
    *const_cast<size_t*>(&runlength) = runlen;

    if( runlength < winlength )
        warning( "SEGAbstract::SetRunAddress: "
                 "Ineffective configuration: Window length > run length." );
}

// -------------------------------------------------------------------------
// Allocate: allocate memory for entropy
//
void SEGAbstract::Allocate( size_t newlen )
{
    if( newlen < 1 )
        throw MYRUNTIME_ERROR( "SEGAbstract::Allocate: Run length is 0." );

    entropies = (float*)malloc( sizeof(float) * newlen );

    if( entropies == NULL )
        throw MYRUNTIME_ERROR( "SEGAbstract::Allocate: Not enough memory." );

    memset( entropies, 0, sizeof(float) * newlen );
}

// -------------------------------------------------------------------------
// Destroy: deallocate memory
//
void SEGAbstract::Destroy()
{
    if( entropies ) {
        free( entropies );
        entropies = NULL;
    }
}

// -------------------------------------------------------------------------
// Run: implement all steps of the algorithm
//
void SEGAbstract::Run()
{
    Initialize();
    FindSegments();
}

// -------------------------------------------------------------------------
// Initialize: initialize and prepare the run for searching low/high entropy
// regions
//
void SEGAbstract::Initialize()
{
    ComputeEntropy();
}

// -------------------------------------------------------------------------
// Entropy: compute segment entropy over the entire run length
//
float SEGAbstract::Entropy() const
{
    Window  window(
            GetVerifier(),
            GetComparator(),
            GetEqComparer(),
            GetRunAddress(),
            GetRunLength(),// !!- GetWinLength(),
            GetAlphabetSize()
    );
    window.SafeShift();
    return window.GetEntropy();
}

// -------------------------------------------------------------------------
// LogProbability: compute log-probability of segment over the entire run
// length
//
float SEGAbstract::LogProbability() const
{
    Window  window(
            GetVerifier(),
            GetComparator(),
            GetEqComparer(),
            GetRunAddress(),
            GetRunLength(),// !!- GetWinLength(),
            GetAlphabetSize()
    );
    window.SafeShift();
    return window.LogProbability();
}

// -------------------------------------------------------------------------
// FindSegments: the main method to be called to find segments of low/high
//     entropy
//
void SEGAbstract::FindSegments()
{
    if( GetRunLength() == 0 )
        return;

    RecurrentSegFind( 0, GetRunLength() - 1 );
}

// -------------------------------------------------------------------------
// RecurrentSegFind: find segments of low/high entropy within the run 
// bounded by left and right
//
void SEGAbstract::RecurrentSegFind( size_t left, size_t right )
{
    size_t  n;
    size_t  locleft, locright;
    size_t  optleft, optright;
    size_t  centre =   GetWinLength()>> 1;
    size_t  adjust = ( GetWinLength() & 1 ) ^ 1;    //for even win lengths, adjust
//     size_t  negadj = ( GetWinLength() & 1 );        //negative adjust
    float   entro;

    if( GetRunLength() <= left || GetRunLength() <= right )
        throw MYRUNTIME_ERROR( "SEGAbstract::RecurrentSegFind: Wrong extent boundaries." );

    if( right < left )
        return;

    for( n = left; n <= right; n++ )
    {
        entro = GetEntropyAt( n );
        if( entro < 0.0f )
            continue;

        if( GetInverted()?
            entro < GetHighEntropy() :
            GetLowEntropy() < entro  )
            continue;

        locleft = locright = n;
        FindExtentAt( n, left, right, &locleft, &locright );

        if( locleft + adjust < centre || GetRunLength() <= locright + centre )
            throw MYRUNTIME_ERROR( "SEGAbstract::RecurrentSegFind: Wrong extent found." );

        //draw out to window boundaries
        optleft = locleft = locleft - centre + adjust;
        optright = (locright += centre);

        OptimizeExtent( &optleft, &optright );

        //if the found optimal extent is beyond the current position,
        //recurse to the left side to find left-hand side segments
        if( n + centre <= optleft && 1 + centre <= optleft ) {
            RecurrentSegFind( locleft + centre - adjust, optleft - 1 - centre );
        }

        //push segment...
        segments.Push( optleft, optright );
        //move next from the leftmost window found
        n = PCMIN( locright - centre, optright + centre - adjust );
        left = n + 1;
    }
}

// -------------------------------------------------------------------------
// ComputeEntropy: compute the entropy of each window within the run
//
void SEGAbstract::ComputeEntropy()
{
    size_t  n, p, r, l;
    size_t  centre = GetWinLength() >> 1;
    size_t  adjust = GetWinLength()? ((GetWinLength() & 1) ^ 1): 0;//for even win lengths, adjust
//     size_t  negadj = ( GetWinLength() & 1 );        //negative adjust
    size_t  invalid = (size_t)-1;               //last invalid position

    Window  window(
            GetVerifier(),
            GetComparator(),
            GetEqComparer(),
            GetRunAddress(),
            GetWinLength(),
            GetAlphabetSize()
    );

    //reset the beginning of entropy calculation to the position of window centre
    for( n = 0; n < GetRunLength() && n < centre - adjust; n++ )
    {
        SetEntropyAt( n, -1.0f );
        //if position is not valid to consider
        if( ! ( *GetValidator())( GetRunAt( n )))
            invalid = n;
    }

    //find the last invalid position within the window
    for( p = n; p < GetRunLength() && p <= n + centre; p++ )
        if( ! ( *GetValidator())( GetRunAt( p )))
            invalid = p;

    for( ;  n < GetRunLength() &&
            n + centre < GetRunLength();
            n++ )
    {
        window.SafeShift();
//         l = n - centre + adjust;
        l = 0;
        if( centre < n + adjust )
            l = n - centre + adjust;
        r = n + centre;

        if( ! ( *GetValidator())( GetRunAt( r )))
            invalid = r;

        //mask windows that contain invalid positions with negative entropy
        if( l <= invalid && invalid <= r )
            SetEntropyAt( n, -1.0f );
        else
            SetEntropyAt( n, window.GetEntropy());
    }

    //reset entropy for the half of positions of the last window
    for( ; n < GetRunLength(); n++ )
        SetEntropyAt( n, -1.0f );
}

// -------------------------------------------------------------------------
// FindExtentAt: find the largest low/high-entropy extent starting from the
// given position
//
void SEGAbstract::FindExtentAt( size_t n, size_t left, size_t right, size_t* locleft, size_t* locright )
{
    if( GetRunLength() <= left || GetRunLength() <= right || right < left )
        throw MYRUNTIME_ERROR( "SEGAbstract::FindExtentAt: Invalid extent boundaries." );

    if( locleft == NULL || locright == NULL )
        return;

    float   entro;
    size_t  p;

    for( p = n; p >= left; p-- ) {
        entro = GetEntropyAt( p );
        if( entro < 0.0f )
            break;

        if( GetInverted()?
            entro < GetLowEntropy() :
            GetHighEntropy() < entro )
            break;

        *locleft = p;

        if( p == 0 )
            break;
    }

    for( p = n; p <= right; p++ ) {
        entro = GetEntropyAt( p );
        if( entro < 0.0f )
            break;
        if( GetInverted()?
            entro < GetLowEntropy() :
            GetHighEntropy() < entro )
            break;
        *locright = p;
    }
}

// -------------------------------------------------------------------------
// OptimizeExtent: find optimal extent within the run segment bounded by
// left and right
//
void SEGAbstract::OptimizeExtent( size_t* locleft, size_t* locright )
{
    if( locleft == NULL || locright == NULL )
        return;

    if( GetRunLength() <= *locleft || GetRunLength() <= *locright || *locright < *locleft )
        throw MYRUNTIME_ERROR( "SEGAbstract::OptimizeExtent: Invalid extent boundaries." );

    float   logprob = 0.0f;
    float   minlogprob = 0.0f;
    size_t  extlength = *locright - *locleft + 1;
    size_t  minlength = ( GetMaxExtentDifference() < extlength  )? extlength - GetMaxExtentDifference(): 1;
    size_t  length;
    size_t  optleft = *locleft;
    size_t  optright = *locright;

    if( GetInverted())
        minlogprob = 100.0f;

    for( length = extlength; length > minlength; length-- )
    {
        Window  window(
                GetVerifier(),
                GetComparator(),
                GetEqComparer(),
                GetRunAt( *locleft ),
                length,
                GetAlphabetSize()
        );

        for( size_t l = *locleft; l + length - 1 <= *locright; l++ )
        {
            window.SafeShift();
            logprob = window.LogProbability();
            if( GetInverted())
                logprob = -logprob;
            if( logprob < minlogprob ) {
                minlogprob = logprob;
                optleft = l;
                optright = l + length - 1;
            }
        }
    }

    *locleft = optleft;
    *locright = optright;
}

// =========================================================================
// PrintSequence: print formatted sequence 
//
void SEGAbstract::PrintSequence( FILE* fp, TResAcceder getres, size_t width )
{
    if( fp == NULL || getres == NULL )
        return;

    size_t  n = 0;
    size_t  p = 0;
    char*   buffer = (char*)malloc( sizeof(char) * ( ((width<GetRunLength())? width: GetRunLength()) + 1 ) );

    if( buffer == NULL )
        throw MYRUNTIME_ERROR( "SEGAbstract::PrintSequence: Not enough memory." );

    for( n = 0; n < GetRunLength(); n += p )
    {
        for( p = 0; p < width && n + p < GetRunLength(); p++ )
            buffer[p] = DehashCode((unsigned char)(*getres)( GetRunAt(n + p)));

        buffer[p] = 0;
        fprintf( fp, "%s%s", buffer, NL );
    }

    fprintf( fp, "%s", NL );
    free( buffer );
}

// -------------------------------------------------------------------------
// PrintSeggedSequence: print formatted sequence with segments found by
//     running the algorithm and masked with Xs
//
void SEGAbstract::PrintSeggedSequence( FILE* fp, TResAcceder getres, size_t width )
{
    if( fp == NULL || getres == NULL )
        return;

    myruntime_error mre;

    char*   residues = (char*)malloc( sizeof(char) * ( GetRunLength() + 1 ));
    char*   buffer = (char*)malloc( sizeof(char) * ( ((width<GetRunLength())? width: GetRunLength()) + 1 ) );

    try{
        if( residues == NULL || buffer == NULL )
            throw MYRUNTIME_ERROR( "SEGAbstract::PrintSeggedSequence: Not enough memory." );

        size_t  n = 0;
        size_t  p = 0;
        size_t  left, right;
        const Segments& loc_segments = GetSegments();

        for( n = 0; n < GetRunLength(); n ++ )
            residues[n] = DehashCode((unsigned char)(*getres)( GetRunAt( n )));

        for( n = 0; n < loc_segments.GetSize(); n++ ) {
            left = loc_segments.GetLeftAt( n );
            right = loc_segments.GetRightAt( n );

            for( p = left; p <= right; p++ )
                if(( *getres )( GetRunAt( n )) != GAP )
                    residues[p] = 'x';//DehashCode( X );
        }

        for( n = 0; n < GetRunLength(); n += p ) {
            memcpy( buffer, residues + n, p = (n+width<GetRunLength())? width: GetRunLength()-n );
            buffer[p] = 0;
            fprintf( fp, "%s%s", buffer, NL );
        }

        fprintf( fp, "%s", NL );

    } catch( myexception const & ex ) {
        mre = ex;
    }

    if( residues )
        free( residues );
    if( buffer )
        free( buffer );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// MaskSequence: mask the sequence intervals found by the SEG algorithm with 
// msym;
// msym, a symbol to maks the segged intervals,
// omitmask, positions not to consider
// omit1 and omit2, symbols not to be masked
//
void SEGAbstract::MaskSequence(
        char* sequence,
        size_t length,
        char msym,
        const char* omitmask,
        char omit1,
        char omit2 ) const
{
    if( sequence == NULL || omitmask == NULL || length == 0 )
        return;

    size_t  nn = 0;//sequence iterator
    size_t  n = 0;
    size_t  p = 0;
    size_t  left, right;
    const Segments& loc_segments = GetSegments();


    for( n = 0; n < loc_segments.GetSize(); n++ ) {
        left = loc_segments.GetLeftAt( n );
        right = loc_segments.GetRightAt( n );

        for( ; p < left && nn < length; nn++ ) {
            if( omitmask[nn] )
                continue;
            p++;
        }

        for( ; p <= right && nn < length; nn++ ) {
            if( omitmask[nn] )
                continue;

            if( sequence[nn] != omit1 && sequence[nn] != omit2 )
                sequence[nn] = msym;
            p++;
        }
    }
}

// =========================================================================
// testing:
// Print: print segments found and the vector of computed entropies
//
void SEGAbstract::Print( FILE* fp )
{
    if( fp == NULL )
        return;

    fprintf( fp, "%sSEGAbstract%s Entropies:%s ", NL, NL, NL );
    for( size_t n = 0; n < GetRunLength(); n++ )
        fprintf( fp, " %6.3f", GetEntropyAt( n ));

    fprintf( fp, "%s", NL );
    segments.Print( fp );
}

}//namespace SEG
