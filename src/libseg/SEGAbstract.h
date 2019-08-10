/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SEGAbstract_h__
#define __SEGAbstract_h__

#include <stdio.h>
#include "liblib/BinarySearchStructure.h"
#include "Segments.h"

namespace SEG {

//count type
typedef size_t  tcount;

// Validator function:
// functions of this type return true if an element is 
// valid for analysis, false otherwise
typedef bool ( *TValidator )( const void* key );
typedef bool ( *TVerifier )( const void* key );

// residue access function; required only for a subclass of problems
typedef size_t ( *TResAcceder )( const void* key );

// comparison function for equality only; returns true if two 
// keys are equal, false otherwise
typedef bool   ( *TEqComparer )( const void* key1, const void* key2 );

// =========================================================================
// Class LinearSearchStructure
//
class LinearSearchStructure: public BinarySearchStructure
{
public:
    LinearSearchStructure( TComparator, TEqComparer, size_t size );
    virtual ~LinearSearchStructure();

    virtual bool Find( const void* key, int* loc = NULL ) const;

private:
    TEqComparer eqcomparer;//comparison function for equality only
};

// =========================================================================
// Class BinaryModifier
//
class BinaryModifier: public BinarySearchStructure
{
public:
    BinaryModifier( TComparator, size_t size );
    virtual ~BinaryModifier();

    bool    IncValue( tcount value );
    bool    DecValue( tcount value );

private:
};

// =========================================================================
// Class Window
//
class Window
{
public:
    Window(
        TVerifier   verffunc,
        TComparator comparator,
        TEqComparer eqcomparer,
        const void* address,
        size_t      len,
        size_t      szalpha = 0
    );

    ~Window();

    size_t      GetLength() const           { return length; }
    size_t      GetAlphabetSize() const     { return szalphabet; }
    float       GetEntropy() const          { return entropy; }

    TVerifier           GetVerifier() const { return verifier; }
    const void*         GetWindowAt( size_t pos ) const;
    tcount              GetCompositionAt( size_t n ) const;
    tcount              GetStateAt( size_t n ) const;

    bool        Uninitialized() const       { return GetEntropy() < 0.0f; }

    void        Initialize();
    void        SafeShift();

    float       LogProbability() const;

    void        Print( FILE* );//testing

protected:
    const void*             GetWindow()                 { return window; }
    SimpleVector*           GetComposition()            { return composition; }
    BinarySearchStructure*  GetCompAddresses()          { return compaddresses; }
    BinarySearchStructure*  GetState()                  { return state; }

    size_t                  GetCompositionSize() const;
    size_t                  GetCompAddressesSize() const;
    size_t                  GetStateSize() const;

    bool        RefreshComposition( const void* value, size_t* address = NULL );

    void        IncCompositionAt( size_t n );//increment composition value by 1
    void        DecCompositionAt( size_t n );//decrement composition value by 1
    void        SetCompositionAt( size_t n, tcount value );
    void        InsertCompositionAt( size_t n, tcount value );

    bool        FindCompAddress( const void* value, size_t* address );
    bool        PushCompAddress( const void* value, size_t* address );

    void        PushStateValue( tcount value );
    void        IncStateValue( tcount value );
    void        DecStateValue( tcount value );

    void        Shift();

    void        SetEntropy( float value ) { entropy = value; }

    void        StateOn();
    void        CompositionOn();
    void        ComputeEntropy();

    float       LogMultinomial() const;
    float       LogNoCompositions( size_t = 0 ) const;

private:
    const TVerifier         verifier;           //verifier of window positions
    const void*             window;             //beginning of actual window
    const size_t            length;             //length of window
    const size_t            szalphabet;         //size of alphabet which can be given or not
    float                   entropy;            //entropy of the window
    SimpleVector*           composition;        //composition vector of window
    BinarySearchStructure*  compaddresses;      //addresses of composition vector
    BinaryModifier*         state;              //state vector (sorted composition) of window
};

// =========================================================================
// _________________________________________________________________________
// Class SEGAbstract
//
class SEGAbstract
{
public:
    SEGAbstract(
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
        size_t  szalpha = 0
    );

    SEGAbstract(
        TValidator  valdfunc,
        TVerifier   verffunc,
        TComparator compfunc,
        TEqComparer eqcomfunc,
        size_t  winlen,
        float   lowent,
        float   highent,
        size_t  maxdiff,
        size_t  szalpha = 0
    );

    virtual ~SEGAbstract();

    void        SetRunAddress( void* address, size_t runlen );

    void        Run();

    void        Initialize();
    void        FindSegments();

    float       Entropy() const;//compute entropy of segment of entire run length
    float       LogProbability() const;//compute log-probability of segment of entire run length

    TValidator          GetValidator() const    { return validator; }
    TVerifier           GetVerifier() const     { return verifier; }
    TComparator         GetComparator() const   { return comparator; }
    TEqComparer         GetEqComparer() const   { return eqcomparer; }
    const void*         GetRunAddress() const   { return runaddress; }
    size_t              GetRunLength() const    { return runlength; }
    size_t              GetWinLength() const    { return winlength; }
    size_t              GetAlphabetSize() const { return szalphabet; }
    float               GetLowEntropy() const   { return lowentropy; }
    float               GetHighEntropy() const  { return highentropy; }
    size_t              GetMaxExtentDifference() const { return maxdiffextent; }

    void        Print( FILE* );//testing

    void        SetLowCSearch() { inverted = false; }
    void        SetHighCSearch() { inverted = true; }

    void        MaskSequence( char*, size_t length, char msym, const char* omitmask, char, char ) const;

protected:
    const void* GetRunAt( size_t pos ) const;

    float       GetEntropyAt( size_t pos ) const;
    void        SetEntropyAt( size_t pos, float value );

    const Segments& GetSegments() const { return segments; }

    bool        GetInverted() const { return inverted; }

    void        Allocate( size_t newlen );
    void        Destroy();

    void        ComputeEntropy();
    void        RecurrentSegFind( size_t left, size_t right );
    void        FindExtentAt( size_t pos, size_t left, size_t right, size_t* locleft, size_t* locright );
    void        OptimizeExtent( size_t* locleft, size_t* locright );

    //-- semi-abstract methods
    void        PrintSequence( FILE*, TResAcceder, size_t width = scszDefSegSeqPrintWidth );
    void        PrintSeggedSequence( FILE*, TResAcceder, size_t width = scszDefSegSeqPrintWidth );

protected:
    static const size_t scszDefSegSeqPrintWidth;//default print width
private:
    static const size_t scszDefNoSegments;//default number of segments

private:
    const TValidator    validator;  //validator function of run elements
    const TVerifier     verifier;   //verifier function for run elements
    const TComparator   comparator; //comparer of run elements
    const TEqComparer   eqcomparer; //comparison function for equality only
    const void*         runaddress; //address of run
    const size_t        runlength;  //length of run
    const size_t        winlength;  //length of window to use
    const size_t        szalphabet; //size of alphabet a run comprises
    const float         lowentropy; //low entropy threshold
    const float         highentropy;//high entropy threshold
    const size_t        maxdiffextent;//maximum difference between extent (segment) length and minimum extent length

    bool                inverted;   //if true, high complexity segment search is in effect

    float*              entropies;  //entropies computed and attributed to window centres
    Segments            segments;   //storage for found segments
};

////////////////////////////////////////////////////////////////////////////
// Class Window INLINES
//
// GetWindowAt: get the address of window element at the given position
//
inline
const void* Window::GetWindowAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !window || length <= pos )
        throw MYRUNTIME_ERROR( "Window::GetWindowAt: Memory access error." );
#endif
    return ( const char* )window + pos * sizeof(void*);
}

// -------------------------------------------------------------------------
// GetCompositionSize: get the size of composition vector
//
inline
size_t Window::GetCompositionSize() const
{
#ifdef __DEBUG__
    if( !composition )
        throw MYRUNTIME_ERROR( "Window::GetCompositionSize: Memory access error." );
#endif
    return composition->GetSize();
}

// -------------------------------------------------------------------------
// GetCompAddressesSize: get the size of vector of composition addresses
//
inline
size_t Window::GetCompAddressesSize() const
{
#ifdef __DEBUG__
    if( !compaddresses )
        throw MYRUNTIME_ERROR( "Window::GetCompAddressesSize: Memory access error." );
#endif
    return compaddresses->GetSize();
}

// -------------------------------------------------------------------------
// GetStateSize: get the size of state vector
//
inline
size_t Window::GetStateSize() const
{
#ifdef __DEBUG__
    if( !state )
        throw MYRUNTIME_ERROR( "Window::GetStateSize: Memory access error." );
#endif
    return state->GetSize();
}

// -------------------------------------------------------------------------
// GetCompositionAt: get composition value of element
//
inline
tcount Window::GetCompositionAt( size_t n ) const
{
#ifdef __DEBUG__
    //number of composition values can be greater than the length of window
    if( !composition )
        throw MYRUNTIME_ERROR( "Window::GetCompositionAt: Memory access error." );
#endif
    return (tcount)composition->GetValueAt(n);
}

// -------------------------------------------------------------------------
// InsertCompositionAt: insert composition value of alphabet element
//
inline
void Window::InsertCompositionAt( size_t n, tcount value )
{
#ifdef __DEBUG__
    if( !composition )
        throw MYRUNTIME_ERROR( "Window::InsertCompositionAt: Memory access error." );
#endif
    composition->InsertValueAt( n, (void*)(value));
}

// -------------------------------------------------------------------------
// SetCompositionAt: set composition value of alphabet element
//
inline
void Window::SetCompositionAt( size_t n, tcount value )
{
#ifdef __DEBUG__
    if( !composition )
        throw MYRUNTIME_ERROR( "Window::SetCompositionAt: Memory access error." );
#endif
    composition->SetValueAt( n, (void*)(value));
}

// -------------------------------------------------------------------------
// IncCompositionAt: increment composition value by 1
//
inline
void Window::IncCompositionAt( size_t n )
{
    SetCompositionAt( n, GetCompositionAt(n)+1 );
}

// -------------------------------------------------------------------------
// DecCompositionAt: decrement composition value by 1
//
inline
void Window::DecCompositionAt( size_t n )
{
    tcount curval = GetCompositionAt(n);
    if( curval < 1 )
        return;
    SetCompositionAt( n, curval-1 );
}

// -------------------------------------------------------------------------
// FindCompAddress: find individual value of alphabet; a position 
// corresponds to the address; return true if found, false otherwise
//
inline
bool Window::FindCompAddress( const void* value, size_t* address )
{
#ifdef __DEBUG__
    if( !compaddresses )
        throw MYRUNTIME_ERROR( "Window::FindCompAddress: Memory access error." );
#endif
    int     location = -1;
    bool    found;

    found = compaddresses->Find( value, &location );

    if( found && location < 0 )
        throw MYRUNTIME_ERROR( "Window::FindCompAddress: "
                               "Failed to determine composition value address." );
    if( found && address )
        *address = (size_t)location;

    return found;
}

// -------------------------------------------------------------------------
// PushCompAddress: insert individual value of alphabet; a position will
// determine the address; return true if the insert was successful, false
// otherwise, but the address will be refreshed anyway
//
inline
bool Window::PushCompAddress( const void* value, size_t* address )
{
#ifdef __DEBUG__
    if( !compaddresses )
        throw MYRUNTIME_ERROR( "Window::PushCompAddress: Memory access error." );
#endif
    int     location = -1;
    bool    inserted;

    inserted = compaddresses->Push( value, &location );

    if( location < 0 )
        throw MYRUNTIME_ERROR( "Window::PushCompAddress: "
                               "Unable to determine composition value address." );
    if( address )
        *address = (size_t)location;

    return inserted;
}

// -------------------------------------------------------------------------
// RefreshComposition: refresh compositional information given the 
// appearance of this new value; return the flag of whether the value is 
// new to the composition
//
inline
bool Window::RefreshComposition( const void* value, size_t* saveaddr )
{
    if( ! (*GetVerifier())( value ))
        return false;

    size_t  address = 0;
    bool    appended = PushCompAddress( value, &address );
    //compositions will be stored continuously as this is provided by
    //addresses' BinarySearchStructure
    if( appended )
        InsertCompositionAt( address, 0 );

    IncCompositionAt( address );

    if( saveaddr )
        *saveaddr = address;

    return appended;
}

// -------------------------------------------------------------------------
// GetStateAt: gets the state value of alphabet element
//
inline
tcount Window::GetStateAt( size_t n ) const
{
#ifdef __DEBUG__
    //number of state values can be greater than the length of window
    if( !state )
        throw MYRUNTIME_ERROR( "Window::GetStateAt: Memory access error." );
#endif
    return ( tcount )state->GetValueAt( n );
}

// -------------------------------------------------------------------------
// PushStateValue: insert state value 
//
inline
void Window::PushStateValue( tcount value )
{
#ifdef __DEBUG__
    if( !state )
        throw MYRUNTIME_ERROR( "Window::PushStateValue: Memory access error." );
#endif
    state->Push((void*)value );
}

// -------------------------------------------------------------------------
// IncStateValue: increase state value by 1
//
inline
void Window::IncStateValue( tcount value )
{
#ifdef __DEBUG__
    if( !state )
        throw MYRUNTIME_ERROR( "Window::IncStateValue: Memory access error." );
#endif
    if( !state->IncValue( value ))
        throw MYRUNTIME_ERROR( "Window::IncStateValue: State value not found." );
}

// -------------------------------------------------------------------------
// DecStateValue: decrease state value by 1
//
inline
void Window::DecStateValue( tcount value )
{
#ifdef __DEBUG__
    if( !state )
        throw MYRUNTIME_ERROR( "Window::DecStateValue: Memory access error." );
#endif
    if( !state->DecValue( value ))
        throw MYRUNTIME_ERROR( "Window::DecStateValue: State value not found." );
}


////////////////////////////////////////////////////////////////////////////
// Class SEGAbstract INLINES
//
// GetRunAt: return the address of element of the run at the given position
//
inline
const void* SEGAbstract::GetRunAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !runaddress || runlength <= pos )
        throw MYRUNTIME_ERROR( "SEGAbstract::GetRunAt: Memory access error." );
#endif
    return (const char*)runaddress + pos * sizeof(void*);
}

// -------------------------------------------------------------------------
// GetEntropyAt: get entropy at the specified run position
//
inline
float SEGAbstract::GetEntropyAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !entropies || runlength <= pos )
        throw MYRUNTIME_ERROR( "SEGAbstract::GetEntropyAt: Memory access error." );
#endif
    return entropies[pos];
}

// -------------------------------------------------------------------------
// SetEntropyAt: set entropy at the given run position
//
inline
void SEGAbstract::SetEntropyAt( size_t pos, float value )
{
#ifdef __DEBUG__
    if( !entropies || runlength <= pos )
        throw MYRUNTIME_ERROR( "SEGAbstract::SetEntropyAt: Memory access error." );
#endif
    entropies[pos] = value;
}

}//namespace SEG

#endif//__SEGAbstract_h__
