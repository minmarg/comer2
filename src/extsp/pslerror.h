/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_pslerror__
#define __extspsl_pslerror__

namespace extspsl {

#define PRIVERROR(MSG) priverror(MSG,__FILE__,__LINE__,__func__)

#define PSL_OK          (   0 )
#define PSL_SUCCESS     (   0 )
#define PSL_ERR_DOMAIN  ( 111 )
#define PSL_ERR_NOPROG  ( 113 )
#define PSL_MAXITERATS  ( 117 )
#define PSL_ERR_DIM     ( 209 )
#define PSL_ERR_MEMORY  ( 211 )
#define PSL_ERR_ADDRESS ( 213 )
#define PSL_ERR_ILLEGAL ( 217 )
#define PSL_ERR_INVALID ( 311 )

#define PSL_ERR_ROUNDOFF    ( 401 )
#define PSL_ERR_UNDERFLOW   ( 403 )
#define PSL_ERR_OVERFLOW    ( 405 )

//gammainc-specific
#define PSL_ERR_GI_XRANGE   ( 501 )
#define PSL_ERR_GI_SMAXIT   ( 503 )
#define PSL_ERR_GI_AMAXIT   ( 505 )
#define PSL_ERR_GI_FMAXIT   ( 507 )

//function to report an error
void priverror( const char*, 
     const char* file, unsigned int line, const char* func );

inline const char* TranslatePSLError( int code )
{
    switch( code ) {
        case PSL_SUCCESS:       return "Converged";
        case PSL_ERR_DOMAIN:    return "Domain error";
        case PSL_ERR_NOPROG:    return "Stopped: No progress towards solution";
        case PSL_MAXITERATS:    return "Maximum number of iterations reached";
        case PSL_ERR_DIM:       return "Vector/matrix dimensions error";
        case PSL_ERR_MEMORY:    return "Not enough memory";
        case PSL_ERR_ADDRESS:   return "Memory access error";
        case PSL_ERR_ILLEGAL:   return "Illegal instructions";
        case PSL_ERR_INVALID:   return "Invalid data";
        case PSL_ERR_ROUNDOFF:  return "Roundoff error";
        case PSL_ERR_UNDERFLOW: return "Underflow";
        case PSL_ERR_OVERFLOW:  return "Overflow";
        //gammainc-specific
        case PSL_ERR_GI_XRANGE: return "Incomplete gamma calculation by series failed: x>>a exceeds range";
        case PSL_ERR_GI_SMAXIT: return "Incomplete gamma calculation by series failed: Maximum number of iterations reached";
        case PSL_ERR_GI_AMAXIT: return "Asymptotic incomplete gamma calculation failed: Maximum number of iterations reached";
        case PSL_ERR_GI_FMAXIT: return "Incomplete gamma calculation by contd. fraction failed: Maximum number of iterations reached";
    }
    return "Unknown";
}

}//namespace extspsl

#endif//__extspsl_pslerror__
