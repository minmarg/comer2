/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedCVS2ScoresAttr_h__
#define __SerializedCVS2ScoresAttr_h__

//attributes of serialized scores
struct SerializedCVS2ScoresAttr {
    SerializedCVS2ScoresAttr( const SerializedCVS2ScoresAttr& attr )
    :   ntotents_(attr.ntotents_), szalloc_(attr.szalloc_),
        nenos_(attr.nenos_), ntbls_(attr.ntbls_), card_(attr.card_), shft_(attr.shft_),
        key1first_(attr.key1first_), value1first_(attr.value1first_), 
        key1last_(attr.key1last_), value1last_(attr.value1last_),
        CVS_loKAPPA0_(attr.CVS_loKAPPA0_), CVS_PowerNU0_(attr.CVS_PowerNU0_), 
        CVS_CTERM_(attr.CVS_CTERM_)
    {}
    SerializedCVS2ScoresAttr()
    :   ntotents_(0), szalloc_(0),
        nenos_(0), ntbls_(0), card_(0), shft_(0),
        key1first_(0.0f), value1first_(0.0f), 
        key1last_(0.0f), value1last_(0.0f),
        CVS_loKAPPA0_(0.0f), CVS_PowerNU0_(0.0f), CVS_CTERM_(0.0f)
    {}
    int ntotents_;//total number of entries in the buffer, included nenos, ntbls, etc.
    int szalloc_;//size (bytes) allocated for scores
    char nenos_;//number of levels for eff. number of observations
    int ntbls_;//number of tables per a pair of eno level values
    int card_;//cardinality for scores
    int shft_;//shift to account for negative scores
    float key1first_, value1first_, key1last_, value1last_;//boundary keys and values from a map
    float CVS_loKAPPA0_;//terms obtained from the CVS object of type _TCTXVECT
    float CVS_PowerNU0_;
    float CVS_CTERM_;
};

#endif//__SerializedCVS2ScoresAttr_h__
