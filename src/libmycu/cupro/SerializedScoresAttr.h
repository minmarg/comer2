/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedScoresAttr_h__
#define __SerializedScoresAttr_h__

//attributes of serialized scores
struct SerializedScoresAttr {
    SerializedScoresAttr( 
            int ntotents, int szalloc, char nplvs, char nenos,
            int card, int ntbls, int nelems,
            float eth1 = 0.0f, float eth2 = 0.0f )
    :   ntotents_(ntotents), szalloc_(szalloc),
        nplvs_(nplvs), nenos_(nenos),
        card_(card), ntbls_(ntbls), nelems_(nelems),
        eth1_(eth1), eth2_(eth2)
    {}
    SerializedScoresAttr( const SerializedScoresAttr& attr )
    :   ntotents_(attr.ntotents_), szalloc_(attr.szalloc_),
        nplvs_(attr.nplvs_), nenos_(attr.nenos_),
        card_(attr.card_), ntbls_(attr.ntbls_), nelems_(attr.nelems_),
        eth1_(attr.eth1_), eth2_(attr.eth2_)
    {}
    SerializedScoresAttr()
    :   ntotents_(0), szalloc_(0),
        nplvs_(0), nenos_(0),
        card_(0), ntbls_(0), nelems_(0),
        eth1_(0.0f), eth2_(0.0f)
    {}
    int ntotents_;//total number of entries in the buffer, included cardinalities, etc.
    int szalloc_;//size (bytes) allocated for scores
    char nplvs_;//number of probability levels
    char nenos_;//number of levels for eff. number of observations
    int card_;//cardinality determinator, # rows of a square score table
    int ntbls_;//number of tables per a pair of probability level values
    int nelems_;//number of entries in one table
    float eth1_;//the first threshold of the effective number of observations [OPTIONAL]
    float eth2_;//the second threshold of the effective number of observations [OPTIONAL]
};

#endif//__SerializedScoresAttr_h__
