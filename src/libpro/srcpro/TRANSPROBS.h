/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __TRANSPROBS_h__
#define __TRANSPROBS_h__


#define MAIN_TRNPRB_LOGFCT 0.4f

//enums
//profile transition states
enum TPTRANS {
    P_MM, P_MI, P_MD,
    P_IM, P_II, P_ID,
    P_DM, P_DI, P_DD,
    P_NSTATES
};
extern const char* gTPTRANS_NAMES[P_NSTATES];

//profile states
enum TPSTATES {
    PS_M,
    PS_I,
    PS_D,
    PS_NSTATES
};

extern const char* gTPSTATES_NAMES[PS_NSTATES];

//number of transitions per state
extern const int gTPTRANS_NTPS;// = 3;

//Dirichlet distribution parameters, alphas, for transition priors
extern float gTPTRANS_ALPHAS[P_NSTATES];

// _________________________________________________________________________
// CLASS _TRANSPROBS
//
class _TRANSPROBS
{
public:
    enum TPriorClass {
        Dirichlet   //Dirichlet prior
    };
    _TRANSPROBS();

    TPriorClass GetClass() const { return class_; }
    void SetClass( TPriorClass value );

    float GetEffNoSequences( int mid ) const;
    void SetEffNoSequences( float nm, float ni, float nd );

    void PME( const float (*obsfreqs)[P_NSTATES]);//posteriorize

    const float (*GetPriors() const)[P_NSTATES] { return &priors_; }
    const float (*GetPMEstimators() const)[P_NSTATES] { return &pmestimators_; }
    const float (*GetLogPMEstimators() const)[P_NSTATES] { return &logpmestimators_; }

protected:
    void InitOnce();
    void Initialize();
    void DirPME( const float(*obsfreqs)[P_NSTATES]);//posteriorize with Dirichlet priors

private:
    TPriorClass class_; //class of priors
    float effnos_[PS_NSTATES];//effective number of sequences for each state
    float sumalphas_[PS_NSTATES];//sum of alpha parameters
    float priors_[P_NSTATES];//prior transition probabilities
    float pmestimators_[P_NSTATES];//transition probability estimators
    float logpmestimators_[P_NSTATES];//log values of estimators
};
extern _TRANSPROBS TRANSPROBS;
// //

#endif//__TRANSPROBS_h__
