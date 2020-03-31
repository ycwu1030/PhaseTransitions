#ifndef TransitionFinder_H
#define TransitionFinder_H

#include "Phases.h"
#include "TraceMin.h"
#include "VTypes.h"
#include "PathDeformation.h"

#define NTempType 2
typedef enum {
    TCRIT = 0,
    TNUCL = 1
} TempType;

struct TransCritical
{
    // double Tcrit;
    // double Tnuc;
    double T[NTempType];
    VD low_vev;
    VD high_vev;
    int low_phase;
    int high_phase;
    int tranorder;
    TempType trantype;
};
TransCritical secondOrderTrans(Phase phase1, Phase phase2, TempType ttype = TCRIT);
typedef std::vector<TransCritical> VTC;
VTC findCriticalTemperatures(MP phases, ScalarFunction f);

struct TransNucleation
{
    VD R;
    VD Phi_1D;
    VD dPhi_1D;
    VVD Phi;
    double action;
    double fRatio;
    // std::vector<VVD> saved_steps;
};
TransNucleation fullTunneling(VVD pts_init, ScalarFunction V_in, dScalarFunction dV_in, double T, int maxiter=20, double fixEndCutoff=0.03, bool save_all_steps=false,int V_spline_samples=100);


#endif