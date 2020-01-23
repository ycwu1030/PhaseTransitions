/*
 * @Description  : Tracing the minima of a potential and find Critical Temperature
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-23 22:12:32
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-23 17:11:52
 */
#ifndef TraceMin_H
#define TraceMin_H
#include "VTypes.h"
#include <string>

struct _traceMinimum_rval
{
    VVD X;
    VD T;
    VVD dXdT;

    VD overX;
    double overT;
};

struct precision_control
{
    double dtabsMax = 20.0; // Control the maximum step size in t, relative to the starting step size;
    double dtfracMax = 0.25;  // Control the maximum step size in t, relative to the temperature of that point;
    double dtmin = 1e-3; // Control the minimum step size in t, relative to the starting step size;
    double deltaX_target = 1e-3; // Control the error in X we can accept;
    double deltaX_tol = 1.2; // Setting the maximum difference between the guess minimum and the truth minimum;
    double minratio_rel = 1e-2; // Setting the ratio between minimum eigenvalue to maximum eigenvalue, such that we treat the minimum one as zero. Relative to the initial ratio. (Assuming that we start at a point not close to saddle point)
    double dtstart = 1e-3; // Control the starting step size in t, relative to the difference between the lowest temperature and the highest temperature
    double tjump = 1e-3; // Control the jump step size in t, relative to the difference between the lowest temperature and the highest temperature. The `jump` means when we found one phase is end (become saddle/maximum point), we jump a step in t (decreasing or increasing depends on we were down-tracing or up-tracing) and try to start another phase
} pre_control;


_traceMinimum_rval traceMinimum(ScalarFunction f, dScalarFunction df_dx, dScalarFunction d2f_dxdt, HM d2f_dx2, VD x0, double t0, double tstop, double dtstart/*, double deltaX_target, double dtabsMax=20.0, double dtfracMax=0.25, double dtmin=1e-3, double deltaX_tol=1.2, double minratio=1e-4*/);

MP traceMultiMin(ScalarFunction f, dScalarFunction df_dx, dScalarFunction d2f_dxdt, HM d2f_dx2, VT points, double tLow, double tHigh, /*double deltaX_target, double dtstart=1e-3, double tjump=1e-3,*/ forbidCrit fc = AllPass);

VVD findApproxLocalMin(ScalarFunction f, VD x0, VD x1, double ti, int n = 100, double edge = 0.05);

void removeRedundantPhases(ScalarFunction f, dScalarFunction df_dx, MP &phases, double xeps=1e-5, double diftol=1e-2);

void _removeRedundantPhase(MP &phases,int index_removed, int index_phase);

struct TransCritical
{
    double Tcrit;
    double Tnuc;
    VD low_vev;
    VD high_vev;
    int low_phase;
    int high_phase;
    int trantype;
};
TransCritical secondOrderTrans(Phase phase1, Phase phase2, std::string = "Tcrit");
typedef std::vector<TransCritical> VTC;
VTC findCriticalTemperatures(MP phases, ScalarFunction f);

#endif