/*
 * @Description  : Tracing the minima of a potential and find Critical Temperature
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-23 22:12:32
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-03 09:34:06
 */
#ifndef TraceMin_H
#define TraceMin_H
#include "VTypes.h"

struct _traceMinimum_rval
{
    VVD X;
    VD T;
    VVD dXdT;

    VD overX;
    double overT;
};

_traceMinimum_rval traceMinimum(ScalarFunction f, dScalarFunction df_dx, dScalarFunction d2f_dxdt, HM d2f_dx2, VD x0, double t0, double tstop, double dtstart, double deltaX_target, double dtabsMax=20.0, double dtfracMax=0.25, double dtmin=1e-3, double deltaX_tol=1.2, double minratio=1e-2);

MP traceMultiMin(ScalarFunction f, dScalarFunction df_dx, dScalarFunction d2f_dxdt, HM d2f_dx2, VT points, double tLow, double tHigh, double deltaX_target, double dtstart=1e-3, double tjump=1e-3, forbidCrit fc = AllPass);

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
TransCritical secondOrderTrans(Phase phase1, Phase phase2, char *str = "Tcrit");
typedef std::vector<TransCritical> VTC;
VTC findCriticalTemperatures(MP phases, ScalarFunction f);

#endif