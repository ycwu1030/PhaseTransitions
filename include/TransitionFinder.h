/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-23 18:00:46
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-26 23:50:05
 */
#ifndef TransitionFinder_H
#define TransitionFinder_H

#include "Phases.h"
#include "TraceMin.h"
#include "VTypes.h"

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
#endif