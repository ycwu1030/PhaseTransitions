#include "TransitionFinder.h"
#include <tuple>
#include <functional>
#include <algorithm>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>

using namespace std;


bool TransCompTc(TransCritical t1, TransCritical t2)
{
    return t1.T[t1.trantype] > t2.T[t2.trantype];
}
struct minipara
{
    ScalarFunction f;
    Phase ph1;
    Phase ph2;
};
double Vcrit_wrap(double T, void *param)
{
    minipara *par = (minipara*)param;
    double v1 = (par->f)((par->ph1).valAt(T),&T);
    double v2 = (par->f)((par->ph2).valAt(T),&T);
    return v1-v2;
}

/*
 * @description: Find all possible critical temperature
 * @param
 *      phases:
 *          The phase map, from traceMultimin (and remove redundant phase)
 *      f:
 *          The potential
 * @return:
 *      A list of critical point
 *          Tcrit,Tnuc:
 *              The critical temperature
 *                  either Tcrit or Tnuc, the other one will be set to -1;
 *          low_vev, high_vev:
 *              The filed value at critical point for low and high temperature phase
 *          low_phase, high_phase:
 *              The key identifying the phase at low and high temperature
 *          trantype:
 *              1: first order
 *              2: second order
 */
VTC findCriticalTemperatures(MP phases, ScalarFunction f)
{
    VTC transitions;
    MP::iterator iter1,iter2;
    SI::iterator iterSI;
    int index1,index2;
    double Tcrit;
    for (iter1 = phases.begin(); iter1 != phases.end(); iter1++)
    {
        for (iter2 = phases.begin(); iter2 != phases.end(); iter2++)
        {
            if (iter1 == iter2)
            {
                continue;
            }
            index1 = iter1->first;
            index2 = iter2->first;
            Phase phase1 = phases.at(index1);
            Phase phase2 = phases.at(index2);
            double tmax = min(phase1.GetTmax(),phase2.GetTmax());
            double tmin = max(phase1.GetTmin(),phase2.GetTmin());
            if (tmin >= tmax)
            {
                // * No overlap, check whether it is second-order PT
                SI ph1SILow = phase1.GetLowTrans();
                if (ph1SILow.find(index2) != ph1SILow.end())
                {
                    transitions.push_back(secondOrderTrans(phase1,phase2,TCRIT));
                }
                continue;
            }

            // For 1-st order PT, just consider phase1->phase2
            if (f(phase1.valAt(tmin),&tmin)-f(phase2.valAt(tmin),&tmin) < 0)
            {
                continue;
            }
            if (f(phase1.valAt(tmax),&tmax)-f(phase2.valAt(tmax),&tmax) > 0)
            {
                continue;
            }
            function<double() > findTc = [&](){
                int status;
                int iter = 0;
                int max_iter = 200;
                const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
                gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
                minipara params = {f, phase1, phase2};
                gsl_function F;
                F.function = &Vcrit_wrap;
                F.params = &params;
                double t_low = tmin;
                double t_high = tmax;
                double t_crit;
                gsl_root_fsolver_set(s, &F, t_low, t_high);

                do
                {
                    iter++;
                    status = gsl_root_fsolver_iterate(s);
                    t_crit = gsl_root_fsolver_root(s);
                    t_low = gsl_root_fsolver_x_lower(s);
                    t_high = gsl_root_fsolver_x_upper(s);
                    status = gsl_root_test_interval (t_low,t_high,0.001,0.001);
                } while (status == GSL_CONTINUE && iter < max_iter);
                
                t_crit = gsl_root_fsolver_root(s);
                gsl_root_fsolver_free(s);
                return t_crit;
            };
            Tcrit = findTc();
            transitions.push_back({{Tcrit,-1},phase2.valAt(Tcrit),phase1.valAt(Tcrit),index2,index1,1,TCRIT});
        }
    }
    sort(transitions.begin(),transitions.end(),TransCompTc);
    return transitions;
}

TransCritical secondOrderTrans(Phase high_phase, Phase low_phase, TempType Ttype)
{
    TransCritical res;
    res.T[0] = -1;
    res.T[1] = -1;
    res.T[Ttype] = 0.5*(high_phase.GetTmin()+low_phase.GetTmax());
    res.low_vev = low_phase.GetXatTmax();
    res.high_vev = high_phase.GetXatTmin();
    res.low_phase = low_phase.GetKey();
    res.high_phase = high_phase.GetKey();
    res.tranorder = 2;
    res.trantype = Ttype;
    return res;
}