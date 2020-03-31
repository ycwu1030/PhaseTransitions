#include "TransitionFinder.h"
#include "Tunneling1D.h"
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

TransNucleation fullTunneling(VVD pts_init, ScalarFunction V_in, dScalarFunction dV_in, double T, int maxiter, double fixEndCutoff, bool save_all_steps,int V_spline_samples)
{
    vector<VVD> saved_steps;
    VVD pts = pts_init;
    VD R;
    VD Phi_1D,phi;
    VD dPhi_1D,dphi;
    double Rerr;
    Deformation_Status deform_info;
    for (size_t num_iter = 0; num_iter < maxiter; num_iter++)
    {
        // cout<<"Starting tunneling step-"<<num_iter+1<<endl;
        // * 1. Interpolate the path by spline
        SplinePath path(pts,V_in,T,V_spline_samples,true);
        // cout<<"\tGot the spline path"<<endl;
        // for (double s = 0; s <= path.GetDistance(); s+= 0.01*path.GetDistance())
        // {
        //     cout<<"\t"<<s<<"\t"<<path.pts_at_dist(s)<<"\t"<<path.V({s},nullptr)<<"\t"<<path.dV({s},nullptr)[0]<<endl;
        // }
        // double s = path.GetDistance();
        // cout<<"\t"<<s<<"\t"<<path.pts_at_dist(s)<<"\t"<<path.V({s},nullptr)<<"\t"<<path.dV({s},nullptr)[0]<<endl;
        
        // * 2. Peform 1-D tunneling along the above path;
        Tunneling1D tunnel1D(0.0,path.GetDistance(),path.V,path.dV,path.d2V);
        // cout<<"\tTry to find 1D profile"<<endl;
        tie(R,Phi_1D,dPhi_1D,Rerr) = tunnel1D.findProfile();
        phi = Phi_1D;
        dphi = dPhi_1D;
        tie(phi,dphi) = tunnel1D.evenlySpacedPhi(phi,dphi,phi.size(),1,false);
        // cout<<"\tGot 1D tunneling results"<<endl;

        dphi.front() = 0;
        dphi.back() = 0;

        // * 3. Deform the path
        pts = path.pts_at_dist(phi);
        Deformation_Spline deform(pts,dphi,dV_in);
        deform_info = deform.deformPath();
        pts = deform.GetPhi();
        // cout<<"\tPath deformed"<<endl;
        // temporally ignore saving the steps

        // * 4. Check convergence
        // * If deformation converged after one step, then we can assume that the path is good.
        // * We don't allow more steps, since during deformPath() we didn't re-calculate the 1D tunneling profile. With more steps, the convergence status should be recalculated.
        if (deform_info == DF_CONVERGED && deform.GetSteps() < 2)
        {
            break;
        }
    }
    Deformation_Spline deform(pts,dphi,dV_in);
    VVD F, Gradient;
    tie(F,Gradient) = deform.forces();
    VD F_mag = pow(F*F,0.5);
    VD Grad_mag = pow(Gradient*Gradient,0.5);
    double F_mag_max = *max_element(F_mag.begin(),F_mag.end());
    double Grad_mag_max = *max_element(Grad_mag.begin(),Grad_mag.end());
    double fRatio = F_mag_max/Grad_mag_max;

    SplinePath path_final(pts,V_in,T,V_spline_samples,true);
    VVD Phi = path_final.pts_at_dist(Phi_1D);

    Tunneling1D tunnel1D_final(0.0,path_final.GetDistance(),path_final.V,path_final.dV,path_final.d2V);
    double action = tunnel1D_final.findAction(R,Phi_1D,dPhi_1D);

    TransNucleation res = {R,Phi_1D,dPhi_1D,Phi,action,fRatio};
    return res;
}