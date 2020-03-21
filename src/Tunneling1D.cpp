#include "Tunneling1D.h"
#include <cmath>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace boost::math;

typedef double (*root_func)(double,void*);
double find_root_gsl_wraper(root_func func, void *params, double x_max, double x_min)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    double r;
    double r_min = x_min;
    double r_max = x_max;

    gsl_function F;
    // struct param_initialConditions params = {this, phi0, dV0, d2V0, phi_absMin, delta_phi_cutoff};

    F.function = func;
    F.params = params;

    gsl_root_fsolver_set(s, &F, r_min, r_max);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r      = gsl_root_fsolver_root(s);
        r_min  = gsl_root_fsolver_x_lower(s);
        r_max  = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(r_min,r_max,1e-8,1e-10);
    } while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fsolver_free(s);
    return r;
}

VD func_for_rkqc(double r, VD y, void *param)
{
    Tunneling1D *mod = (Tunneling1D*) param;
    return mod->equationOfMotion(r,y);
}
Tunneling1D::Tunneling1D()
{
    V = nullptr;
    dV = nullptr;
    d2V = nullptr;

    _rk_calculator.SetDOF(2);
    _rk_calculator.SetODE(func_for_rkqc);
    _rk_calculator.SetParams(this);

    phi_absMin = NAN;
    phi_metaMin = NAN;
    phi_bar = NAN;

    rscale = NAN;

    Spatial_Dim = NAN;
    alpha = NAN;

    phi_eps_rel = 1e-3;
    phi_eps_abs = NAN;
}
void Tunneling1D::SetMinima(double absMin, double metaMin)
{
    phi_absMin = absMin;
    phi_metaMin = metaMin;
}
void Tunneling1D::SetPotential(ScalarFunction potential)
{
    V = potential;
}
void Tunneling1D::SetdPotential(dScalarFunction dpotential)
{
    dV = dpotential;
}
void Tunneling1D::SetHM(HM d2potential)
{
    d2V = d2potential;
}
void Tunneling1D::SetPhiAtBarrier(double phibar)
{
    phi_bar = phibar;
}
void Tunneling1D::SetSpatialDim(double dim)
{
    Spatial_Dim = dim;
    alpha = Spatial_Dim - 1;
}
void Tunneling1D::SetPrecision(double eps_rel)
{
    phi_eps_rel = eps_rel;
    // if (isnan(phi_absMin) || isnan(phi_metaMin))
    // {
    //     phi_eps_abs = NAN;
    //     cout<<"Please Set the Location of the Minimum First"<<endl;
    // }
    // else
    // {
    //     phi_eps_abs = phi_eps_rel * abs(phi_absMin - phi_metaMin);
    // }
}
double Tunneling1D::dV_from_absMin(double delta_phi)
{
    double phi = phi_absMin + delta_phi;
    double dV_f = dV({phi},0)[0];

    double dV_d = d2V({phi},0)[0][0] * delta_phi;

    double blend_factor = exp(-pow(delta_phi/phi_eps_abs,2));

    return dV_f*(1-blend_factor) + dV_d*blend_factor;
}
double Tunneling1D::findBarrierLocation()
{
    double phi_tol = abs(phi_absMin-phi_metaMin)*1e-12;
    double V_meta = V({phi_metaMin},0);
    double phiH = phi_metaMin;
    double phiL = phi_absMin;
    double phiM = (phiH + phiL)/2;

    double V0;
    while (abs(phiH-phiL) > phi_tol)
    {
        V0 = V({phiM},0);
        if (V0 > V_meta)
        {
            phiH = phiM;
        }
        else
        {
            phiL = phiM;
        }
        phiM = (phiH + phiL)/2;
    }
    return phiM;
}
double func_for_findRScale(double x, void *params)
{
    Tunneling1D * mod = (Tunneling1D*)params;
    return -mod->VvalatX(x);
}
double Tunneling1D::findRScale()
{
    if (isnan(phi_bar)) findBarrierLocation();
    double phi_tol = abs(phi_bar - phi_metaMin)*1e-6;
    double x1 = min(phi_bar,phi_metaMin);
    double x2 = max(phi_bar,phi_metaMin);

    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);

    double low = x1, high = x2;
    double pred = (low + high)/2;

    gsl_function F;
    F.function = &func_for_findRScale;
    F.params = this;

    gsl_min_fminimizer_set(s, &F, pred, low, high);

    do
    {
        iter++;
        status = gsl_min_fminimizer_iterate(s);
        pred = gsl_min_fminimizer_x_minimum(s);
        low = gsl_min_fminimizer_x_lower(s);
        high = gsl_min_fminimizer_x_upper(s);

        status = gsl_min_test_interval(low,high,phi_tol,0);
    } while (status == GSL_CONTINUE && iter < max_iter);
    
    double phi_bar_top = pred;
    gsl_min_fminimizer_free(s);

    if (phi_bar_top + phi_tol > x2 || phi_bar_top - phi_tol < x1)
    {
        cout<<"In findRScale: No Barrier for the potential: can't find the top position."<<endl;
    }
    

    double Vtop = VvalatX(phi_bar_top) - VvalatX(phi_metaMin);
    double xtop = phi_bar_top - phi_metaMin;

    if (Vtop <= 0)
    {
        cout<<"In findRScale: No Barrier for the potential: non-positive barrier height."<<endl;
    }
    // This is the `time'-scale for a harmonic oscillation (Approximating the potential by a quadratic potential).
    return abs(xtop)/sqrt(abs(2*Vtop));
}
tuple<double, double> Tunneling1D::exactSolution(double r, double phi0, double dV, double d2V)
{
    // ! Find phi(r) (and dphi(r)) given phi(0) assuming a quadractic potential
    double beta = sqrt(abs(d2V));
    double beta_r = beta*r;
    double nu = (alpha - 1)/2;

    double phi = 0;
    double dphi = 0;
    double tmp;
    if (beta_r < 1e-2)
    {
        // Using the expansion to approximate the Bessel function
        double s = d2V>0?1:-1;
        for (int k = 1; k < 4; k++)
        {
            tmp = pow(beta_r/2,2*k-2)*pow(s,k)/(tgamma((k+1)*1.0)*tgamma(k+1+nu));
            phi += tmp;
            dphi += tmp*(2*k);
        }
        phi *= tgamma(nu+1)*r*r*dV*s/4;
        dphi *= tgamma(nu+1)*r*dV*s/4;
        phi += phi0;
    }
    else if (d2V > 0)
    {
        phi = (tgamma(nu+1)/pow(beta_r/2,nu)*cyl_bessel_i(nu,beta_r)-1)*dV/d2V;
        dphi = -nu/r/pow(beta_r/2,nu)*cyl_bessel_i(nu,beta_r);
        dphi += beta/2/pow(beta_r/2,nu)*(cyl_bessel_i(nu-1,beta_r)+cyl_bessel_i(nu+1,beta_r));
        dphi *= tgamma(nu+1)*dV/d2V;
        phi += phi0;
    }
    else
    {
        phi = (tgamma(nu+1)/pow(beta_r/2,nu)*cyl_bessel_j(nu,beta_r)-1)*dV/d2V;
        dphi = -nu/r/pow(beta_r/2,nu)*cyl_bessel_j(nu,beta_r);
        dphi += beta/2/pow(beta_r/2,nu)*(cyl_bessel_j(nu-1,beta_r)-cyl_bessel_j(nu+1,beta_r));
        dphi *= tgamma(nu+1)*dV/d2V;
        phi += phi0;
    }
    return make_tuple(phi,dphi);    
}
struct param_initialConditions
{
    Tunneling1D *tun;
    double phi0;
    double dV0;
    double d2V0;

    double phi_absMin;
    double delta_phi_cutoff;
};
double func_for_initialConditions(double r, void *params)
{
    param_initialConditions *mod = (param_initialConditions*)params;
    double phi0 = mod->phi0;
    double dV0 = mod->dV0;
    double d2V0 = mod->d2V0;
    double phi_absMin = mod->phi_absMin;
    double delta_phi_cutoff = mod->delta_phi_cutoff;
    double phir,dphir;
    tie(phir,dphir) = (mod->tun)->exactSolution(r,phi0,dV0,d2V0);
    return abs(phir - phi_absMin) - abs(delta_phi_cutoff);
}
tuple<double, double, double> Tunneling1D::initialConditions(double delta_phi0, double rmin, double delta_phi_cutoff)
{
    /* 
    * Find the initial conditions for the ODE integration.
    * 
    * The instanton equations of motion are singular at `r=0`, 
    * so we need to start the integration at some larger radius. 
    * This function finds the value `r0` such that `phi(r0) = phi_cutoff`.
    * If there is no such value, it returns the intial conditions at `rmin`.
    */
   
    double phi0 = phi_absMin + delta_phi0;
    double dV0 = dV_from_absMin(delta_phi0);
    double d2V0 = d2V({phi0},0)[0][0];

    double phi_rmin, dphi_rmin;
    tie(phi_rmin, dphi_rmin) = exactSolution(rmin, phi0, dV0, d2V0);
    if (abs(phi_rmin - phi_absMin) > abs(delta_phi_cutoff))
    {
        return make_tuple(rmin, phi_rmin, dphi_rmin);
    }
    if (sign(dphi_rmin) != sign(delta_phi0))
    {
        return make_tuple(rmin, phi_rmin, dphi_rmin);
    }
    
    double r_cur = rmin;
    double r_last = rmin;

    double phi, dphi;

    while (std::isfinite(r_cur))
    {
        r_last = r_cur;
        r_cur *= 10;
        tie(phi, dphi) = exactSolution(r_cur, phi0, dV0, d2V0);
        if (abs(phi - phi_absMin) > abs(delta_phi_cutoff))
        {
            break;
        }
    }
    struct param_initialConditions params = {this, phi0, dV0, d2V0, phi_absMin, delta_phi_cutoff};
    double r = find_root_gsl_wraper(&func_for_initialConditions,&params,r_cur,r_last);

    tie(phi,dphi) = exactSolution(r,phi0,dV0,d2V0);
    return make_tuple(r,phi,dphi);    
}
VD Tunneling1D::equationOfMotion(double r, VD y)
{
    VD res(2);
    res[0] = y[1];
    res[1] = dV({y[0]},0)[0]-alpha*y[1]/r;
    return res;
}
struct cubic_param
{
    double y0;
    double dy0;
    double y1;
    double dy1;
    double diff;
};

double cubicInterpolation(double x, void *param)
{
    cubic_param* mod = (cubic_param*)param;
    double mt = 1-x;
    double c3 = mod->y1;
    double c2 = mod->y1 - mod->dy1/3.0;
    double c1 = mod->y0 + mod->dy0/3.0;
    double c0 = mod->y0;
    return c0*pow(mt,3) + 3*c1*mt*mt*x + 3*c2*mt*x*x + c3*pow(x,3) - mod->diff;
}
tuple<double, VD, CONVERGENCETYPE> Tunneling1D::integrateProfile(double r0, VD y0, double dr0, double epsfrac, double epsabs, double drmin, double rmax)
{
    VD y_final_value = {phi_metaMin,0};
    VD y_diff;
    double dr_guess = dr0;
    double dr_did,dr_next;
    double r = r0;
    VD y = y0;
    VD dydr = equationOfMotion(r,y);
    double r_cache;
    VD y_cache;
    VD dydr_cache;
    VD y_scale;
    VD y_inter(2);
    int ysign = sign(y0[0]-phi_metaMin);
    rmax += r0;

    CONVERGENCETYPE convergQ = NONE;
    cubic_param inter_param;
    double x;
    while (true)
    {
        y_scale = abs(y)+abs(dydr*dr_guess);
        r_cache = r;
        y_cache = y;
        dydr_cache = dydr;
        _rk_calculator._RKQC_SingleStep(r_cache,y_cache,dydr_cache,dr_guess,epsabs,y_scale,dr_did,dr_next);
        dydr_cache = equationOfMotion(r_cache,y_cache);

        y_diff = abs(y_cache-y_final_value);
        if ( y_diff[0] < epsabs && y_diff[1] < epsabs)
        {
            convergQ = CONVERGED;
            break;
        }
        
        if (y_cache[1]*ysign > 0)
        {
            // This means the `ball` is heading back, so it will never reach the desired point.
            convergQ = UNDERSHOOT;
            inter_param = {y[1],dydr[1]*dr_did,y_cache[1],dydr_cache[1]*dr_did,0};
            x = find_root_gsl_wraper(&cubicInterpolation,&inter_param,1,0);
            r += dr_did*x;
            y_inter[1] = cubicInterpolation(x,&inter_param);
            inter_param = {y[0],dydr[0]*dr_did,y_cache[0],dydr_cache[0]*dr_did,0};
            y_inter[0] = cubicInterpolation(x,&inter_param);
            y = y_inter;
            break;
        }

        if ((y_cache[0]-phi_metaMin)*ysign<0)
        {
            // Already passing the desired ending point
            convergQ = OVERSHOOT;
            inter_param = {y[0],dydr[0]*dr_did,y_cache[0],dydr_cache[0]*dr_did,phi_metaMin};
            x = find_root_gsl_wraper(&cubicInterpolation,&inter_param,1,0);
            r += dr_did*x;
            inter_param = {y[1],dydr[1]*dr_did,y_cache[1],dydr_cache[1]*dr_did,0};
            y_inter[1] = cubicInterpolation(x,&inter_param);
            inter_param = {y[0],dydr[0]*dr_did,y_cache[0],dydr_cache[0]*dr_did,0};
            y_inter[0] = cubicInterpolation(x,&inter_param);
            y = y_inter;
            break;
        }

        r = r_cache;
        y = y_cache;
        dydr = dydr_cache;
        dr_guess = dr_next;
    }
    y_diff = abs(y-y_final_value);
    if ( y_diff[0] < epsabs && y_diff[1] < epsabs)
    {
        convergQ = CONVERGED;
    }
    return make_tuple(r,y,convergQ);
}

tuple<VD, VD, VD, double> integrateAndSaveProfile(VD R, VD y0, double dr, double epsfrac, double epsabs, double drmin)
{
    int N = R.size();
    double r0 = R[0];
    VVD Yout(N,VD(y0.size(),0));
    Yout[0] = y0;
    
}