#include "Tunneling1D.h"
#include <cmath>
#include <iostream>
// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_min.h>
// #include <gsl/gsl_roots.h>
#include "GSL_Wraper.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace boost::math;

// typedef double (*root_func)(double,void*);
// double find_root_gsl_wraper(root_func func, void *params, double x_max, double x_min)
// {
//     int status;
//     int iter = 0, max_iter = 100;
//     const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
//     gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
//     double r;
//     double r_min = x_min;
//     double r_max = x_max;

//     gsl_function F;
//     // struct param_initialConditions params = {this, phi0, dV0, d2V0, phi_absMin, delta_phi_cutoff};

//     F.function = func;
//     F.params = params;

//     gsl_root_fsolver_set(s, &F, r_min, r_max);

//     do
//     {
//         iter++;
//         status = gsl_root_fsolver_iterate(s);
//         r      = gsl_root_fsolver_root(s);
//         r_min  = gsl_root_fsolver_x_lower(s);
//         r_max  = gsl_root_fsolver_x_upper(s);
//         status = gsl_root_test_interval(r_min,r_max,1e-8,1e-10);
//     } while (status == GSL_CONTINUE && iter < max_iter);
    
//     gsl_root_fsolver_free(s);
//     return r;
// }

VD func_for_rkqc(double r, VD y, void *param)
{
    Tunneling1D *mod = (Tunneling1D*) param;
    return mod->equationOfMotion(r,y);
}
Tunneling1D::Tunneling1D(double absMin, double metaMin, ScalarFunction V_, dScalarFunction dV_, HM d2V_, double dim, double phi_eps_rel_)
{
    V = V_;
    dV = dV_;
    d2V = d2V_;

    _rk_calculator.SetDOF(2);
    _rk_calculator.SetODE(func_for_rkqc);
    _rk_calculator.SetParams(this);

    phi_absMin = absMin;
    phi_metaMin = metaMin;
    phi_bar = findBarrierLocation();
    // cout<<"Barrier at "<<phi_bar<<endl;

    rscale = findRScale();
    // cout<<"r-scale: "<<rscale<<endl;

    Spatial_Dim = dim;
    alpha = dim-1;

    phi_eps_rel = phi_eps_rel_;
    phi_eps_abs = phi_eps_rel*abs(phi_absMin-phi_metaMin);
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
    double T = 0;
    double dV_f = dV({phi},&T)[0];

    double dV_d = d2V({phi},&T)[0][0] * delta_phi;

    double blend_factor = exp(-pow(delta_phi/phi_eps_abs,2));

    return dV_f*(1-blend_factor) + dV_d*blend_factor;
}
double Tunneling1D::findBarrierLocation()
{
    double phi_tol = abs(phi_absMin-phi_metaMin)*1e-12;
    double T = 0;
    double V_meta = V({phi_metaMin},&T);
    double phiH = phi_metaMin;
    double phiL = phi_absMin;
    double phiM = (phiH + phiL)/2;

    double V0;
    while (abs(phiH-phiL) > phi_tol)
    {
        V0 = V({phiM},&T);
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

    double phi_bar_top = find_min_arg_gsl_wraper(func_for_findRScale,this,x2,x1,phi_tol);

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
tuple<double, double> Tunneling1D::exactSolution(double r, double phi0, double dV_, double d2V_)
{
    // ! Find phi(r) (and dphi(r)) given phi(0) assuming a quadractic potential
    double beta = sqrt(abs(d2V_));
    double beta_r = beta*r;
    double nu = (alpha - 1.0)/2.0;

    double phi = 0;
    double dphi = 0;
    double tmp;
    if (beta_r < 1e-2)
    {
        // Using the expansion to approximate the Bessel function
        double s = d2V_>0?1:-1;
        for (int k = 1; k < 4; k++)
        {
            tmp = pow(beta_r/2,2*k-2)*pow(s,k)/(tgamma((k+1)*1.0)*tgamma(k+1+nu));
            phi += tmp;
            dphi += tmp*(2*k);
        }
        phi *= tgamma(nu+1)*r*r*dV_*s/4;
        dphi *= tgamma(nu+1)*r*dV_*s/4;
        phi += phi0;
    }
    else if (d2V_ > 0)
    {
        // cout<<"beta_r: "<<beta_r<<endl;
        try
        {
            phi = (tgamma(nu+1)/pow(beta_r/2,nu)*cyl_bessel_i(nu,beta_r)-1)*dV_/d2V_;
            dphi = -nu/r/pow(beta_r/2,nu)*cyl_bessel_i(nu,beta_r);
            dphi += beta/2/pow(beta_r/2,nu)*(cyl_bessel_i(nu-1,beta_r)+cyl_bessel_i(nu+1,beta_r));
            dphi *= tgamma(nu+1)*dV_/d2V_;
            phi += phi0;
        }
        catch(const boost::wrapexcept<std::overflow_error>& e)
        {
            // std::cerr << e.what() << '\n';
            // just ignore the overflow
            phi = INFINITY;
            dphi = INFINITY;
        }
    }
    else
    {
        phi = (tgamma(nu+1)/pow(beta_r/2,nu)*cyl_bessel_j(nu,beta_r)-1)*dV_/d2V_;
        dphi = -nu/r/pow(beta_r/2,nu)*cyl_bessel_j(nu,beta_r);
        dphi += beta/2/pow(beta_r/2,nu)*(cyl_bessel_j(nu-1,beta_r)-cyl_bessel_j(nu+1,beta_r));
        dphi *= tgamma(nu+1)*dV_/d2V_;
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
    std::tie(phir,dphir) = (mod->tun)->exactSolution(r,phi0,dV0,d2V0);
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
   
    double T = 0;
    double phi0 = phi_absMin + delta_phi0;
    double dV0 = dV_from_absMin(delta_phi0);
    double d2V0 = d2V({phi0},&T)[0][0];

    double phi_rmin, dphi_rmin;
    std::tie(phi_rmin, dphi_rmin) = exactSolution(rmin, phi0, dV0, d2V0);
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
    r_last = r_cur;
    r_cur *= 10;
    while (std::isfinite(r_cur))
    {
        std::tie(phi, dphi) = exactSolution(r_cur, phi0, dV0, d2V0);
        if (!std::isfinite(phi))
        {
            r_cur = (r_last + r_cur)/2.0;
            continue;
        }
        if (abs(phi - phi_absMin) > abs(delta_phi_cutoff))
        {
            break;
        }
        r_last = r_cur;
        r_cur *= 10;
    }
    struct param_initialConditions params = {this, phi0, dV0, d2V0, phi_absMin, delta_phi_cutoff};
    // cout<<"Before root finding"<<endl;
    // cout<<"r_cur="<<r_cur<<" r_last="<<r_last<<endl;
    double r = find_root_gsl_wraper(&func_for_initialConditions,&params,r_cur,r_last);
    // cout<<"After root finding"<<endl;

    std::tie(phi,dphi) = exactSolution(r,phi0,dV0,d2V0);
    return make_tuple(r,phi,dphi);    
}
VD Tunneling1D::equationOfMotion(double r, VD y)
{
    VD res(2);
    double T = 0;
    res[0] = y[1];
    res[1] = dV({y[0]},&T)[0]-alpha*y[1]/r;
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
        // cout<<"\t\t----"<<endl;
        // cout<<"\t\t"<<r_cache<<"  "<<y_cache[0]<<"  "<<y_cache[1]<<"  "<<dydr_cache[0]<<"  "<<dydr_cache[1]<<"  "<<dr_guess<<"  "<<epsabs<<endl;
        _rk_calculator._RKQC_SingleStep(r_cache,y_cache,dydr_cache,dr_guess,epsabs,y_scale,dr_did,dr_next);
        dydr_cache = equationOfMotion(r_cache,y_cache);

        y_diff = abs(y_cache-y_final_value);
        // cout<<"\t\t"<<r_cache<<"  "<<y_cache[0]<<"  "<<y_cache[1]<<"  "<<ysign<<"  "<<dr_did<<"  "<<dr_next<<endl;
        // cout<<"\t\t\t"<<y_diff[0]<<"/"<<epsabs<<"  "<<y_diff[1]<<"/"<<epsabs<<endl;
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

tuple<VD, VD, VD, double> Tunneling1D::integrateAndSaveProfile(VD R, VD y0, double dr, double epsfrac, double epsabs, double drmin)
{
    int N = R.size();
    double r0 = R[0];
    VVD Yout(y0.size(),VD(N,0));
    Yout[0][0] = y0[0];
    Yout[1][0] = y0[1];
    VD dydr0 = equationOfMotion(r0,y0);
    double Rerr = NAN;

    int i = 1;
    double r = r0;
    VD y = y0;
    VD dydr = dydr0;
    double r_cache;
    VD y_cache;
    VD dydr_cache;
    double dr_guess = dr;
    double dr_did,dr_next;
    VD y_scale(2);
    cubic_param inter_param;
    while (i<N)
    {
        y_scale = abs(y)+abs(dydr*dr_guess);
        r_cache = r;
        y_cache = y;
        dydr_cache = dydr;
        _rk_calculator._RKQC_SingleStep(r_cache,y_cache,dydr_cache,dr_guess,epsabs,y_scale,dr_did,dr_next);
        if (dr_did < drmin)
        {
            y_cache = y + (y_cache-y)*drmin/dr_did;
            dr_did = drmin;
            dr_next = drmin;
            r_cache = r + dr_did;
            if (!(isnan(Rerr)))
            {
                Rerr = r_cache;
            }
        }
        dydr_cache = equationOfMotion(r_cache,y_cache);
        if (r < R[i] && R[i] <= r_cache)
        {
            while (i < N && r < R[i] && R[i] <= r_cache)
            {
                double x = (R[i]-r)/dr_did;
                inter_param = {y[0], dr_did*dydr[0], y_cache[0], dr_did*dydr_cache[0], 0};
                Yout[0][i] = cubicInterpolation(x, &inter_param);
                inter_param = {y[1], dr_did*dydr[1], y_cache[1], dr_did*dydr_cache[1], 0};
                Yout[1][i] = cubicInterpolation(x, &inter_param);
                i += 1;
            }   
        }

        r = r_cache;
        y = y_cache;
        dydr = dydr_cache;
        dr_guess = dr_next;
    }
    
    return make_tuple(R,Yout[0],Yout[1],Rerr);
}
tuple<VD,VD,VD,double> Tunneling1D::findProfile(double xguess,double xtol,double phitol,double thinCutoff,int npoints,double rmin, double rmax, int max_interior_pts)
{
    double xmin = xtol*10;
    double xmax = INFINITY;
    double x;
    if (!isnan(xguess))
    {
        x = xguess;
    }
    else
    {
        x = - log(abs((phi_bar-phi_absMin)/(phi_metaMin-phi_absMin)));
    }
    // cout<<"Starting point: x = "<<x<<endl;
    double xincrease = 5.0;

    rmin *= rscale;
    double dr0 = rmin;
    double drmin = rmin*1e-2;
    rmax *= rscale;

    double delta_phi = phi_metaMin - phi_absMin;
    double epsabs = abs(delta_phi*phitol);
    double epsfrac = phitol;
    double delta_phi_cutoff = thinCutoff*delta_phi;

    double rf = NAN;
    double delta_phi0;
    double r0_,phi0,dphi0;
    double r0;
    VD y0;
    VD yf;
    CONVERGENCETYPE ctype;
    // cout<<"Starting of the shooting: "<<endl;
    while (true)
    {
        // cout<<"--------"<<endl;
        delta_phi0 = exp(-x)*delta_phi;
        // cout<<"\tdelta_phi0="<<delta_phi0<<endl;
        std::tie(r0_,phi0,dphi0) = initialConditions(delta_phi0,rmin,delta_phi_cutoff);
        // cout<<"\tInitial condition: r0="<<r0_<<"  phi0="<<phi0<<" dphi0="<<dphi0<<endl;
        if ( !std::isfinite(r0_) || !std::isfinite(x))
        {
            if (isnan(rf))
            {
                cerr<<"Failed to retrieve initial conditions on the first try"<<endl;
            }
            break;
        }
        r0 = r0_;
        y0 = {phi0,dphi0};
        std::tie(rf,yf,ctype) = integrateProfile(r0,y0,dr0,epsfrac,epsabs,drmin,rmax);
        if (ctype == CONVERGED)
        {
            break;
        }
        else if (ctype == UNDERSHOOT)
        {
            xmin = x;
            x = std::isfinite(xmax)?(xmin+xmax)/2:x*xincrease;
        }
        else if (ctype == OVERSHOOT)
        {
            xmax = x;
            x = (xmin+xmax)/2;
        }
        
        if (xmax-xmin < xtol)
        {
            break;
        }   
    }
    
    VD R(npoints);
    for (size_t i = 0; i < npoints; i++)
    {
        R[i] = r0 + i*(rf-r0)/(npoints-1);
    }
    VD Phi_ex;
    VD dPhi_ex;
    double Rerr;
    std::tie(R,Phi_ex,dPhi_ex,Rerr)=integrateAndSaveProfile(R,y0,dr0,epsfrac,epsabs,drmin);

    VD R_int;
    if (max_interior_pts < 0)
    {
        max_interior_pts = R.size()/2;
    }
    if (max_interior_pts > 0)
    {
        double dx0 = R[1]-R[0];
        if (R[0]/dx0 <= max_interior_pts)
        {
            int n = ceil(R[0]/dx0);
            for (size_t i = 0; i < n; i++)
            {
                R_int.push_back(0+i*(R[0])/(n));
            }
        }
        else
        {
            int n = max_interior_pts;
            double a = (R[0]/dx0 - n)*2/(n*(n+1));
            for (size_t i = 0; i < n; i++)
            {
                int k = n-i;
                R_int.push_back(R[0]-dx0*(k + a*k*(k+1)/2));
            }
            R_int[0] = 0.0;
        }  
    }
    VD Phi_int(R_int.size(),0);
    VD dPhi_int(R_int.size(),0);
    Phi_int[0] = phi_absMin + delta_phi0;
    dPhi_int[0] = 0.0;
    double dV_ = dV_from_absMin(delta_phi0);
    double d2V_ = d2V({Phi_int[0]},0)[0][0];
    for (size_t i = 1; i < R_int.size(); i++)
    {
        std::tie(Phi_int[i],dPhi_int[i]) = exactSolution(R_int[i],Phi_int[0],dV_,d2V_);
    }
    VD R_final(R_int);
    VD Phi_final(Phi_int);
    VD dPhi_final(dPhi_int);

    R_final.insert(R_final.end(),R.begin(),R.end());
    Phi_final.insert(Phi_final.end(),Phi_ex.begin(),Phi_ex.end());
    dPhi_final.insert(dPhi_final.end(),dPhi_ex.begin(),dPhi_ex.end());
    
    return make_tuple(R_final,Phi_final,dPhi_final,Rerr);
    
}
double Tunneling1D::findAction(VD R, VD Phi, VD dPhi)
{
    int N = R.size();
    double Sphere_area = pow(M_PI,Spatial_Dim/2.0)/tgamma(Spatial_Dim/2.0);
    VD area = pow(R,alpha)*Sphere_area;
    VD integrand(N);
    double T = 0;
    for (size_t i = 0; i < N; i++)
    {
        integrand[i] = (pow(dPhi[i],2)/2 + V({Phi[i]},&T) - V({phi_metaMin},&T))*area[i];
    }
    double S = Simpson(R,integrand);

    // For the bulk inside the bubble interior
    double volume = pow(R[0],Spatial_Dim)*pow(M_PI,Spatial_Dim/2.0)/tgamma(Spatial_Dim/2.0+1.0);
    S += volume*(V({Phi[0]},&T)-V({phi_metaMin},&T));
    return S;
}
std::tuple<VD, VD> Tunneling1D::evenlySpacedPhi(VD phi, VD dphi, int npoint, int k, bool fixAbs)
{
    if (fixAbs)
    {
        phi.insert(phi.begin(),phi_absMin);
        phi.insert(phi.end(),phi_metaMin);
        dphi.insert(dphi.begin(),0.0);
        dphi.insert(dphi.end(),0.0);
    }
    else
    {
        phi.insert(phi.end(),phi_metaMin);
        dphi.insert(dphi.end(),0.0);
    }
    
    // Sort phi in increasing order
    VVD fullPhi = transpose({phi,dphi});
    // cout<<"fullPhi dim: ("<<fullPhi.size()<<","<<fullPhi[0].size()<<")"<<endl;
    sort(fullPhi.begin(),fullPhi.end(),[](VD x1, VD x2){return x1[0]<x2[0];});
    VVD::iterator iter = unique(fullPhi.begin(),fullPhi.end(),[](VD x1, VD x2){return x1[0]==x2[0];});
    fullPhi.resize(distance(fullPhi.begin(),iter));
    fullPhi=transpose(fullPhi);
    // cout<<fullPhi[0]<<endl;

    GSL_Spline_Inter inter;
    inter.SetData(&fullPhi[1],&fullPhi[0]);

    VD p;
    if (fixAbs)
    {
        p = linspace(phi_absMin,phi_metaMin,npoint);
    }
    else
    {
        p = linspace(phi[0],phi_metaMin,npoint);
    }

    return make_tuple(p,inter.valAt(p));
}