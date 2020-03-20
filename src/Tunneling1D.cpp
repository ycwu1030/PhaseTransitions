#include "Tunneling1D.h"
#include <cmath>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace boost::math;

Tunneling1D::Tunneling1D()
{
    V = nullptr;
    dV = nullptr;
    d2V = nullptr;

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