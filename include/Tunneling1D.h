#ifndef TUNNELING_1D_H
#define TUNNELING_1D_H

#include "VTypes.h"
#include "RungeKutta.h"
#include <cmath>
#include <tuple>

enum CONVERGENCETYPE
{
    UNDERSHOOT = -1,
    CONVERGED = 0,
    OVERSHOOT = 1,
    NONE      = -9
};

/*
 * This class is used to calculate the tunneling profile for only one field case
*/
class Tunneling1D
{
private:
    double phi_absMin; // The absolute minimum
    double phi_metaMin; // The meta minimum
    double phi_bar; // The position of the barrier

    // * The potential functions, include first and second derivatives
    ScalarFunction V;
    dScalarFunction dV;
    HM d2V;

    // * These parameters control the precision of the solution
    double phi_eps_rel;
    double phi_eps_abs;

    // * The spatial dimension, alpha = Spatial_Dim - 1;
    // * dim = 4 for T = 0
    // * dim = 3 for T is non-zero
    double Spatial_Dim;
    double alpha;

    // * This is the scale of r
    double rscale;

    // * This is the RungeKutta Integrator
    RungeKutta _rk_calculator;

    // * dV is quite useful when calculating the solution quite close to r = 0, so we need to get it more accurately. 
    double dV_from_absMin(double delta_phi);

    // * Find the position of the barrier
    double findBarrierLocation();

    // * Find the r-scale
    // * Which is calculated to be the `time`-scale of the harmonic oscillation at the barrier top (approximate by -k phi^2)
    double findRScale();

    // * Since the ODE for tunneling calculation is singluar at r = 0, so we need to start at some finite radius.
    // * return the starting r and the phi, dphi/dr at that point
    std::tuple<double, double, double> initialConditions(double delta_phi0, double rmin, double delta_phi_cutoff);


    // * Integrate the ODE, stop whenever we can determine it is overshoot, undershoot or converged.
    // * return the final radius, the phi, dphi/dr at that point, and the status
    std::tuple<double, VD, CONVERGENCETYPE> integrateProfile(double r0, VD y0, double dr0, double epsfrac, double epsabs, double drmin, double rmax);

    // * The same as above, but also save the profile
    std::tuple<VD, VD, VD, double> integrateAndSaveProfile(VD R, VD y0, double dr, double epsfrac, double epsabs, double drmin);

public:
    Tunneling1D(double absMin, double metaMin, ScalarFunction V_, dScalarFunction dV_, HM d2V_, double dim = 3, double phi_eps_rel_ = 1e-3);
    ~Tunneling1D(){};

    void SetMinima(double absMin, double metaMin);
    void SetPotential(ScalarFunction potential);
    void SetdPotential(dScalarFunction dpotential);
    void SetHM(HM d2potential);

    void SetPhiAtBarrier(double phibar);
    void SetSpatialDim(double dim);
    void SetPrecision(double eps_rel);

    double VvalatX(double x) {return V({x},0);}


    // * Provide the exact solution around starting point
    // * By assuming that the potention around that point is in a quadratic form
    std::tuple<double, double> exactSolution(double r, double phi0, double dV, double d2V);

    // * The EOM/ODE equation used for RungeKutta ODE integration
    // * The original ODE is phi'' + alpha/r phi' = dV/dphi
    // * Assume y[0] = phi, y[1] = dphi
    // * Then:
    // *  y[0]' = y[1]
    // *  y[1]' = dV(y[0]) - alpha/r y[1]
    VD equationOfMotion(double r, VD y);


    // * Calculate the bubble profile by calling integrateProfile repeatedly.
    std::tuple<VD,VD,VD,double> findProfile(double xguess=NAN,double xtol=1e-4,double phitol=1e-4,double thinCutoff=0.01,int npoints=500,double rmin=1e-4, double rmax=1e4, int max_interior_pts = -1);

    // * Calculate the Euclidean action for the instanton
    double findAction(VD R, VD Phi, VD dPhi);

    // * From the findProfile itself, the point is kind of evenly spaced in r, (due to the precision control in runge-kutta, the space in r between points can vary, but it is roughly evenly spaced)
    // * For 1D problem itself, it is enough.
    // * But when linking 1D solution to path deformation, if phi is not evenly spaced, when extend the path to the two local minima, we need too much steps, which will tremendously increase the computation time. So for that purpose, we want to make the list evenly spaced in phi.
    std::tuple<VD, VD> evenlySpacedPhi(VD phi, VD dphi, int npoint = 100, int k = 1, bool fixAbs = true);
};


#endif //TUNNELING_1D_H