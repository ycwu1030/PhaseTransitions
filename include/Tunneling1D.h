#ifndef TUNNELING_1D_H
#define TUNNELING_1D_H

#include "VTypes.h"
#include "RungeKutta.h"
#include <tuple>

enum CONVERGENCETYPE
{
    UNDERSHOOT = -1,
    CONVERGED = 0,
    OVERSHOOT = 1,
    NONE      = -9
};

class Tunneling1D
{
private:
    double phi_absMin;
    double phi_metaMin;
    double phi_bar;
    ScalarFunction V;
    dScalarFunction dV;
    HM d2V;

    double phi_eps_rel;
    double phi_eps_abs;

    double Spatial_Dim;
    double alpha;

    double rscale;

    RungeKutta _rk_calculator;

    double dV_from_absMin(double delta_phi);
    double findBarrierLocation();
    double findRScale();
    std::tuple<double, double, double> initialConditions(double delta_phi0, double rmin, double delta_phi_cutoff);
    std::tuple<double, VD, CONVERGENCETYPE> integrateProfile(double r0, VD y0, double dr0, double epsfrac, double epsabs, double drmin, double rmax);
    std::tuple<VD, VD, VD, double> integrateAndSaveProfile(VD R, VD y0, double dr, double epsfrac, double epsabs, double drmin);

public:
    Tunneling1D();
    ~Tunneling1D(){};

    void SetMinima(double absMin, double metaMin);
    void SetPotential(ScalarFunction potential);
    void SetdPotential(dScalarFunction dpotential);
    void SetHM(HM d2potential);

    void SetPhiAtBarrier(double phibar);
    void SetSpatialDim(double dim);
    void SetPrecision(double eps_rel);

    double VvalatX(double x) {return V({x},0);}

    std::tuple<double, double> exactSolution(double r, double phi0, double dV, double d2V);

    VD equationOfMotion(double r, VD y);
};


#endif //TUNNELING_1D_H