#ifndef TUNNELING_1D_H
#define TUNNELING_1D_H

#include "VTypes.h"
#include <tuple>


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

    double dV_from_absMin(double delta_phi);
    double findBarrierLocation();
    double findRScale();
    std::tuple<double, double> exactSolution(double r, double phi0, double dV, double d2V);

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
};


#endif //TUNNELING_1D_H