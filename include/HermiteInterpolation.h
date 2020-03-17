/*
 * @Description  : Containing the class doing Hermite Interpolation
 */
#ifndef Hermite_Interpolation_H
#define Hermite_Interpolation_H
#include "VTypes.h"

// ! Not working, need to be improved.
class HermiteInterpolation1D
{
private:
    VD X;// Dimension NDim
    VD F;// Dimension NDim
    VD dF;// Dimension NDim
    int NDim;

    double* Z; //Dimension 2NDim
    double** Q; // Dimension (2NDim) x (2NDim)
public:
    HermiteInterpolation1D(VD _X, VD _F, VD _dF);
    ~HermiteInterpolation1D();
    double valAt(double x);
};


#endif