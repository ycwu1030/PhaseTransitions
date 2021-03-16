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

    // double* Z; //Dimension 2NDim
    // double** Q; // Dimension (2NDim) x (2NDim)
    inline double H00(double t);
    inline double H10(double t);
    inline double H01(double t);
    inline double H11(double t);

    int Find_Index(double x); // Find the index i, such that  X[i]<= x < X[i+1], using binary search.
public:
    HermiteInterpolation1D(VD _X, VD _F, VD _dF);
    ~HermiteInterpolation1D();
    double valAt(double x);
};

double HermiteInterpolation1D::H00(double t)
{
    return (1+2*t)*(1-t)*(1-t);
}

double HermiteInterpolation1D::H10(double t)
{
    return t*(1-t)*(1-t);
}

double HermiteInterpolation1D::H01(double t)
{
    return t*t*(3-2*t);
}

double HermiteInterpolation1D::H11(double t)
{
    return t*t*(t-1);
}
#endif