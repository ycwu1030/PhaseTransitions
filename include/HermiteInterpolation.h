/*
 * @Description  : Containing the class doing Hermite Interpolation
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-22 12:45:50
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-03 09:32:41
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