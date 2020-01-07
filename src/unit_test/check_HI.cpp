/*
 * @Description  : Unit test for Hermite Interpolation
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-22 19:30:17
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-03 10:07:16
 */

#include "HermiteInterpolation.h"
#include <cmath>

double F(double x)
{
    return 1/x;
}
double dF(double x)
{
    return -1/x/x;
}

int main(int argc, char const *argv[])
{
    double Pi = 3.14159265;
    int N = 100;
    VD X;
    VD FF;
    VD DFF;
    double xx;
    for (size_t i = 0; i < N; i++)
    {
        xx = 1+i/((double)N-1)*2*Pi;
        X.push_back(xx);
        FF.push_back(F(xx));
        DFF.push_back(dF(xx));
    }
    HermiteInterpolation1D inter(X,FF,DFF);
    for (size_t i = 0; i < N; i++)
    {
        printf("%d At %f: %f, %f, %f\n",i,X[i],FF[i], DFF[i], inter.valAt(X[i])-FF[i]);
    }
    return 0;
}
