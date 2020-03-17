/*
 * @Description  : Hermite Interpolation
 */
#include "HermiteInterpolation.h"
using namespace std;
// * Using https://www3.nd.edu/~zxu2/acms40390hw/Hermite.cpp Failed
// * Good for NDim is not that large.
// * Has problem for NDim is large.
HermiteInterpolation1D::HermiteInterpolation1D(VD _X, VD _F, VD _dF)
{
    X = _X;
    F = _F;
    dF = _dF;
    NDim = X.size();
    Z = new double[2*NDim];
    Q = new double*[2*NDim];
    for (size_t i = 0; i < NDim; i++)
    {
        Q[2*i] = new double[2*NDim];
        Q[2*i+1] = new double[2*NDim];
        Z[2*i] = X[i];
        Z[2*i+1] = X[i];
        Q[2*i][0] = F[i];
        Q[2*i+1][0] = F[i];
        Q[2*i+1][1] = dF[i];
        if (i!=0)
        {
            Q[2*i][1] = (Q[2*i][0]-Q[2*i-1][0])/(Z[2*i]-Z[2*i-1]);
        } 
    }
    int K = 2*(NDim-1) + 1;
    for (size_t i = 2; i <= K; i++)
    {
        for (size_t j = 2; j <= i; j++)
        {
            Q[i][j] = (Q[i][j-1]-Q[i-1][j-1])/(Z[i]-Z[i-j]);
        }
    }
}

HermiteInterpolation1D::~HermiteInterpolation1D()
{
    delete[] Z;
    for (size_t i = 0; i < 2*NDim; i++)
    {
        delete[] Q[i];
    }
    delete[] Q;
}

double HermiteInterpolation1D::valAt(double x)
{
    int K = 2*NDim-1;
    int J;
    double S = Q[K][K]*(x-Z[K-1]);
    for (size_t i = 2; i <= K; i++)
    {
        J = K - i + 1;
        S = (S+Q[J][J])*(x-Z[J-1]);
    }
    S += Q[0][0];
    return S;
}