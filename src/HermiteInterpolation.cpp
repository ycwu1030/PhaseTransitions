/*
 * @Description  : Hermite Interpolation
 */
#include "HermiteInterpolation.h"
#include <iostream>
using namespace std;
// * Using https://www3.nd.edu/~zxu2/acms40390hw/Hermite.cpp Failed
// * Good for NDim is not that large.
// * Has problem for NDim is large.
// * So, I use the method described in Wiki: https://en.wikipedia.org/wiki/Cubic_Hermite_spline
HermiteInterpolation1D::HermiteInterpolation1D(VD _X, VD _F, VD _dF)
{
    X = _X;
    F = _F;
    dF = _dF;
    NDim = X.size();
    // Z = new double[2*NDim];
    // Q = new double*[2*NDim];
    // for (size_t i = 0; i < NDim; i++)
    // {
    //     Q[2*i] = new double[2*NDim];
    //     Q[2*i+1] = new double[2*NDim];
    //     Z[2*i] = X[i];
    //     Z[2*i+1] = X[i];
    //     Q[2*i][0] = F[i];
    //     Q[2*i+1][0] = F[i];
    //     Q[2*i+1][1] = dF[i];
    //     if (i!=0)
    //     {
    //         Q[2*i][1] = (Q[2*i][0]-Q[2*i-1][0])/(Z[2*i]-Z[2*i-1]);
    //     } 
    // }
    // int K = 2*(NDim-1) + 1;
    // for (size_t i = 2; i <= K; i++)
    // {
    //     for (size_t j = 2; j <= i; j++)
    //     {
    //         Q[i][j] = (Q[i][j-1]-Q[i-1][j-1])/(Z[i]-Z[i-j]);
    //     }
    // }
}

int HermiteInterpolation1D::Find_Index(double x)
{
    // Find the index i, such that  X[i]<= x < X[i+1], using binary search.
    // if x < X[0] return -1
    // if x > X[-1] return NDim
    if (x<X[0])
    {
        return -1;
    }
    if (x>X[NDim-1])
    {
        return NDim;
    }
    if (x==X[NDim-1])
    {
        return NDim-2;
    }
    int indmin = 0;
    int indmax = NDim-2;
    int ind = (indmax+indmin)/2;
    while (true)
    {
        // cout<<x<<" ["<<X[ind]<<","<<X[ind+1]<<"]"<<endl;
        if (x<X[ind])
        {
            indmax = ind;
            ind = (indmax+indmin)/2;
            continue;
        }
        if (x>=X[ind+1])
        {
            indmin = ind+1;
            ind = (indmax+indmin)/2;
            continue;
        }
        break;
    }
    return ind;
}
HermiteInterpolation1D::~HermiteInterpolation1D()
{
    // delete[] Z;
    // for (size_t i = 0; i < 2*NDim; i++)
    // {
    //     delete[] Q[i];
    // }
    // delete[] Q;
}

double HermiteInterpolation1D::valAt(double x)
{
    // int K = 2*NDim-1;
    // int J;
    // double S = Q[K][K]*(x-Z[K-1]);
    // for (size_t i = 2; i <= K; i++)
    // {
    //     J = K - i + 1;
    //     S = (S+Q[J][J])*(x-Z[J-1]);
    // }
    // S += Q[0][0];
    // return S;
    int ind = Find_Index(x);
    if (ind < 0 || ind == NDim)
    {
        return 0;
    }
    double xmax = X[ind+1];
    double xmin = X[ind];
    double t = (x-xmin)/(xmax-xmin);
    return H00(t)*F[ind]+H10(t)*(xmax-xmin)*dF[ind]+H01(t)*F[ind+1]+H11(t)*(xmax-xmin)*dF[ind+1];
}