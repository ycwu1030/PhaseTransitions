/*
 * @Description  : Unit test for tracing single minimum
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-30 20:37:46
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-03-24 11:51:30
 */
#include "TraceMin.h"
#include "Phases.h"
#include <iostream>
#include <cmath>

using namespace std;

double Potential(VD X, void *param)
{
    double T = *((double *)param);
    double mu2h = -pow(120,2);
    double mu2s = -pow(230,2);
    double lamh = 1.4;
    double lams = 2.1;
    double lamhs = 1.5;
    double ch = 0.66;
    double cs = 0.72;
    double h = X[0];
    double s = X[1];
    double res = 0;
    res += (mu2h + ch*T*T)*h*h;
    res += (mu2s + cs*T*T)*s*s;
    res += lamh*pow(h,4);
    res += lams*pow(s,4);
    res += lamhs*pow(h*s,2);
    return res;
}
VD dPotential_dX(VD X, void *param)
{
    double T = *((double *)param);
    double mu2h = -pow(120,2);
    double mu2s = -pow(230,2);
    double lamh = 1.4;
    double lams = 2.1;
    double lamhs = 1.5;
    double ch = 0.66;
    double cs = 0.72;
    double h = X[0];
    double s = X[1];
    VD res(2);
    res[0] = (mu2h + ch*T*T)*h*2 + 4*lamh*pow(h,3) + 2*lamhs*s*s*h;
    res[1] = (mu2s + cs*T*T)*s*2 + 4*lams*pow(s,3) + 2*lamhs*s*h*h;
    return res;
}
VD d2Potential_dXdT(VD X, void *param)
{
    double T = *((double *)param);
    double mu2h = -pow(120,2);
    double mu2s = -pow(230,2);
    double lamh = 1.4;
    double lams = 2.1;
    double lamhs = 1.5;
    double ch = 0.66;
    double cs = 0.72;
    double h = X[0];
    double s = X[1];
    VD res(2);
    res[0] = 4*ch*T*h;
    res[1] = 4*cs*T*s;
    return res;
}
VVD d2Potential_dX2(VD X, void *param)
{
    double T = *((double *)param);
    double mu2h = -pow(120,2);
    double mu2s = -pow(230,2);
    double lamh = 1.4;
    double lams = 2.1;
    double lamhs = 1.5;
    double ch = 0.66;
    double cs = 0.72;
    double h = X[0];
    double s = X[1];
    VVD res;
    VD resh(2);
    resh[0] = (mu2h + ch*T*T)*2 + 12*lamh*pow(h,2) + 2*lamhs*s*s;
    resh[1] = 4*lamhs*s*h;
    res.push_back(resh);
    VD ress(2);
    ress[0] = 4*lamhs*s*h;
    ress[1] = (mu2s + cs*T*T)*2 + 12*lams*pow(s,2) + 2*lamhs*h*h;
    res.push_back(ress);
    return res;
}

int main(int argc, char const *argv[])
{
    VD x0(2);
    x0[0] = 0.0;
    x0[1] = 112.229;
    _traceMinimum_rval res = traceMinimum(Potential,dPotential_dX,d2Potential_dXdT, d2Potential_dX2,x0,0,500,0.5);
    cout<<res.T.size()<<endl;
    Phase ph(0,res.X,res.T,res.dXdT);
    cout<<ph.repr()<<endl;
    return 0;
}
