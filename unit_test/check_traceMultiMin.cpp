/*
 * @Description  : Unit test for tracing multi-minima
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-01 13:14:03
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-03 10:07:56
 */

#include "TraceMin.h"
#include "Phases.h"
#include <iostream>
#include <cmath>

using namespace std;

double Potential(VD X, double T)
{
    double mu2h = -pow(120,2);
    double mu2s = -pow(160,2);
    double lamh = 1.4;
    double lams = 1.6;
    double lamhs = 2.5;
    double ch = 0.66;
    double cs = 0.72;
    double h = X[0];
    double s = X[1];
    double res = 0;
    res += (mu2h + ch*T*T)/2.0*h*h;
    res += (mu2s + cs*T*T)/2.0*s*s;
    res += lamh*pow(h,4)/4.0;
    res += lams*pow(s,4)/4.0;
    res += lamhs*pow(h*s,2)/2.0;
    return res;
}
VD dPotential_dX(VD X, double T)
{
    double mu2h = -pow(120,2);
    double mu2s = -pow(160,2);
    double lamh = 1.4;
    double lams = 1.6;
    double lamhs = 2.5;
    double ch = 0.66;
    double cs = 0.72;
    double h = X[0];
    double s = X[1];
    VD res(2);
    res[0] = (mu2h + ch*T*T)*h + lamh*pow(h,3) + lamhs*s*s*h;
    res[1] = (mu2s + cs*T*T)*s + lams*pow(s,3) + lamhs*s*h*h;
    return res;
}
VD d2Potential_dXdT(VD X, double T)
{
    double mu2h = -pow(120,2);
    double mu2s = -pow(160,2);
    double lamh = 1.4;
    double lams = 1.6;
    double lamhs = 2.5;
    double ch = 0.66;
    double cs = 0.72;
    double h = X[0];
    double s = X[1];
    VD res(2);
    res[0] = 2*ch*T*h;
    res[1] = 2*cs*T*s;
    return res;
}
VVD d2Potential_dX2(VD X, double T)
{
    double mu2h = -pow(120,2);
    double mu2s = -pow(160,2);
    double lamh = 1.4;
    double lams = 1.6;
    double lamhs = 2.5;
    double ch = 0.66;
    double cs = 0.72;
    double h = X[0];
    double s = X[1];
    VVD res;
    VD resh(2);
    resh[0] = (mu2h + ch*T*T) + 3*lamh*pow(h,2) + lamhs*s*s;
    resh[1] = 2*lamhs*s*h;
    res.push_back(resh);
    VD ress(2);
    ress[0] = 2*lamhs*s*h;
    ress[1] = (mu2s + cs*T*T) + 3*lams*pow(s,2) + lamhs*h*h;
    res.push_back(ress);
    return res;
}

int main(int argc, char const *argv[])
{
    VD x0(2);
    x0[0] = 101.419;
    x0[1] = 0.0;
    VD x1(2);
    x1[0] = 0.0;
    x1[1] = 126.491;
    VT points;
    points.push_back(make_tuple(x0,0.0));
    points.push_back(make_tuple(x1,0.0));
    MP AllPhases = traceMultiMin(Potential,dPotential_dX,d2Potential_dXdT,d2Potential_dX2,points,0,1000,0.1);
    MP::iterator iter;
    for ( iter = AllPhases.begin(); iter != AllPhases.end(); iter++)
    {
        int key = iter->first;
        cout<<key<<" "<<(AllPhases.at(key)).repr()<<endl;
    }
    removeRedundantPhases(Potential,dPotential_dX,AllPhases);
    for ( iter = AllPhases.begin(); iter != AllPhases.end(); iter++)
    {
        int key = iter->first;
        cout<<key<<" "<<(AllPhases.at(key)).repr()<<endl;
    }
    return 0;
}
