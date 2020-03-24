/*
 * @Description  : Unit test for ApproxLocalMin
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-31 12:49:52
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-03-24 11:47:34
 */

#include "TraceMin.h"
#include "Phases.h"
#include <iostream>
#include <cmath>

using namespace std;

double Potential(VD X, void *T)
{
    double h = X[0];
    double s = X[1];
    return sin(h);
}
int main(int argc, char const *argv[])
{
    double Pi = 3.141592653589793;
    VD x0(2);
    VD x1(2);
    x0[0] = Pi*3.0/2.0;
    x0[1] = 0.0;
    x1[0] = x0[0] + Pi*4;
    VVD res = findApproxLocalMin(Potential,x0,x1,0,100);
    for (int i = 0; i < res.size(); i++)
    {
        cout<<"x_mid_"<<i<<": "<<endl;
        cout<<"\t(";
        for (int j = 0; j < res[i].size(); j++)
        {
            cout<<" "<<res[i][j]<<" ";
        }
        cout<<")"<<endl;
    }
    return 0;
}
