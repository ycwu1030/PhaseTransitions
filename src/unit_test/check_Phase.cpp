/*
 * @Description  : Unit test for the Phase Class
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-22 22:23:51
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-16 16:14:57
 */
#include "Phases.h"
#include <iostream>
#include <cmath>

using namespace std;

double F0(double x)
{
    return 1/x;
}
double F1(double x)
{
    return sin(x);
}
double dF0(double x)
{
    return -1/x/x;
}
double dF1(double x)
{
    return cos(x);
}

int main(int argc, char const *argv[])
{
    if (argc < 2)
    {
        cout<<"Usage: "<<argv[0]<<" number_of_points"<<endl;
        return -1;
    }
    
    double Pi = 3.14159265;
    int N = atoi(argv[1]);
    VD X;
    VVD FF;
    VVD DFF;
    double xx;
    for (size_t i = 0; i < N; i++)
    {
        if (N == 1) {xx = 1.0;}
        else {xx = 1+i/((double)N-1)*2*Pi;}
        X.push_back(xx);
        VD atF,atDF;
        atF.push_back(F0(xx));
        atF.push_back(F1(xx));
        FF.push_back(atF);
        atDF.push_back(dF0(xx));
        atDF.push_back(dF1(xx));
        DFF.push_back(atDF);
    }
    Phase ph(0,FF,X,DFF);
    for (size_t i = 0; i < N; i++)
    {
        printf("%d At %f: [",i,X[i]);
        for (size_t j = 0; j < 2; j++)
        {
            printf(" %f ",FF[i][j]);
        }
        printf("], [");
        VD inter = ph.valAt(X[i]);
        for (size_t j = 0; j < 2; j++)
        {
            printf(" %f ",inter[j]);
        }
        printf("]\n");
    }
    cout<<ph.repr()<<endl;
    Phase ph1(ph);
    cout<<ph1.repr()<<endl;
    return 0;
}
