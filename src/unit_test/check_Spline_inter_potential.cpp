#include <cmath>
#include "PathDeformation.h"
#include <iostream>
#include <fstream>

using namespace std;

/* Number of data points to fit */
#define N 150
double Potential(VD X, void *param)
{
    double T = *((double *)param);
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

int main(int argc, char const *argv[])
{
    const int n = N;

    // ! Generating the point linking two local minima
    // ! (0,126.313), (101.186,0) @ T = 10 GeV
    
    
    ofstream fdata("Spline_Potential_Data.dat");
    fdata<<"phi0\tphi1\ts\tV"<<endl;
    double T = 10;
    VD min1 = {0,126.313};
    VD min2 = {101.186,0};
    VVD pts;
    VD dphi = (min2-min1)/(n-1);
    VD dist;
    double phi0;
    double phi1;
    double s_max = 0;//sqrt(dist*dist);
    for (int i = 0; i < n; i++)
    {
        phi0 = sqrt(i)*min2[0]/sqrt(n-1);
        phi1 = (n-1-i)*(n-1-i)*min1[1]/pow(n-1,2);
        pts.push_back({phi0,phi1});
        dist=i==0?pts[i]-min1:pts[i] - pts[i-1];
        s_max += sqrt(dist*dist);
        fdata<<pts[i][0]<<"\t"<<pts[i][1]<<"\t"<<s_max<<"\t"<<Potential(pts[i],&T)<<endl;
    }
    
    fdata.close();

    
    SplinePath spt(pts,Potential,T);
    double V,dV,ddV;
    VD pts_pre;
    ofstream fout("Spline_Potential_Inter.dat");
    fout<<"s\tV\tdV\td2V\tphi0\tphi1"<<endl;
    for (double s = -0.1*s_max; s <= s_max*1.1; s+=1.5)
    {
        pts_pre = spt.pts_at_dist(s);
        V = spt.V({s},nullptr);
        dV = spt.dV({s},nullptr)[0];
        ddV = spt.d2V({s},nullptr)[0][0];
        fout<<s<<"\t"<<V<<"\t"<<dV<<"\t"<<ddV<<"\t"<<pts_pre<<endl;
    }
    fout.close();
    
    return 0;
}

