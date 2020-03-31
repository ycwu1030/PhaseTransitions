#include <cmath>
#include "GSL_Wraper.h"
#include <iostream>
#include <fstream>

using namespace std;

/* Number of data points to fit */
#define N 100

/* Number of fit coefficients */
#define NCOEFFS 10

#define K 4 // Cubic order
/* NBREAK = NCOEFFS + 2 - K */
#define NBREAK (NCOEFFS + 2 - K)

VVD pathDeriv(VVD pts)
{
    VVD res;
    int n = pts.size();
    if (n >= 5)
    {
        res = deriv14_const_dx(pts);
    }
    else if (n > 2)
    {
        res = VVD(n);
        res[0] = -1.5*pts[0] + 2.0*pts[1] - 0.5*pts[2];
        for (int i = 1; i < n-1; i++)
        {
            res[i] = (pts[i+1]-pts[i-1])/2.0;
        }
        res[n-1] = 1.5*pts[n-1] - 2.0*pts[n-2] + 0.5*pts[n-3];
    }
    else
    {
        res = VVD(n);
        for (size_t i = 0; i < n; i++)
        {
            res[n] = pts[n-1]-pts[0];
        }
    }
    return res;
}

int main(int argc, char const *argv[])
{
    const size_t n = N;
    const size_t ncoeffs = NCOEFFS;
    const size_t nbreak = NBREAK;
    
    double phiX0 = 200;
    double phiY0 = 140;

    VVD phi;
    for (size_t i = 0; i < 100; i++)
    {
        double angle = i/99.0*M_PI_2;
        phi.push_back({phiX0*cos(angle),phiY0*sin(angle)});
    }
    VVD dphi = pathDeriv(phi);
    VD p_dist = cumtrapz(pow(dphi*dphi,0.5));
    double L = p_dist.back();
    VD t = p_dist/L;
    VVD phi_lin;
    VD DelPhi = phi.back()-phi.front();
    for (size_t i = 0; i < t.size(); i++)
    {
        phi_lin.push_back(phi.front()+t[i]*DelPhi);
    }
    // cout<<t<<endl;
    GSL_BSpline_Fit fitter0(phi,t,K,ncoeffs);
    GSL_BSpline_Fit fitter1(phi-phi_lin,t,K,ncoeffs);
    GSL_BSpline_Fit fitter2(phi,p_dist,K,ncoeffs);
    
    VVD dphi0,dphi1,dphi2;
    VVD d2phi0,d2phi1,d2phi2;
    for (size_t i = 0; i < t.size(); i++)
    {
        VD _dphi,_d2phi;
        tie(_dphi,_d2phi) = fitter0.derivAt(t[i]);
        dphi0.push_back(_dphi);
        d2phi0.push_back(_d2phi);

        tie(_dphi,_d2phi) = fitter1.derivAt(t[i]);
        dphi1.push_back(_dphi + DelPhi);
        d2phi1.push_back(_d2phi);

        tie(_dphi,_d2phi) = fitter2.derivAt(p_dist[i]);
        dphi2.push_back(_dphi);
        d2phi2.push_back(_d2phi);
    }

    VD dphi0_dist = pow(dphi0*dphi0,0.5);
    VD dphi1_dist = pow(dphi1*dphi1,0.5);
    VD dphi2_dist = pow(dphi2*dphi2,0.5);

    ofstream output("BSplineFitter.dat");
    output<<"t\tpdist\tphix\tphiy\tdphix0\tdphiy0\td2phix0\td2phiy0\tdist0\tdphix1\tdphiy1\td2phix1\td2phiy1\tdist1\tdphix2\tdphiy2\td2phix2\td2phiy2\tdist2"<<endl;
    for (size_t i = 0; i < t.size(); i++)
    {
        output<<t[i]<<"\t"<<p_dist[i]<<"\t"<<phi[i]<<"\t"<<dphi0[i]<<"\t"<<d2phi0[i]<<"\t"<<dphi0_dist[i]<<"\t"<<dphi1[i]<<"\t"<<d2phi1[i]<<"\t"<<dphi1_dist[i]<<"\t"<<dphi2[i]<<"\t"<<d2phi2[i]<<"\t"<<dphi2_dist[i]<<endl;
    }
    output.close();
    

    return 0;
}

