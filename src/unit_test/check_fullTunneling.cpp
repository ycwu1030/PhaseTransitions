#include "TransitionFinder.h"
#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;

class Potential
{
private:
    double _c;
    double _fx;
    double _fy;

    void SetUpPotential();

public:
    Potential(){};
    Potential(double c, double fx, double fy){_c=c; _fx=fx; _fy=fy; SetUpPotential();}
    ~Potential(){};

    ScalarFunction V;
    dScalarFunction dV;

};

void Potential::SetUpPotential()
{
    V = [&](VD X, void* param){
        double x = X[0];
        double y = X[1];
        double r1 = x*x + _c*y*y;
        double r2 = _c*(x-1)*(x-1) + (y-1)*(y-1);
        double r3 = _fx*(pow(x,4)/4.0 - pow(x,3)/3.0);
        r3 += _fy*(pow(y,4)/4.0 - pow(y,3)/3.0);
        return r1*r2+r3;
    };
    dV = [&](VD X, void* param){
        double x = X[0];
        double y = X[1];
        double r1 = x*x + _c*y*y;
        double r2 = _c*(x-1)*(x-1) + (y-1)*(y-1);
        double dr1dx = 2*x;
        double dr1dy = 2*_c*y;
        double dr2dx = 2*_c*(x-1);
        double dr2dy = 2*(y-1);
        double dVdx = r1*dr2dx + dr1dx*r2 + _fx*x*x*(x-1);
        double dVdy = r1*dr2dy + dr1dy*r2 + _fy*y*y*(y-1);
        VD res = {dVdx, dVdy};
        return res;
    };
}


int main(int argc, char const *argv[])
{
    auto start = chrono::high_resolution_clock::now();
    for (size_t j = 0; j < 10; j++)
    {
    
    VVD pts_init(2);
    pts_init[0] = {1,1};
    pts_init[1] = {0,0};
    Potential Thin(5.0,0.0,2.0);
    TransNucleation thin_nucl = fullTunneling(pts_init,Thin.V,Thin.dV,0);
    ofstream fthin("fullTunneling_thin.dat");
    fthin<<"r\tphix\tphiy\tphir"<<endl;
    for (size_t i = 0; i < thin_nucl.R.size(); i++)
    {
        fthin<<thin_nucl.R[i]<<"\t"<<thin_nucl.Phi[i]<<"\t"<<thin_nucl.Phi_1D[i]<<endl;
    }
    fthin.close();


    Potential Thick(5.0,0.0,80.0);
    TransNucleation thick_nucl = fullTunneling(pts_init,Thick.V,Thick.dV,0);
    ofstream fthick("fullTunneling_thick.dat");
    fthick<<"r\tphix\tphiy\tphir"<<endl;
    for (size_t i = 0; i < thick_nucl.R.size(); i++)
    {
        fthick<<thick_nucl.R[i]<<"\t"<<thick_nucl.Phi[i]<<"\t"<<thick_nucl.Phi_1D[i]<<endl;
    }
    fthick.close();
    
    }
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout<<" Elapsed time: "<<elapsed.count() << "s"<<endl;
    return 0;
}

