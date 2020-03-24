#include "Tunneling1D.h"
#include <iostream>
#include <fstream>

using namespace std;

double Vthin(VD Phi, double T)
{
    double phi = Phi[0];
    return 0.25*pow(phi,4) - 0.49*pow(phi,3) + 0.235*pow(phi,2);
}
VD dVthin(VD Phi, double T)
{
    double phi = Phi[0];
    return {pow(phi,3) - 3*0.49*pow(phi,2) + 0.47*phi};
}
VVD d2Vthin(VD Phi, double T)
{
    double phi = Phi[0];
    VVD res(1,VD(1,0));
    res[0][0] = 3*pow(phi,2) - 2*3*0.49*phi + 0.47;
    return res;
}
double Vthick(VD Phi, double T)
{
    double phi = Phi[0];
    return 0.25*pow(phi,4) - 0.4*pow(phi,3) + 0.1*pow(phi,2);
}
VD dVthick(VD Phi, double T)
{
    double phi = Phi[0];
    return {pow(phi,3) - 1.2*pow(phi,2) + 0.2*phi};
}
VVD d2Vthick(VD Phi, double T)
{
    double phi = Phi[0];
    VVD res(1,VD(1,0));
    res[0][0] = 3*pow(phi,2) - 2.4*phi + 0.2;
    return res;
}

int main(int argc, char const *argv[])
{
    Tunneling1D solver_thin(1.0,0.0,Vthin,dVthin,d2Vthin);
    VD R;
    VD Phi;
    VD dPhi;
    double Rerr;
    tie(R,Phi,dPhi,Rerr) = solver_thin.findProfile();
    ofstream output("Thin_Profile.dat");
    for (size_t i = 0; i < R.size(); i++)
    {
        output << R[i] << "\t" << Phi[i] << "\t" << dPhi[i] << endl;
    }
    output.close();
    Tunneling1D solver_thick(1.0,0.0,Vthick,dVthick,d2Vthick);
    tie(R,Phi,dPhi,Rerr) = solver_thick.findProfile();
    ofstream output1("Thick_Profile.dat");
    for (size_t i = 0; i < R.size(); i++)
    {
        output1 << R[i] << "\t" << Phi[i] << "\t" << dPhi[i] << endl;
    }
    output1.close();
    return 0;
}

