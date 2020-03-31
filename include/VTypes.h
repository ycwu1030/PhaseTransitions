#ifndef VTypes_H
#define VTypes_H
#include <vector>
#include <map>
#include <tuple>
#include <set>
#include <functional>
#include <iostream>

class Phase;
typedef std::vector<double> VD;
typedef std::vector<std::vector<double> > VVD;
typedef std::set<int> SI;
typedef std::map<int, Phase> MP;
typedef std::vector<std::tuple<VD, double> > VT;
typedef std::vector<std::tuple<double, double, VD, int> > VnT;
typedef std::function<double(VD,void*)> ScalarFunction;
typedef std::function<VD(VD,void*)> dScalarFunction;
typedef std::function<VVD(VD,void*)> HM;
typedef std::function<bool(VD)> forbidCrit;
// typedef double (*ScalarFunction)(VD,double);
// typedef VD (*dScalarFunction)(VD,double);
// typedef VVD (*HM)(VD,double);
// typedef bool (*forbidCrit)(VD);

inline bool AllPass(VD x){return false;}
std::vector<double> abs(const std::vector<double> &input);
std::vector<double> pow(const std::vector<double> &input, double power);

std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs);
std::vector<double> operator+(const std::vector<double> &lhs, const double &cons);
std::vector<double> operator+(const double &cons, const std::vector<double> &rhs);
std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs);


std::vector<double> operator-(const std::vector<double> &lhs, const std::vector<double> &rhs);
std::vector<double> operator-(const std::vector<double> &lhs, const double &rhs);
std::vector<double> operator-(const double &lhs, const std::vector<double> &rhs);
std::vector<double> operator-(const std::vector<double> &rhs);

std::vector<double> operator*(const std::vector<double> &lhs, const double &s);
std::vector<double> operator*(const double &s, const std::vector<double> &rhs);
double operator*(const std::vector<double> &lhs, const std::vector<double> &rhs); // Scalar Product

VD operator/(const VD &lhs, const VD &rhs); // elementary-wise divide
VD operator/(const VD &lhs, const double &s);

std::ostream& operator<<(std::ostream& out, const VD& s);

VVD operator*(const VVD &lhs, const double &s);
VVD operator*(const double &s, const VVD &rhs);
VVD operator/(const VVD &lhs, const double &s);
VVD operator+(const VVD &lhs, const VVD &rhs);
VVD operator-(const VVD &lhs, const VVD &rhs);
VD operator*(const VVD &lhs, const VVD &rhs);

double Simpson(VD X, VD Y);

VD cumtrapz(VD pts, VD x = {}, double dx = 1.0, double initial=0.0);

VVD transpose(const VVD &mat);

VD linspace(double start, double end, int n);
VD linspace(double start, double end, double dx);

template<class T>
T deriv14_const_dx(T y, double dx = 1.0)
{
    int N = y.size();
    T dy;
    dy.push_back(-25.0*y[0]+48.0*y[1]-36.0*y[2]+16.0*y[3]-3*y[4]);
    dy.push_back(-3.0*y[0]-10.0*y[1]+18.0*y[2]-6.0*y[3]+y[4]);
    for (int i = 2; i < N-2; i++)
    {
        dy.push_back(y[i-2]-8.0*y[i-1]+8.0*y[i+1]-y[i+2]);
    }
    dy.push_back(3.0*y[N-1]+10.0*y[N-2]-18.0*y[N-3]+6.0*y[N-4]-y[N-5]);
    dy.push_back(25.0*y[N-1]-48.0*y[N-2]+36.0*y[N-3]-16.0*y[N-4]+3.0*y[N-5]);

    return dy/(12.0*dx);
}

#endif