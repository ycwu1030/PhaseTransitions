#ifndef VTypes_H
#define VTypes_H
#include <vector>
#include <map>
#include <tuple>
#include <set>
#include <iostream>

class Phase;
typedef std::vector<double> VD;
typedef std::vector<std::vector<double> > VVD;
typedef std::set<int> SI;
typedef std::map<int, Phase> MP;
typedef std::vector<std::tuple<VD, double> > VT;
typedef std::vector<std::tuple<double, double, VD, int> > VnT;
typedef double (*ScalarFunction)(VD,double);
typedef VD (*dScalarFunction)(VD,double);
typedef VVD (*HM)(VD,double);
typedef bool (*forbidCrit)(VD);

inline bool AllPass(VD x){return false;}
std::vector<double> abs(const vector<double> &input);

std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs);
std::vector<double> operator+(const std::vector<double> &lhs, const double &cons);
std::vector<double> operator+(const double &cons, const std::vector<double> &rhs);
std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs);


std::vector<double> operator-(const std::vector<double> &lhs, const std::vector<double> &rhs);


std::vector<double> operator*(const std::vector<double> &lhs, const double &s);
std::vector<double> operator*(const double &s, const std::vector<double> &rhs);
double operator*(const std::vector<double> &lhs, const std::vector<double> &rhs); // Scalar Product

VD operator/(const VD &lhs, const VD &rhs); // elementary-wise divide
VD operator/(const VD &lhs, const double &s);

std::ostream& operator<<(std::ostream& out, const VD& s);

VVD operator/(const VVD &lhs, const double &s);

#endif