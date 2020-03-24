#include <VTypes.h>
#include <cmath>
using namespace std;

vector<double> abs(const vector<double> &input)
{
    vector<double> res;
    for (int i = 0; i < input.size(); ++i)
    {
        res.push_back(abs(input[i]));
    }
    return res;
}

vector<double> operator+(const vector<double> &lhs, const vector<double> &rhs){
    vector<double> res;
    for (int i = 0; i < lhs.size(); i++)
    {
        res.push_back(lhs[i]+rhs[i]);
    }
    return res;
}

vector<double> operator+(const vector<double> &lhs, const double &cons){
    vector<double> res;
    for (int i = 0; i < lhs.size(); i++)
    {
        res.push_back(lhs[i]+cons);
    }
    return res;
}

vector<double> operator+(const double &cons, const vector<double> &rhs){
    vector<double> res;
    for (int i = 0; i < rhs.size(); i++)
    {
        res.push_back(cons+rhs[i]);
    }
    return res;
}


vector<double> operator-(const vector<double> &lhs, const vector<double> &rhs){
    vector<double> res;
    for (int i = 0; i < lhs.size(); i++)
    {
        res.push_back(lhs[i]-rhs[i]);
    }
    return res;
}


vector<double> operator*(const vector<double> &lhs, const double &s)
{
    vector<double> res;
    for (int i = 0; i < lhs.size(); i++)
    {
        res.push_back(lhs[i]*s);
    }
    return res;
}


vector<double> operator*(const double &s, const vector<double> &rhs)
{
    vector<double> res;
    for (int i = 0; i < rhs.size(); i++)
    {
        res.push_back(rhs[i]*s);
    }
    return res;
}


double operator*(const vector<double> &lhs, const vector<double> &rhs)
{
    double res = 0;
    for (int i = 0; i < lhs.size(); i++)
    {
        res += lhs[i]*rhs[i];
    }
    return res;
}
vector<double> operator/(const vector<double> &lhs, const double &s)
{
    vector<double> res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]/s);
    }
    return res;
}

vector<double> operator/(const vector<double> &lhs, const vector<double> &rhs)
{
    vector<double> res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]/rhs[i]);
    }
    return res;
}
VVD operator/(const VVD &lhs, const double &s)
{
    VVD res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]/s);
    }
    return res;
}

ostream& operator<<(ostream& out, const VD& s)
{
    for (size_t i = 0; i < s.size()-1; i++)
    {
        out<<s[i]<<"\t";
    }
    out<<s[s.size()-1];
    return out;
}

double Simpson(VD X, VD Y)
{
    if (X.size() != Y.size())
    {
        cerr<<"In Simpson Integral: The size of X and Y didn't match."<<endl;
        return NAN;
    }
    int N = X.size();
    int index = (N%2==0)?N/2-1:(N-1)/2;

    double res = 0;
    double x1,x2,x3;
    double y1,y2,y3;
    for (size_t i = 0; i < index; i++)
    {
        x1 = X[2*i];
        x2 = X[2*i+1];
        x3 = X[2*i+2];
        y1 = Y[2*i];
        y2 = Y[2*i+1];
        y3 = Y[2*i+2];
        res += (x3-x1)*(2*x1-3*x2+x3)/(x1-x2)/6.0*y1;
        res += pow(x3-x1,3)/(x1-x2)/(x2-x3)/6.0*y2;
        res += (x1-x3)*(x1-3*x2+2*x3)/(x2-x3)/6.0*y3;
    }
    if (N%2==0)
    {
        res += (X[N-1]-X[N-2])*(Y[N-1]+Y[N-2])/2.0;
    }
    
    return res;
}

VD pow(const VD &input, double power)
{
    int N = input.size();
    VD res(N);
    for (size_t i = 0; i < N; i++)
    {
        res[i] = pow(input[i],power);
    }
    return res;
}