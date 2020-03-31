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
VD operator-(const VD &lhs, const double &rhs)
{
    VD res;
    for (size_t i = 0; i < lhs.size(); i++)
    {
        res.push_back(lhs[i]-rhs);
    }
    return res;
}
VD operator-(const double &lhs, const VD &rhs)
{
    VD res;
    for (size_t i = 0; i < rhs.size(); i++)
    {
        res.push_back(lhs-rhs[i]);
    }
    return res;
}
VD operator-(const VD &rhs)
{
    return 0-rhs;
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
VVD operator*(const VVD &lhs, const double &s)
{
    VVD res;
    for (size_t i = 0; i < lhs.size(); i++)
    {
        res.push_back(lhs[i]*s);
    }
    return res;
}
VVD operator*(const double &s, const VVD &rhs)
{
    VVD res;
    for (size_t i = 0; i < rhs.size(); i++)
    {
        res.push_back(rhs[i]*s);
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
VVD operator+(const VVD &lhs, const VVD &rhs)
{
    VVD res;
    for (size_t i = 0; i < lhs.size(); i++)
    {
        res.push_back(lhs[i]+rhs[i]);
    }
    return res;
}
VVD operator-(const VVD &lhs, const VVD &rhs)
{
    VVD res;
    for (size_t i = 0; i < lhs.size(); i++)
    {
        res.push_back(lhs[i]-rhs[i]);
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
VD operator*(const VVD &lhs, const VVD &rhs)
{
    VD res;
    for (size_t i = 0; i < lhs.size(); i++)
    {
        res.push_back(lhs[i]*rhs[i]);
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

VVD transpose(const VVD &mat)
{
    int row_size = mat.size();
    if (row_size==0)
    {
        return VVD(0);
    }
    int col_size = mat[0].size();
    VVD res(col_size,VD(row_size,0));
    for (int i = 0; i < row_size; i++)
    {
        for (int j = 0; j < col_size; j++)
        {
            res[j][i] = mat[i][j];
        }
    }
    return res;
}

VD linspace(double start, double end, int n)
{
    VD res;
    double dx;
    double x_end;
    if (n == 1)
    {
        dx = end-start;
        x_end = start;
    }
    else
    {
        dx = (end-start)/(n-1);
        x_end = end;
    }
    for (int i = 0; i < n; i++)
    {
        res.push_back(start+i*dx);
    }
    res.back() = x_end;
    return res;
}
VD linspace(double start, double end, double dx)
{
    VD res;
    double x_cur = start;
    while (x_cur<=end)
    {
        res.push_back(x_cur);
        x_cur+=dx;
    }
    return res;
}
// template<class T>
// T deriv14_const_dx(T y, double dx)
// {
//     int N = y.size();
//     T dy;
//     dy.push_back(-25.0*y[0]+48.0*y[1]-36.0*y[2]+16.0*y[3]-3*y[4]);
//     dy.push_back(-3.0*y[0]-10.0*y[1]+18.0*y[2]-6.0*y[3]+y[4]);
//     for (int i = 2; i < N-2; i++)
//     {
//         dy.push_back(y[i-2]-8.0*y[i-1]+8.0*y[i+1]-y[i+2]);
//     }
//     dy.push_back(3.0*y[N-1]+10.0*y[N-2]-18.0*y[N-3]+6.0*y[N-4]-y[N-5]);
//     dy.push_back(25.0*y[N-1]-48.0*y[N-2]+36.0*y[N-3]-16.0*y[N-4]+3.0*y[N-5]);

//     return dy/(12.0*dx);
// }
VD cumtrapz(VD pts, VD x, double dx, double initial)
{
    int Nsize=pts.size();
    bool default_dx = (x.size()!=Nsize);
    VD res(Nsize);
    res[0] = initial;
    double dx_used;
    for (size_t i = 1; i < Nsize; i++)
    {
        dx_used = default_dx?dx:x[i]-x[i-1];
        res[i] = res[i-1] + (pts[i]+pts[i-1])*dx_used/2.0;
    }
    return res;
}