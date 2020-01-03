/*
 * @Description  : overload operators
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-30 18:55:43
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-03 10:06:23
 */
#include <VTypes.h>
using namespace std;


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