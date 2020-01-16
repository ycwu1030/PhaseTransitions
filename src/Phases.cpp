/*
 * @Description  : Phase class
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-21 21:08:54
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-16 16:11:10
 */
#include "Phases.h"
#include <iostream>
#include <sstream>
#include <algorithm>


using namespace std;

bool VVDSort(VD l1, VD l2)
{
    return l1[0]<l2[0];
}
int Phase::NPhases = 0;

Phase::Phase()
{
    acc = gsl_interp_accel_alloc();
}

Phase::Phase(int _key, VVD _X, VD _T, VVD _dXdT)
{
    acc = gsl_interp_accel_alloc();
    SetPhase(_key,_X,_T,_dXdT);
    SetUpInterpolation();
    NPhases++;
}

Phase::Phase(const Phase &ph)
{
    acc = gsl_interp_accel_alloc();
    key = ph.key;
    DimX = ph.DimX;
    DimData = ph.DimData;
    X = ph.X;
    T = ph.T;
    dXdT = ph.dXdT;
    low_trans = ph.low_trans;
    high_trans = ph.high_trans;
    SetUpInterpolation();
}

Phase::Phase(const Phase &ph1, const Phase &ph2)
{
    // ! When using this constructor, please make sure one of the two phases are redundant.
    acc = gsl_interp_accel_alloc();
    key = NPhases;
    NPhases++;
    double tmin1 = ph1.GetTmin();
    double tmin2 = ph2.GetTmin();
    double tmax1 = ph1.GetTmax();
    double tmax2 = ph2.GetTmax();
    double tmin = max(tmin1,tmin2);
    double tmax = min(tmax1,tmax2);
    if (tmin >= tmax && tmin1 < tmin2)
    {
        DimX = ph1.DimX;
        DimData = ph1.DimData+ph2.DimData;
        X = VVD(ph1.X);
        dXdT = VVD(ph1.dXdT);
        for (int i = 0; i < DimX; i++)
        {
            X[i].insert(X[i].end(),ph2.X[i].begin()+1,ph2.X[i].end());
            dXdT[i].insert(dXdT[i].end(),ph2.dXdT[i].begin()+1,ph2.dXdT[i].end());
        }
        T = VD(ph1.T);
        T.insert(T.end(),ph2.T.begin(),ph2.T.end());
        low_trans = ph1.low_trans;
        low_trans.insert(ph2.low_trans.begin(),ph2.low_trans.end());
        high_trans = ph1.high_trans;
        high_trans.insert(ph2.high_trans.begin(),ph2.high_trans.end());
    }
    else if (tmin >= tmax && tmin2 <= tmin1)
    {
        DimX = ph2.DimX;
        DimData = ph1.DimData+ph2.DimData;
        X = VVD(ph2.X);
        dXdT = VVD(ph2.dXdT);
        for (int i = 0; i < DimX; i++)
        {
            X[i].insert(X[i].end(),ph1.X[i].begin()+1,ph1.X[i].end());
            dXdT[i].insert(dXdT[i].end(),ph1.dXdT[i].begin()+1,ph1.dXdT[i].end());
        }
        T = VD(ph2.T);
        T.insert(T.end(),ph1.T.begin(),ph1.T.end());
        low_trans = ph2.low_trans;
        low_trans.insert(ph1.low_trans.begin(),ph1.low_trans.end());
        high_trans = ph2.high_trans;
        high_trans.insert(ph1.high_trans.begin(),ph1.high_trans.end());
    }
    else if (tmin1>=tmin2 && tmax1<=tmax2)
    {
        DimX = ph2.DimX;
        DimData = ph2.DimData;
        X = ph2.X;
        T = ph2.T;
        dXdT = ph2.dXdT;
        low_trans = ph2.low_trans;
        high_trans = ph2.high_trans;
    }
    else if (tmin2>tmin1 && tmax2<tmax1)
    {
        DimX = ph1.DimX;
        DimData = ph1.DimData;
        X = ph1.X;
        T = ph1.T;
        dXdT = ph1.dXdT;
        low_trans = ph1.low_trans;
        high_trans = ph1.high_trans;
    } 
    else if (tmin1>=tmin2 && tmax1>=tmax2)
    {
        DimX = ph2.DimX;
        X = ph2.X;
        T = ph2.T;
        dXdT = ph2.dXdT;
        for (int i = 0; i < ph1.DimData; i++)
        {
            if (ph1.T[i]>tmax)
            {
                T.push_back(ph1.T[i]);
                for (int j = 0; j < DimX; j++)
                {
                    X[j].push_back(ph1.X[j][i]);
                    dXdT[j].push_back(ph1.dXdT[j][i]);
                }
            }
        }
        DimData = T.size();
        low_trans = ph2.low_trans;
        low_trans.insert(ph1.low_trans.begin(),ph1.low_trans.end());
        high_trans = ph2.high_trans;
        high_trans.insert(ph1.high_trans.begin(),ph1.high_trans.end());
    }
    else if (tmin1<tmin2 && tmax1<tmax2)
    {
        DimX = ph1.DimX;
        X = ph1.X;
        T = ph1.T;
        dXdT = ph1.dXdT;
        for (int i = 0; i < ph2.DimData; i++)
        {
            if (ph2.T[i]>tmax)
            {
                T.push_back(ph2.T[i]);
                for (int j = 0; j < DimX; j++)
                {
                    X[j].push_back(ph2.X[j][i]);
                    dXdT[j].push_back(ph2.dXdT[j][i]);
                }
            }
        }
        DimData = T.size();
        low_trans = ph1.low_trans;
        low_trans.insert(ph2.low_trans.begin(),ph2.low_trans.end());
        high_trans = ph1.high_trans;
        high_trans.insert(ph2.high_trans.begin(),ph2.high_trans.end());
    }
    else
    {
        DimX = ph1.DimX;
        DimData = ph1.DimData;
        X = ph1.X;
        T = ph1.T;
        dXdT = ph1.dXdT;
        low_trans = ph1.low_trans;
        high_trans = ph1.high_trans;
    }
    SetUpInterpolation();
}

void Phase::SetPhase(int _key, VVD _X, VD _T, VVD _dXdT)
{
    key = _key;
    DimX = _X[0].size();
    DimData = _T.size();
    VVD _TXdXdT;
    for (size_t i = 0; i < DimData; i++)
    {
        VD _current;
        _current.push_back(_T[i]);
        for (size_t j = 0; j < _X[i].size(); j++)
        {
            _current.push_back(_X[i][j]);
        }
        for (size_t j = 0; j < _dXdT[i].size(); j++)
        {
            _current.push_back(_dXdT[i][j]);
        }
        _TXdXdT.push_back(_current);
    }
    sort(_TXdXdT.begin(),_TXdXdT.end(),VVDSort);
    for (size_t i = 0; i < DimX; i++)
    {
        X.push_back(vector<double>(DimData,0.0));
        dXdT.push_back(vector<double>(DimData,0.0));
    }
    for (size_t i = 0; i < _TXdXdT.size(); i++)
    {
        T.push_back(_TXdXdT[i][0]);
        for (size_t j = 0; j < DimX; j++)
        {
            X[j][i] = _TXdXdT[i][j+1];
        }
        for (size_t j = 0; j < DimX; j++)
        {
            dXdT[j][i] = _TXdXdT[i][j+1+DimX];   
        }
    }
}

void Phase::SetUpInterpolation()
{
    FreeInterpolation();
    for (size_t i = 0; i < DimX; i++)
    {
        if (DimData >= 3)
        {
            inters.push_back(gsl_spline_alloc(gsl_interp_steffen,DimData));
            gsl_spline_init(inters[i],T.data(),X[i].data(),DimData);
        }
        else if (DimData >= 2)
        {
            inters.push_back(gsl_spline_alloc(gsl_interp_linear,DimData));
            gsl_spline_init(inters[i],T.data(),X[i].data(),DimData);
        }
        // ! if DimData == 1, no need to interpolate, just return its value;  
    }
}

void Phase::FreeInterpolation()
{
    for (size_t i = 0; i < inters.size(); i++)
    {
        gsl_spline_free(inters[i]);
    }
    inters.clear();
}

Phase::~Phase()
{
    FreeInterpolation();
    gsl_interp_accel_free(acc); 
}

VD Phase::valAt(double _T)
{
    VD res(DimX,0);
    if (_T<T[0])
    {
        cout<<"Error[Phase::valAt]: Required Temperature "<<_T<<" is less than the minimum temperature "<<T[0]<<". Evaluate at "<<T[0]<<" instead."<<endl;
        _T = T[0];
    }
    if (_T>T[DimData-1])
    {
        cout<<"Error[Phase::valAt]: Required Temperature "<<_T<<" is larger than the maximum temperature "<<T[DimData-1]<<". Evaluate at "<<T[DimData-1]<<" instead."<<endl;
        _T = T[DimData-1];
    }
    for (size_t i = 0; i < DimX; i++)
    {
        if (DimData >= 2)
        {
            res[i]=gsl_spline_eval(inters[i],_T,acc);
        }
        else
        {
            res[i]=X[i][DimData-1];
        }
    }
    return res;
}

void Phase::addLinkFrom(Phase *other_phase)
{
    if (T[0] >= other_phase->GetTmax())
    {
        addLowTrans(other_phase->GetKey());
        other_phase->addHighTrans(key);
    }
    if (T[DimData-1] <= other_phase->GetTmin())
    {
        addHighTrans(other_phase->GetKey());
        other_phase->addLowTrans(key);
    }   
}

string Phase::repr()
{
    char sT[200];
    string sX, sdX;
    stringstream ssX, ssdX;
    stringstream tmp;
    string res;
    sprintf(sT,"[%.3f,...,%.3f]",T[0],T[DimData-1]);
    ssX<<"[["<<X[0][0];
    ssdX<<"[["<<dXdT[0][0];
    for (size_t i = 1; i < DimX; i++)
    {
        ssX<<","<<X[i][0];
        ssdX<<","<<dXdT[i][0];
    }
    ssX<<"],...,["<<X[0][DimData-1];
    ssdX<<"],...,["<<dXdT[0][DimData-1];
    for (size_t i = 1; i < DimX; i++)
    {
        ssX<<","<<X[i][DimData-1];
        ssdX<<","<<dXdT[i][DimData-1];
    }
    ssX<<"]]";
    ssdX<<"]]";
    ssX>>sX;
    ssdX>>sdX;
    tmp<<"Phase(key="<<key<<",X="<<sX<<",T="<<sT<<",dXdT="<<sdX<<")";
    tmp>>res;
    // printf("%s",res.c_str());
    return res;
}