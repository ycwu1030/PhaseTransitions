#include "RungeKutta.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

#define POW_GROW (-0.2)
#define POW_SHRINK (-0.25)
#define SAFETY (0.9)
#define MAXSTEPS (10000)
#define TINY (1e-30)

bool abs_compare(double i, double j) {return abs(i)<abs(j);}

RungeKutta::RungeKutta()
{
    SetDOF();
    _derivs = nullptr;
    _param = nullptr;
}
RungeKutta::RungeKutta(int DOF)
{
    SetDOF(DOF);
    _derivs = nullptr;
    _param = nullptr;
}
void RungeKutta::SetDOF(int DOF)
{
    _DOF = DOF;
}
void RungeKutta::SetBound(double x_begin, double x_end, VD BOUND)
{
    _x_begin = x_begin;
    _x_end = x_end;
    _BOUND_CONDITION = BOUND;
}
void RungeKutta::SetODE(F_ODE derivs)
{
    _derivs = derivs;
}
void RungeKutta::SetParams(void *param)
{
    _param = param;
}
void RungeKutta::_RESET()
{
    _X.clear();
    _Y.clear();
    _dYdX.clear();
    _X.push_back(_x_begin);
    _Y.push_back(_BOUND_CONDITION);
}
void RungeKutta::_RK4_SingleStep(double X_CUR, VD Y_CUR, VD dY_CUR, double step_size, VD &Y_NEXT)
{
    double half_step = step_size/2;

// The first step, just use the current dY/dX 
    VD dY_Step1 = step_size*dY_CUR;

// The second step, half_step in x, half dY_Step1 in Y
    VD dY_Step2 = step_size*_derivs(X_CUR+half_step,Y_CUR+dY_Step1/2,_param);

// The third step, half_step in x, half dY_Step2 in Y
    VD dY_Step3 = step_size*_derivs(X_CUR+half_step,Y_CUR+dY_Step2/2,_param);

// The fourth step, full step in x, full dY_Step3 in Y
    VD dY_Step4 = step_size*_derivs(X_CUR+step_size,Y_CUR+dY_Step3,_param);
    
    Y_NEXT = Y_CUR + dY_Step1/6 + dY_Step2/3 + dY_Step3/3 + dY_Step4/6;
}

void RungeKutta::_RKQC_SingleStep(double &X, VD &Y, VD dY, double step_size_guess, double eps, VD Y_Scale, double &step_size_did, double &step_size_further)
{
    double step_size = step_size_guess;
    double half_step_size;
    double error_max = 0;
    double min_step_size = 1e-5*step_size_guess; // Since we will adapt the step size according to the estimated error, we need to terminate such operate at some point, otherwise, we will stuck at here for some case.
    double max_step_size; // We will enlarge the step size, when we are within the precision. But we need to limit it within a reasonable range.
    double step_size_temp;

    double x_cache = X;
    VD y_cache = Y;
    VD dy_cache = dY;
    VD y_temp;
    VD Delta_y;
    VD error_temp;

    while (true)
    {
        // Take two half steps
        // cout<<"Step Size: "<<step_size<<endl;
        half_step_size = step_size/2.0;
        _RK4_SingleStep(x_cache,y_cache,dy_cache,half_step_size,y_temp);
        X = x_cache + half_step_size;
        dY = _derivs(X,y_temp,_param);
        _RK4_SingleStep(X,y_temp,dY,half_step_size,Y);

        // Take one large step
        X = x_cache + step_size;
        _RK4_SingleStep(x_cache,y_cache,dy_cache,step_size,y_temp);
        Delta_y = Y - y_temp;
        error_temp = abs(Delta_y/Y_Scale);
        // error_temp = abs(y_temp);
        error_max = *max_element(error_temp.begin(),error_temp.end());
        error_max /= eps;
        if (error_max <= 1.0)
        {
            step_size_did = step_size;
            step_size_temp = SAFETY*step_size*exp(POW_GROW*log(error_max));
            max_step_size = 4*step_size_did;
            step_size_further = min(step_size_temp,max_step_size);
            break;
        }
        step_size_temp = step_size*SAFETY*exp(POW_SHRINK*log(error_max));
        if (abs(step_size_temp) < abs(min_step_size))
        {
            step_size_did = step_size;
            step_size_further = step_size;
            break;
        }
        step_size = step_size_temp;
    }
    Y = Y + Delta_y/15;
}

void RungeKutta::ODEINTEGRAL(double step_start,double eps)
{
    _RESET();
    double x = _X[0];
    VD y = _Y[0];
    VD dydx = _derivs(x,y,_param);
    _dYdX.push_back(dydx);
    double step_size = (_x_end > _x_begin)? abs(step_start):-abs(step_start);
    double step_size_did;
    double step_size_next;
    for (size_t nstp = 0; nstp < MAXSTEPS; ++nstp)
    {
        _Y_SCALE = abs(y)+abs(dydx*step_size);
        for (size_t i = 0; i < _DOF; i++)
        {
            _Y_SCALE[i] = min(1.0,_Y_SCALE[i]);
        }
        // Check whether the step size is too big, and already pass the end point
        if ( (step_size > 0 && x+step_size > _x_end) || (step_size < 0 && x+step_size < _x_end) )
        {
            step_size = _x_end - x;
        }
        _RKQC_SingleStep(x,y,dydx,step_size,eps,_Y_SCALE,step_size_did,step_size_next);
        _X.push_back(x);
        _Y.push_back(y);
        dydx=_derivs(x,y,_param);
        _dYdX.push_back(dydx);
        if ((x - _x_end)*(_x_end-_x_begin)>=0)
        {
            // We are finished;
            return;
        }
        // step_size = abs(step_size_next) > 1e-3*abs(_x_end-_x_begin)?step_size_next:1e-3*(_x_end-_x_begin);
        step_size = step_size_next;
        if (abs(step_size_next) > 1e-1*abs(_x_end-_x_begin))
        {
            step_size = 1e-1*(_x_end-_x_begin);
        }
        if (abs(step_size_next) < 1e-5*abs(_x_end-_x_begin))
        {
            step_size = 1e-5*(_x_end-_x_begin);
        }
    }
    cout<<"TAKE TOO MANY STEPS"<<endl;
}

void RungeKutta::PrintSolution()
{
    cout<<"The Solution is:"<<endl;
    cout<<"x\t";
    for (size_t i = 0; i < _DOF; i++)
    {
        cout<<"y_"<<i<<"\t";
    }
    for (size_t i = 0; i < _DOF; i++)
    {
        cout<<"dy_"<<i<<"/dx"<<"\t";
    }
    cout<<endl;
    for (size_t i = 0; i < _X.size(); i++)
    {
        cout<<_X[i]<<"\t";
        for (size_t j = 0; j < _DOF; j++)
        {
            cout<<_Y[i][j]<<"\t";
        }
        for (size_t j = 0; j < _DOF; j++)
        {
            cout<<_dYdX[i][j]<<"\t";
        }
        cout<<endl;
    }
}
void RungeKutta::DumpSolution(string filename)
{
    ofstream output(filename.c_str());
    output<<"The Solution is:"<<endl;
    output<<"x\t";
    for (size_t i = 0; i < _DOF; i++)
    {
        output<<"y_"<<i<<"\t";
    }
    for (size_t i = 0; i < _DOF; i++)
    {
        output<<"dy_"<<i<<"/dx"<<"\t";
    }
    output<<endl;
    for (size_t i = 0; i < _X.size(); i++)
    {
        output<<_X[i]<<"\t";
        for (size_t j = 0; j < _DOF; j++)
        {
            output<<_Y[i][j]<<"\t";
        }
        for (size_t j = 0; j < _DOF; j++)
        {
            output<<_dYdX[i][j]<<"\t";
        }
        output<<endl;
    }
}