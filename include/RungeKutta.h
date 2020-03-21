/*
 * @Description  : Following Numerical Recipes in C (1988 version), Section 15.1 and 15.2
*/
#ifndef RungeKutta_H
#define RungeKutta_H

#include "VTypes.h"
#include <string>

typedef VD (*F_ODE)(double x, VD y, void *param);

class RungeKutta
{
private:
    int _DOF;
    double _x_begin, _x_end;
    VD _X; // Position, dimension depends on whether we use the adaptive stepsize control.
    VVD _Y; // The vector (as the same dimension as _X) of vector (dimension _DOF) of Y
    VVD _dYdX; // The vector (as the same dimension as _X) of vector (dimension _DOF) of dY/dX
    VD _BOUND_CONDITION;  // The boundary condition at the starting point
    VD _Y_SCALE; // The scale 
    F_ODE _derivs; // The ODE functions, arguments are the position and the Y values 
    void *_param; // The parameters that will be passed to _derivs

    void _RESET(); // Whenever we change the DOF, boundary condition or the ODE function, we need to reset X,Y,dYdX and be prepared for re-do the RK iteration;
    void _RK4_SingleStep(double X_CUR, VD Y_CUR, VD dY_CUR, double step_size, VD &dY_NEXT); // This is the usual 4th-order Runge-Kutta method, take one step forward.

public:
    RungeKutta();
    RungeKutta(int DOF);
    ~RungeKutta(){};

    void SetDOF(int DOF=1);
    void SetBound(double x_begin, double x_end, VD BOUND); // This is used for usual ODE problem where the boundary conditions are given only at x_begin;
    void SetODE(F_ODE derivs);
    void SetParams(void *param);

    void _RKQC_SingleStep(double &X, VD &Y, VD dY, double step_size_guess, double eps, VD Y_Scale, double &step_size_did, double &step_size_further); // This is the Runge-Kutta method with quality controlled, which will achieve 5-th order accuracy. (adaptive stepsize)
    void ODEINTEGRAL(double step_start, double eps=1e-6);

    void PrintSolution();
    void DumpSolution(std::string filename);

};

#endif