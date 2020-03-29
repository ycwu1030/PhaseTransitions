#ifndef GSL_WRAPER_H
#define GSL_WRAPER_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>
#include "VTypes.h"

typedef double (*root_func)(double,void*);
typedef double (*multimin_func)(const gsl_vector*, void*);
double find_root_gsl_wraper(root_func func, void *params, double x_max, double x_min, double epsabs=1e-8, double epsfrac=1e-10);
double find_min_arg_gsl_wraper(root_func func, void *params, double x_max, double x_min, double epsabs=1e-8, double epsfrac=0);
VD find_multimin_arg_gsl_wraper(multimin_func func, void *params, VD start, double epsabs=1e-2);

class GSL_BSpline_Fit
{
private:
    VVD _Y;
    VD _X;

    int _K;
    int _NCOEFFS;
    int _NBREAKS;  
public:
    GSL_BSpline_Fit(VVD Y, VD X, int k, int ncoeffs);
    ~GSL_BSpline_Fit();
};

class GSL_Spline_Inter
{
private:
    int _size;
    VD *_Y; // I just keep the point to the data, but not own the data
    VD *_X;

    // Since I don't own the data, so I can't guarantee the original data is still there. So I keep some of them first.
    double _xmin;
    double _xmax;
    double _ymin;
    double _ymax;
    double _yaver;

    gsl_interp_accel *_acc = nullptr;
    gsl_spline *_spline = nullptr;

public:
    GSL_Spline_Inter();
    ~GSL_Spline_Inter();

    void SetData(VD *Y, VD *X);
    double valAt(double x,int deri = 0);

};

class GSL_Multi_Spline_Inter
{
private:
    int _NDim;
    int _size;
    // I just keep the point to the data, but not own the data
    // The first index of _Y is the same as _X, the second index of _Y is the field index, so it is not convenience for interpolation. So I will transpose it to _Y_transposed
    VVD *_Y; 
    VD *_X;

    VVD _Y_transposed; // I keep the transposed data;

    // Since I don't own the data, so I can't guarantee the original data is still there. So I keep some of them first.
    double _xmin;
    double _xmax;
    VD _ymin;
    VD _ymax;
    VD _yaver;

    gsl_interp_accel *_acc = nullptr;
    std::vector<gsl_spline*> _splines;

public:
    GSL_Multi_Spline_Inter();
    ~GSL_Multi_Spline_Inter();

    void SetData(VVD *Y, VD *X);
    VD valAt(double x,int deri = 0);

};

#endif //GSL_WRAPER_H