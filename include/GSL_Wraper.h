#ifndef GSL_WRAPER_H
#define GSL_WRAPER_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
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
    VVD _Y_transposed;
    VD _X;
    int _NFields;

    int _K;
    int _NCOEFFS;
    int _NBREAKS;
    int _NDataPoints;  

    gsl_bspline_workspace *_bw = nullptr;
    gsl_multifit_linear_workspace *_mw = nullptr;

    gsl_vector* _B = nullptr; // Store the coefficients at each x;
    gsl_matrix* _dB = nullptr; // Store the coefficients at each x
    gsl_matrix* _XC = nullptr; // Store the coefficients at all x, which then will be used for fitting
    std::vector<gsl_vector*> _Cs; // The coefficient from fitting
    // gsl_vector *_Weight; // The weight in fitting
    std::vector<gsl_matrix*> _COVs; // The covariant matrix

    void Fitting();

public:
    GSL_BSpline_Fit(int k, int ncoeffs);
    GSL_BSpline_Fit(VVD Y, VD X, int k, int ncoeffs);
    ~GSL_BSpline_Fit();

    void SetDataX(VD X);
    void UpdateDataY(VVD Y);
    VD valAt(double x);
    VVD valAt(VD X);
    std::tuple<VD,VD> derivAt(double x); // first and second derivative;
    std::tuple<VVD,VVD> derivAt(VD X);
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
    VD valAt(VD X);

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