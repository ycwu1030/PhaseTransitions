#include "GSL_Wraper.h"
#include <iostream>
#include <algorithm>

using namespace std;

double find_root_gsl_wraper(root_func func, void *params, double x_max, double x_min, double epsabs, double epsfrac)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    double r;
    double r_min = x_min;
    double r_max = x_max;

    gsl_function F;
    // struct param_initialConditions params = {this, phi0, dV0, d2V0, phi_absMin, delta_phi_cutoff};

    F.function = func;
    F.params = params;

    gsl_root_fsolver_set(s, &F, r_min, r_max);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r      = gsl_root_fsolver_root(s);
        r_min  = gsl_root_fsolver_x_lower(s);
        r_max  = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(r_min,r_max,epsabs,epsfrac);
    } while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fsolver_free(s);
    return r;
}

double find_min_arg_gsl_wraper(root_func func, void *params, double x_max, double x_min, double epsabs, double epsfrac)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);

    double low = x_min, high = x_max;
    double pred = (low + high)/2;

    gsl_function F;
    F.function = func;
    F.params = params;

    gsl_min_fminimizer_set(s, &F, pred, low, high);

    do
    {
        iter++;
        status = gsl_min_fminimizer_iterate(s);
        pred = gsl_min_fminimizer_x_minimum(s);
        low = gsl_min_fminimizer_x_lower(s);
        high = gsl_min_fminimizer_x_upper(s);

        status = gsl_min_test_interval(low,high,epsabs,epsfrac);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_min_fminimizer_free(s);
    return pred;
}

VD find_multimin_arg_gsl_wraper(multimin_func func, void *params, VD x_start, double epsabs)
{
    int status;
    double size;
    size_t iter = 0, max_iter = 100;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *x;
    gsl_vector *step_size;

    int NDim = x_start.size();

    x = gsl_vector_alloc(NDim);
    step_size = gsl_vector_alloc(NDim);
    for (int i = 0; i < NDim; i++)
    {
        gsl_vector_set(x,i,x_start[i]);
        gsl_vector_set(step_size,i,1.0); // The initial step size is 1.0
    }

    gsl_multimin_function minex_func;
    minex_func.n = NDim;
    minex_func.f = func;
    minex_func.params = params;

    s = gsl_multimin_fminimizer_alloc(T, NDim);
    gsl_multimin_fminimizer_set(s, &minex_func, x, step_size);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status)
        {
            break;
        }
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size,epsabs);

    } while (status == GSL_CONTINUE && iter < max_iter);
    
    VD res(NDim,0);

    if (status == GSL_SUCCESS || status == GSL_ENOPROG)
    {
        for (int i = 0; i < NDim; i++)
        {
            res[i] = gsl_vector_get(s->x,i);
        }
    }

    gsl_vector_free(x);
    gsl_vector_free(step_size);
    gsl_multimin_fminimizer_free(s);
    return res;

}

GSL_Spline_Inter::GSL_Spline_Inter()
{
    _acc = gsl_interp_accel_alloc();
}
GSL_Spline_Inter::~GSL_Spline_Inter()
{
    gsl_interp_accel_free(_acc);
    if (!_spline)
    {
        gsl_spline_free(_spline);
    }
}
void GSL_Spline_Inter::SetData(VD *Y, VD *X)
{
    _Y = Y;
    _X = X;
    _xmin = X->front();
    _xmax = X->back();
    _ymin = Y->front();
    _ymax = Y->back();
    _yaver = (_ymin + _ymax)/2.0;
    _size = X->size();
    if (!_spline)
    {
        gsl_spline_free(_spline);
    }
    if (_size >= 3)
    {
        _spline = gsl_spline_alloc(gsl_interp_cspline,_size);
    }
    else if (_size >= 2)
    {
        _spline = gsl_spline_alloc(gsl_interp_linear,_size);
    }
    gsl_spline_init(_spline,_X->data(),_Y->data(),_size);
}
double GSL_Spline_Inter::valAt(double xi, int deri)
{
    if (_size < 2)
    {
        return _yaver;
    }
    // ! GSL does not support extrapolation, so we extrapolate it using the simplest method (mirror about the end point) by ourself
    // ! But we don't extrapolate too much, otherwise, the user should provide more points to do interpolate
    double x_pick = xi;
    if (xi < _xmin)
    {
        x_pick = _xmin + (_xmin - xi);
        x_pick = x_pick > _xmax?_xmax:x_pick;
    }
    if (xi > _xmax)
    {
        x_pick = _xmax - (xi-_xmax);
        x_pick = x_pick < _xmin?_xmin:x_pick;
    }
    double res;
    switch (deri)
    {
        case 0:
            res = gsl_spline_eval(_spline,x_pick,_acc);
            break;
        case 1:
            return gsl_spline_eval_deriv(_spline,x_pick,_acc);
        case 2:
            res = gsl_spline_eval_deriv2(_spline,x_pick,_acc);
            if (xi < _xmin || xi > _xmax)
            {
                return -res;
            }
            else
            {
                return res;
            }   
        default:
            cout<<"[Warning: GSL_Spline_Inter::valAt] Wrong derivative order ("<<deri<<"), which should be <=2. Use 0 instead"<<endl;
            res = gsl_spline_eval(_spline,x_pick,_acc);
            break;
    }
    if (xi < _xmin)
    {
        res = _ymin - (res - _ymin);
    }
    if (xi > _xmax)
    {
        res = _ymax + (_ymax - res);
    }
    return res;
}
VD GSL_Spline_Inter::valAt(VD X)
{
    VD res;
    for (size_t i = 0; i < X.size(); i++)
    {
        res.push_back(valAt(X[i]));
    }
    return res;
}

GSL_Multi_Spline_Inter::GSL_Multi_Spline_Inter()
{
    _acc = gsl_interp_accel_alloc();
}
GSL_Multi_Spline_Inter::~GSL_Multi_Spline_Inter()
{
    gsl_interp_accel_free(_acc);
    for (size_t i = 0; i < _splines.size(); i++)
    {
        gsl_spline_free(_splines.at(i));
    }
    _splines.clear();
}
void GSL_Multi_Spline_Inter::SetData(VVD *Y, VD *X)
{
    _Y = Y;
    _X = X;
    _size = X->size();
    _NDim = (Y->at(0)).size();
    _xmin = X->front();
    _xmax = X->back();
    _ymin = Y->front();
    _ymax = Y->back();
    _yaver = (_ymin + _ymax)/2.0;
    _Y_transposed = transpose(*Y);
    if (_splines.size()!=0)
    {
        for (size_t i = 0; i < _splines.size(); i++)
        {
            gsl_spline_free(_splines.at(i));
        }
        _splines.clear();
    }

    for (size_t i = 0; i < _NDim; i++)
    {
        if (_size >= 3)
        {
            _splines.push_back(gsl_spline_alloc(gsl_interp_cspline,_size));

        }
        else if (_size >= 2)
        {
            _splines.push_back(gsl_spline_alloc(gsl_interp_linear,_size));
        }
        gsl_spline_init(_splines[i],_X->data(),_Y_transposed[i].data(),_size);
    }
}

VD GSL_Multi_Spline_Inter::valAt(double xi, int deri)
{
    if (_size < 2)
    {
        return _yaver;
    }
    // ! GSL does not support extrapolation, so we extrapolate it using the simplest method (mirror about the end point) by ourself
    // ! But we don't extrapolate too much, otherwise, the user should provide more points to do interpolate
    double x_pick = xi;
    if (xi < _xmin)
    {
        x_pick = _xmin + (_xmin - xi);
        x_pick = x_pick > _xmax?_xmax:x_pick;
    }
    if (xi > _xmax)
    {
        x_pick = _xmax - (xi-_xmax);
        x_pick = x_pick < _xmin?_xmin:x_pick;
    }
    VD res(_NDim);
    switch (deri)
    {
        case 0:
            for (int i = 0; i < _NDim; i++)
            {
                res[i] = gsl_spline_eval(_splines[i],x_pick,_acc);
            }
            break;
        case 1:
            for (int i = 0; i < _NDim; i++)
            {
                res[i] = gsl_spline_eval_deriv(_splines[i],x_pick,_acc);
            }
            return res;
        case 2:
            for (int i = 0; i < _NDim; i++)
            {
                res[i] = gsl_spline_eval_deriv2(_splines[i],x_pick,_acc);
            }
            if (xi<_xmin || xi > _xmax)
            {
                return -res;
            }
            else
            {
                return res;
            }
        default:
            cout<<"[Warning: GSL_Spline_Inter::valAt] Wrong derivative order ("<<deri<<"), which should be <=2. Use 0 instead"<<endl;
            for (int i = 0; i < _NDim; i++)
            {
                res[i] = gsl_spline_eval(_splines[i],x_pick,_acc);
            }
            break;
    }
    if (xi < _xmin)
    {
        res = _ymin - (res - _ymin);
    }
    if (xi > _xmax)
    {
        res = _ymax + (_ymax - res);
    }
    return res;
}
GSL_BSpline_Fit::GSL_BSpline_Fit(int k, int ncoeffs)
{
    _K = k;
    _NCOEFFS = ncoeffs;
    _NBREAKS = _NCOEFFS + 2 - _K;

    _bw = gsl_bspline_alloc(_K,_NBREAKS);
    _B = gsl_vector_alloc(_NCOEFFS);
    _dB = gsl_matrix_alloc(_NCOEFFS,3);
}
GSL_BSpline_Fit::GSL_BSpline_Fit(VVD Y, VD X, int k, int ncoeffs)
{
    _K = k;
    _NCOEFFS = ncoeffs;
    _NBREAKS = _NCOEFFS + 2 - _K;

    _bw = gsl_bspline_alloc(_K,_NBREAKS);
    _B = gsl_vector_alloc(_NCOEFFS);
    _dB = gsl_matrix_alloc(_NCOEFFS,3);

    SetDataX(X);
    UpdateDataY(Y);
}
GSL_BSpline_Fit::~GSL_BSpline_Fit()
{
    gsl_bspline_free(_bw);
    gsl_vector_free(_B);
    gsl_matrix_free(_dB);
    if(!_XC) gsl_matrix_free(_XC);
    if(!_mw) gsl_multifit_linear_free(_mw);
    for (size_t i = 0; i < _Cs.size(); i++)
    {
        gsl_vector_free(_Cs[i]);
    }
    _Cs.clear();
    for (size_t i = 0; i < _COVs.size(); i++)
    {
        gsl_matrix_free(_COVs[i]);
    }
    _COVs.clear();
}
void GSL_BSpline_Fit::SetDataX(VD X)
{
    _X = X;
    _NDataPoints = X.size();
    if (!_XC)
    {
        gsl_matrix_free(_XC);
    }
    _XC = gsl_matrix_alloc(_NDataPoints,_NCOEFFS);

    auto res = minmax_element(_X.begin(),_X.end());
    // cout<<"min: "<<*res.first<<"  max: "<<*res.second<<endl;
    gsl_bspline_knots_uniform(*res.first,*res.second,_bw);

    for (size_t i = 0; i < _NDataPoints; i++)
    {
        gsl_bspline_eval(_X[i],_B,_bw);
        for (size_t j = 0; j < _NCOEFFS; j++)
        {
            gsl_matrix_set(_XC,i,j,gsl_vector_get(_B,j));
        }
    }
    if (!_mw)
    {
        gsl_multifit_linear_free(_mw);
    }
    _mw = gsl_multifit_linear_alloc(_NDataPoints,_NCOEFFS);
}
void GSL_BSpline_Fit::UpdateDataY(VVD Y)
{
    _Y = Y;
    _Y_transposed = transpose(Y);
    _NFields = _Y_transposed.size();

    for (size_t i = 0; i < _Cs.size(); i++)
    {
        gsl_vector_free(_Cs[i]);
    }
    _Cs.clear();
    for (size_t i = 0; i < _COVs.size(); i++)
    {
        gsl_matrix_free(_COVs[i]);
    }
    _COVs.clear();
    Fitting();
}
void GSL_BSpline_Fit::Fitting()
{
    double chisq;
    for (size_t i = 0; i < _NFields; i++)
    {
        _Cs.push_back(gsl_vector_alloc(_NCOEFFS));
        _COVs.push_back(gsl_matrix_alloc(_NCOEFFS,_NCOEFFS));
        gsl_vector_view gv = gsl_vector_view_array(_Y_transposed[i].data(),_NDataPoints);
        gsl_multifit_linear(_XC,&gv.vector,_Cs[i],_COVs[i],&chisq,_mw);
    }
}
VD GSL_BSpline_Fit::valAt(double x)
{
    gsl_bspline_eval(x, _B, _bw);
    VD res(_NFields);
    double yerr;
    for (size_t i = 0; i < _NFields; i++)
    {
        gsl_multifit_linear_est(_B,_Cs[i],_COVs[i],&res[i],&yerr);
    }
    return res;
}
VVD GSL_BSpline_Fit::valAt(VD X)
{
    VVD res;
    for (size_t i = 0; i < X.size(); i++)
    {
        res.push_back(valAt(X[i]));
    }
    return res;
}
tuple<VD,VD> GSL_BSpline_Fit::derivAt(double x)
{
    gsl_vector_view dBi;
    gsl_bspline_deriv_eval(x,2,_dB,_bw);
    VD dY(_NFields);
    VD d2Y(_NFields);
    double yerr;
    for (size_t i = 0; i < _NFields; i++)
    {
        dBi = gsl_matrix_column(_dB,1);
        gsl_multifit_linear_est(&dBi.vector,_Cs[i],_COVs[i],&dY[i],&yerr);
        dBi = gsl_matrix_column(_dB,2);
        gsl_multifit_linear_est(&dBi.vector,_Cs[i],_COVs[i],&d2Y[i],&yerr);
    }
    return make_tuple(dY,d2Y);
}
tuple<VVD,VVD> GSL_BSpline_Fit::derivAt(VD X)
{
    VVD dY;
    VVD d2Y;
    for (size_t i = 0; i < X.size(); i++)
    {
        VD dYi,d2Yi;
        tie(dYi,d2Yi)=derivAt(X[i]);
        dY.push_back(dYi);
        d2Y.push_back(d2Yi);
    }
    return make_tuple(dY,d2Y);
}
