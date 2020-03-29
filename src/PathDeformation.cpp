#include "PathDeformation.h"
#include "RungeKutta.h"

VVD _pathDeriv(VVD pts)
{
    VVD res;
    int N = pts.size();
    if (N >= 5)
    {
        res = deriv14_const_dx(pts);
    }
    else if (N > 2)
    {
        res = VVD(N);
        res[0] = -1.5*pts[0] + 2.0*pts[1] - 0.5*pts[2];
        for (int i = 1; i < N-1; i++)
        {
            res[i] = (pts[i+1]-pts[i-1])/2.0;
        }
        res[N-1] = 1.5*pts[N-1] - 2.0*pts[N-2] + 0.5*pts[N-3];
    }
    else
    {
        res = VVD(N);
        for (size_t i = 0; i < N; i++)
        {
            res[N] = pts[N-1]-pts[0];
        }
    }
    return res;
}
struct _param_extend_to_minima
{
    VD p0;
    VD dp0;
    ScalarFunction V;
    double T;
};

double _func_extend_to_minima(const gsl_vector *x, void *param)
{
    _param_extend_to_minima *mod = (_param_extend_to_minima*)param;
    VD p0 = mod->p0;
    VD dp0 = mod->dp0;
    double xi = gsl_vector_get(x,0);
    VD px = p0 + xi*dp0;
    return mod->V(px,&(mod->T));
}
VD _func_refine_dist(double x, VD y, void *param)
{
    SplinePath* mod = (SplinePath*)param;
    VD dpdx = mod->pts_at_dist(x,1);
    double d_dist = sqrt(dpdx*dpdx);
    return {d_dist};
}
SplinePath::SplinePath(VVD pts, ScalarFunction V_, double T, int V_spline_samples, bool extend_to_minima, bool re_eval_distances)
{
    int N_pts = pts.size();
    
    // * 1. Find derivatives along the path formed by pts. (Such that later we can determine the length of the path)
    VVD dpts = _pathDeriv(pts);
    
    // * 2. Extend the path to minima:
    VD xmin_v;
    VD minima;
    double xmin;
    int nx;
    double dx;
    if (extend_to_minima)
    {
        // * Extend the front of the path
        _param_extend_to_minima front = {pts[0], dpts[0], V_, T};
        xmin_v = find_multimin_arg_gsl_wraper(_func_extend_to_minima,&front,{0.0},1e-6);
        xmin = xmin_v[0]>0.0?0.0:xmin_v[0];
        minima = pts[0]+xmin*dpts[0];
        nx = ceil(abs(xmin)-0.5);
        if (nx > 0)
        {
            dx = xmin/nx;
            for (size_t i = 0; i < nx; i++)
            {
                VD p = pts[0] + dx*dpts[0];
                pts.insert(pts.begin(),p);   
            } 
            pts[0] = minima; // Force to be at the correct place;
        }

        // * Extend at the end of the path
        _param_extend_to_minima end = {pts.back(),dpts.back(),V_,T};
        xmin_v = find_multimin_arg_gsl_wraper(_func_extend_to_minima,&end,{0.0},1e-6);
        xmin = xmin_v[0]<0.0?0.0:xmin_v[0];
        minima = pts.back()+xmin*dpts.back();
        nx = ceil(abs(xmin)-0.5);
        if (nx > 0)
        {
            dx = xmin/nx;
            for (size_t i = 0; i < nx; i++)
            {
                VD p = pts.back() + dx*dpts.back();
                pts.insert(pts.end(),p);
            }
            pts.back() = minima;
        }
        dpts = _pathDeriv(pts);
    }

    // * 3. Get the spline from distance to field points:
    VD p_dist = cumtrapz(pow(dpts*dpts,0.5));

    _L = p_dist.back(); // The last one is the total length of the path

    // ! Interpolate the pts (as y) vs. p_dist (as x) using spline
    _path_Inter.SetData(&pts,&p_dist);

    // * 4. Re-evaluate the distance 
    if (re_eval_distances)
    {
        RungeKutta _rk_solver(1);
        _rk_solver.SetODE(_func_refine_dist);
        _rk_solver.SetParams(this);
        VVD _dist_tmp = _rk_solver.ODEINTEGRAL(p_dist,{0.0},1e-1,_L*1e-8);
        p_dist = transpose(_dist_tmp)[0];
        _L = p_dist.back();
        _path_Inter.SetData(&pts,&p_dist);
    }

    // * 5. Interpolate the potential and its derivatives
    VD x = linspace(0,_L,V_spline_samples);
    // Extend a little bit at two ends to model the end points more accurately
    VD x_ext = linspace(x[1],_L*0.2,x[1]);
    VD x_end = _L + x_ext;
    x_ext = -x_ext;
    VD x_front(x_ext.rbegin(),x_ext.rend());
    x.insert(x.begin(),x_front.begin(),x_front.end());
    x.insert(x.end(),x_end.begin(),x_end.end());
    VD y;
    for (int i = 0; i < x.size(); i++)
    {
        y.push_back(V_(pts_at_dist(x[i]),&T));
    }
    _V_Inter.SetData(&y,&x);
    V = [&](VD x, void *param){return _V_Inter.valAt(x[0]);};
    dV = [&](VD x, void *param){VD res = {_V_Inter.valAt(x[0],1)}; return res;};
    d2V = [&](VD x, void *param){VVD res; res.push_back({_V_Inter.valAt(x[0],2)}); return res;}; 
    
}

VD SplinePath::pts_at_dist(double x, int deri)
{
    return _path_Inter.valAt(x,deri);
}