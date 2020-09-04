#include "PathDeformation.h"
#include "RungeKutta.h"
#include <algorithm>

using namespace std;

VVD _pathDeriv(VVD pts)
{
    VVD res;
    int N = pts.size();
    // cout<<"In _pathDeriv: N="<<N<<endl;
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
        for (int i = 0; i < N; i++)
        {
            // cout<<"i="<<i<<endl;
            res[i] = pts[N-1]-pts[0];
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
    // cout<<"N_pts: "<<N_pts<<endl;
    // * 1. Find derivatives along the path formed by pts. (Such that later we can determine the length of the path)
    VVD dpts = _pathDeriv(pts);
    // for (size_t i = 0; i < dpts.size(); i++)
    // {
    //     cout<<"("<<dpts[i]<<")"<<endl;
    // }
    
    // * 2. Extend the path to minima:
    VD xmin_v;
    VD minima;
    double xmin;
    int nx;
    double dx;
    // cout<<""
    if (extend_to_minima)
    {
        // * Extend the front of the path
        _param_extend_to_minima front = {pts[0], dpts[0], V_, T};
        // cout<<"Start finding minimum at front"<<endl;
        // cout<<"From: p=("<<pts[0]<<") with dp=("<<dpts[0]<<")"<<endl;
        xmin_v = find_multimin_arg_gsl_wraper(_func_extend_to_minima,&front,{0.0},1e-6);
        xmin = xmin_v[0]>0.0?0.0:xmin_v[0];
        minima = pts[0]+xmin*dpts[0];
        // cout<<"Need to extend to "<<xmin<<endl;
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
        // cout<<"Start finding minimum at end"<<endl;
        // cout<<"From: p=("<<pts.back()<<") with dp=("<<dpts.back()<<")"<<endl;
        xmin_v = find_multimin_arg_gsl_wraper(_func_extend_to_minima,&end,{0.0},1e-6);
        xmin = xmin_v[0]<0.0?0.0:xmin_v[0];
        minima = pts.back()+xmin*dpts.back();
        // cout<<"Need to extend to "<<xmin<<endl;
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
VVD SplinePath::pts_at_dist(VD X)
{
    VVD res;
    for (size_t i = 0; i < X.size(); i++)
    {
        res.push_back(pts_at_dist(X[i]));
    }
    return res;
}
Deformation_Spline::Deformation_Spline(VVD phi, VD dphidr, dScalarFunction dV, int ncoeffs, int K, double v2min, bool fix_start, bool fix_end, bool save_all_steps):_fitter(K,ncoeffs)
{
    VVD dphi = _pathDeriv(phi);
    VD p_dist = cumtrapz(pow(dphi*dphi,0.5));
    _L = p_dist.back();
    _t = p_dist/_L;

    // VVD phi_lin = phi.front() + (phi.back()-phi.front())*_t;
    _fitter.SetDataX(_t);
    _fitter.UpdateDataY(phi);

    _phi = phi;
    _dV = dV;
    _v2 = pow(dphidr,2);
    _fix_start = fix_start;
    _fix_end = fix_end;
    _save_all_steps = save_all_steps;
    _num_steps = 0;
}
Deformation_Spline::~Deformation_Spline()
{
}
tuple<VVD,VVD> Deformation_Spline::forces()
{
    VVD phi = _phi;
    VD t = _t;
    VVD dphi;
    VVD d2phi;
    for (size_t i = 0; i < t.size(); i++)
    {
        VD _dphi,_d2phi;
        tie(_dphi,_d2phi) = _fitter.derivAt(t[i]);
        dphi.push_back(_dphi);
        d2phi.push_back(_d2phi);
    }
    VD dphi_sq = dphi*dphi;
    VD dphiXd2phi = dphi*d2phi;
    VVD dphids;
    VVD d2phids2;
    for (size_t i = 0; i < t.size(); i++)
    {
        dphids.push_back(dphi[i]/sqrt(dphi_sq[i]));
        d2phids2.push_back(d2phi[i]/dphi_sq[i] - dphi[i]*dphiXd2phi[i]/pow(dphi_sq[i],2));
    }
    VVD dV;
    VVD F_norm;
    VD dV_perp;
    for (size_t i = 0; i < t.size(); i++)
    {
        dV.push_back(_dV(phi[i],nullptr));
        dV_perp = dV[i] - (dV[i]*dphids[i])*dphids[i];
        F_norm.push_back(d2phids2[i]*_v2[i] - dV_perp);
    }
    if (_fix_start)
    {
        F_norm.front() = VD(F_norm.front().size(),0);
    }
    if (_fix_end)
    {
        F_norm.back() = VD(F_norm.back().size(),0);
    }
    return make_tuple(F_norm,dV);
}
tuple<double,bool,double> Deformation_Spline::step(double last_step, double max_step, double min_step, double reverseCheck, double stepIncrease, double stepDecrease, bool checkAfterFit)
{
    VVD F,dV;
    tie(F,dV) = forces();
    VD F_mag = pow(F*F,0.5);
    double F_mag_max = *max_element(F_mag.begin(),F_mag.end());
    VD dV_mag = pow(dV*dV,0.5);
    double dV_mag_max = *max_element(dV_mag.begin(),dV_mag.end());
    double fRatio1 = F_mag_max/dV_mag_max;

    F = F*_L/dV_mag_max;

    double stepsize = last_step;
    bool step_reversed = false;
    VVD phi = _phi;
    if (reverseCheck < 1 && _F_last.size() != 0)
    {
        VD FdotFlast = F*_F_last;
        int reversed = count_if(FdotFlast.begin(),FdotFlast.end(),[](double x){return x<0;});
        if (reversed > (FdotFlast.size())*reverseCheck)
        {
            // Too many points change the sign, we need to reverse to last situation
            // Only if last time, the stepsize is not the minimum one
            // (If it is already minimum stepsize, reverse it helps nothing)
            if (stepsize > min_step)
            {
                step_reversed = true;
                phi = _phi_last;
                F = _F_last;
                stepsize = last_step/stepDecrease;
            }
        }
        else
        {
            // No that many points change the sign, so we want to increase the step size a bit to accelerate the procedure
            stepsize = last_step*stepIncrease;
        }
    }
    if (stepsize > max_step)
    {
        stepsize = max_step;
    }
    if (stepsize < min_step)
    {
        stepsize = min_step;
    }
    
    _phi_last = phi;
    _F_last = F;

    if (_save_all_steps)
    {
        _phi_list.push_back(phi);
        _F_list.push_back(F);
    }

    // * Deform the points by one step
    phi = phi + F*stepsize;

    _fitter.UpdateDataY(phi);
    
    phi = _fitter.valAt(_t);
    _phi = phi;

    VVD Ffit = (phi - _phi_last)/stepsize;
    VD Ffit_mag = pow(Ffit*Ffit,0.5);
    double fRatio2 = *max_element(Ffit_mag.begin(),Ffit_mag.end())/_L;
    
    double fRatio = checkAfterFit?fRatio2:fRatio1;
    return make_tuple(stepsize,step_reversed,fRatio);
}
Deformation_Status Deformation_Spline::deformPath(double start_step,double fRatioConv,double converge_0, double fRatioIncrease, int maxiter, function<bool(Deformation_Spline*)> callback)
{
    double minfRatio = INFINITY;
    int minfRatio_index = 0;
    VVD minfRatio_phi;
    double stepsize = start_step;
    bool step_reversed;
    Deformation_Status deformation_info;
    double fRatio;
    while (true)
    {
        _num_steps += 1;
        tie(stepsize,step_reversed,fRatio) = step(stepsize);
        if (!callback(this))
        {
            deformation_info = FAILCALLBACK;
            break;
        }

        // Check if the deformation converged
        if (fRatio < fRatioConv || (_num_steps == 1 && fRatio < converge_0*fRatioConv))
        {
            deformation_info = DF_CONVERGED;
            break;
        }

        // Update the cache storing the information of minium fRatio
        if (fRatio < minfRatio)
        {
            minfRatio = fRatio;
            minfRatio_index = _num_steps;
            minfRatio_phi = _phi;
        }
        
        // Check if the convergence is getting worse, if so, just abort before thing getting pretty bad.
        if (!step_reversed && fRatio > fRatioIncrease * minfRatio)
        {
            _phi = minfRatio_phi;
            _fitter.UpdateDataY(_phi);
            _phi_list.resize(minfRatio_index);
            _F_list.resize(minfRatio_index);
            deformation_info = DIVERGENT;
            break;
        }   
    }
    return deformation_info;
}
