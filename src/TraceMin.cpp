/*
 * @Description  : Trace the minima of a potential and find critical temperature
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-29 16:15:53
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-23 17:17:47
 */
#include <iostream>
#include "TraceMin.h"
#include "Phases.h"
#include <tuple>
#include <functional>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>

using namespace std;

/*
 * @description: Convert from VD to gsl_vector_view
 * @param
 *      vec:
 * @return: 
 *      gsl_vector_view
*/
gsl_vector_view Get_GSL_Vector_View(VD &vec)
// * Using reference in order to keep the scope of vec_data from vec, otherwise, the vector_view will loss the value.
{
    int Ndim = vec.size();
    double *vec_data = vec.data();
    return gsl_vector_view_array(vec_data,Ndim);
}
/*
 * @description: Get the minimum/maximum eigenvalue of a matrix
 * @param:
 *      mat: VVD
 *          The matrix which we want to get the eigenvalue
 *      sort:
 *          The gsl_eigen_sort_t sort method that sort the eigenvalue
 * @return: 
 *      eigmin:
 *          The minimum eigenvalue (in the sence of the sort method `sort`)
 *      eigmax:
 *          The maximum eigenvalue
*/
void Get_Matrix_Eigen(VVD mat_vvd, gsl_eigen_sort_t sort, double &eigmin, double &eigmax)
{
    int Ndim = mat_vvd.size();
    double *M0_data = new double[Ndim*Ndim];
    for (int i = 0; i < Ndim; i++)
    {
        for (int j = 0; j < Ndim; j++)
        {
            M0_data[i*Ndim+j] = mat_vvd[i][j];
        }
    }
    gsl_matrix_view mat = gsl_matrix_view_array(M0_data,Ndim,Ndim);
    gsl_vector *eval = gsl_vector_alloc(Ndim);
    gsl_matrix *evec = gsl_matrix_alloc(Ndim,Ndim);//Just used for sorting
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(Ndim);
    gsl_eigen_symmv(&mat.matrix,eval,evec,w);
    gsl_eigen_symmv_sort(eval,evec,sort);
    eigmin = gsl_vector_get(eval,0);
    eigmax = gsl_vector_get(eval,Ndim-1);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_eigen_symmv_free(w);
    delete []M0_data;
}

// * The structure used to pass to following wrap function.
struct fdf_for_min
{
    ScalarFunction f;
    dScalarFunction df_dx;
    double t;
};

/*
 * @description: The wrap function that need to be minimized. The true function that is minimized (as well as other parameters) is passed through param, which is of struct type fdf_for_min.
 * @param:
 *      x:
 *          The argument of the function
 *      param:
 *          The struct fdf_for_min
 * @return: 
 *      The value of the function @ x and other parameters in `param`
 */
double f_min_wrap(const gsl_vector *x, void *param)
{
    fdf_for_min *pp = (fdf_for_min *)param;
    int Ndim = x->size;
    double t = pp->t;
    VD X(Ndim);
    for (int i = 0; i < Ndim; i++)
    {
        X[i] = gsl_vector_get(x,i);
    }
    return (pp->f)(X,t);
}
void df_min_wrap(const gsl_vector *x, void *param, gsl_vector *g)
{
    fdf_for_min *pp = (fdf_for_min *)param;
    int Ndim = x->size;
    double t = pp->t;
    VD X(Ndim);
    for (int i = 0; i < Ndim; i++)
    {
        X[i] = gsl_vector_get(x,i);
    }
    VD res = (pp->df_dx)(X,t);
    for (int i = 0; i < Ndim; i++)
    {
        gsl_vector_set(g,i,res[i]);
    }
}
void fdf_min_wrap(const gsl_vector *x, void *param, double *f, gsl_vector *g)
{
    *f = f_min_wrap(x,param);
    df_min_wrap(x,param,g);
}
/*
 * @description: This function trace one single local minimum
 * @param 
 *      f: 
 *          The scalar potential function `f(x,t)` which needs to be
 *          minimized. The input will be of the same type as `(x0,t0)`
 *      df_dx, d2f_dxdt, d2f_dx2:
 *          Functions which return derivatives of `f(x)`.
 *          `df_dx` return the derivative of `f(x)` repect to `x`
 *          `d2f_dxdt` return the derivative of the gradient of 
 *          `f(x)` with respect to `t`
 *          `d2f_dx2` return the Hessian matrix of `f(x)`
 *      x0:
 *          The initial value for the fields
 *      t0:
 *          The initial value for the temperature
 *      tstop:
 *          Stop the trace when `t` reaches `tstop`
 *      dtstart:
 *          Initial stepsize in `t`
 *      deltaX_target:
 *          The target error in x at each step. Determines the
 *          stepsize in t by extrapolation from last error.
 *      dtabsMax: default = 20
 *      dtfracMax: default = 0.25
 *          The largest stepsize in t will be the LARGEST of
 *          ``abs(dtstart)*dtabsMax`` and ``t*dtfracMax``.
 *      dtmin: default = 1e-3
 *          The smallest stepsize we'll allow before assuming the
 *          transition ends,
 *          relative to `dtstart`
 *      deltaX_tol: default = 1.2
 *          ``deltaX_tol*deltaX_target`` gives the maximum error in x
 *          before we want to shrink the stepsize and recalculate the
 *          minimum
 *      minratio: default = 1e-4
 *          The smallest ratio between smallest and largest eigenvalues
 *          in the Hessian Matrix before treating the smallest
 *          eigenvalue as zero (and thus signaling a saddle point and
 *          the end of the minimum)
 * @return: _traceMinimum_rval structure 
 *      X, T, dXdT:
 *          Arrays of the minimum at different values of t,
 *          and its derivative with respect to t.
 *      overX:
 *          The point beyond which the phase seems to disappear
 *      overT:
 *          The t-value beyond which the phase seems to disappear
*/
_traceMinimum_rval traceMinimum(ScalarFunction f, dScalarFunction df_dx, dScalarFunction d2f_dxdt, HM d2f_dx2, VD x0, double t0, double tstop, double dtstart/*, double deltaX_target, double dtabsMax, double dtfracMax, double dtmin, double deltaX_tol, double minratio*/)
{
#if VERBOSE == 1
    printf("traceMinimum t0 = %.6f\n",t0);
#endif
    int Ndim = x0.size();
    VVD M0_VVD = d2f_dx2(x0,t0);
    double min_abs_eigen, max_abs_eigen;
    Get_Matrix_Eigen(M0_VVD,GSL_EIGEN_SORT_ABS_ASC,min_abs_eigen,max_abs_eigen);

    // This determine when we treat the matrix is singular
    double minratio = abs(min_abs_eigen)/abs(max_abs_eigen) * pre_control.minratio_rel;

    // This function determine the dx/dt at minimum point, and also check whether we encounter saddle/maximum point
    function<tuple<VD, bool>(VD,double) > dxmindt = [&](VD x, double t){
        VVD M1_VVD = d2f_dx2(x,t);
        double *M0_data = new double[Ndim*Ndim];
        for (int i = 0; i < Ndim; i++)
        {
            for (int j = 0; j < Ndim; j++)
            {
                M0_data[i*Ndim+j] = M1_VVD[i][j];
            }
        }
        gsl_matrix_view M = gsl_matrix_view_array(M0_data,Ndim,Ndim);
        VD b_VD = d2f_dxdt(x,t);
        gsl_vector_view b = Get_GSL_Vector_View(b_VD);
        
        gsl_vector *dxdt_tmp = gsl_vector_alloc (Ndim);
        int s;
        gsl_permutation * p = gsl_permutation_alloc (Ndim);
        gsl_linalg_LU_decomp (&M.matrix, p, &s);
        gsl_linalg_LU_solve (&M.matrix, p, &b.vector, dxdt_tmp);
        double min_eigen, max_eigen;
        Get_Matrix_Eigen(M1_VVD,GSL_EIGEN_SORT_VAL_ASC,min_eigen,max_eigen);
        bool isneg = (min_eigen<0 || min_eigen/max_eigen < minratio);
        VD dxdt(Ndim);
        for (int i = 0; i < Ndim; i++)
        {
            dxdt[i] = -gsl_vector_get(dxdt_tmp,i);
        }
        gsl_vector_free(dxdt_tmp);
        gsl_permutation_free(p);
        delete []M0_data;
        return make_tuple(dxdt, isneg);
    };

    const double xeps = pre_control.deltaX_target*1e-2;
    function<VD(VD, double)> fmin = [&](VD x, double t){
        const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
        gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, Ndim);
        gsl_multimin_function_fdf minex_func;
        gsl_vector_view X = gsl_vector_view_array(x.data(),Ndim);
        minex_func.n = Ndim;
        minex_func.f = f_min_wrap;
        minex_func.df = df_min_wrap;
        minex_func.fdf = fdf_min_wrap;
        fdf_for_min par = {.f = f, .df_dx = df_dx, .t = t};
        minex_func.params = (void *)&par;
        gsl_multimin_fdfminimizer_set(s, &minex_func, &X.vector, xeps*1e-1,xeps);
        int iter = 0;
        int status;
        do
        {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate(s);
            if (status)
            {
                break;
            }
            status = gsl_multimin_test_gradient(s->gradient, xeps);
            
        } while (status == GSL_CONTINUE && iter < 100);
        
        VD res(x);
        if (status == GSL_SUCCESS || status == GSL_ENOPROG)
        {
            #if VERBOSE == 2
            cout<<"status: "<<status<<"; [";
            #endif
            for (int i = 0; i < Ndim; i++)
            {
                res[i] = gsl_vector_get(s->x,i);
                #if VERBOSE == 2
                cout<<res[i]<<",";
                #endif
            }
            #if VERBOSE == 2
            cout<<"]"<<endl;
            #endif
        }

        gsl_multimin_fdfminimizer_free(s);
        return res;
    };

//  1.2 usr_set   1.2       *  usr_set
    const double deltaX_tol = pre_control.deltaX_tol * pre_control.deltaX_target;
//  initial step in t
    const double tscale = abs(dtstart);
//                          20     * initial step in t
    const double dtabsMax = pre_control.dtabsMax*tscale;

//          1e-3  * initial step in t
    const double dtmin = pre_control.dtmin * tscale;

    VD x = fmin(x0,t0); // Re-find the minimum point to make sure it is the local minimum.
    double t = t0;
    double dt = dtstart;
    double xerr = 0.0;
    VD dxdt;
    bool negeig;
    tie(dxdt,negeig) = dxmindt(x,t);

    VVD X; // Store the minimum points
    VD T; // Store the temperature at each point
    VVD dXdT; // Store the dX/dT at each point

    // ! push in the first point
    X.push_back(x);
    T.push_back(t);
    dXdT.push_back(dxdt);
    VD overX;
    double overT = -11;

    double tnext;
    VD xnext;
    VD dxdt_next;
    // int index=0;
    #if VERBOSE == 2
    int index = 0;
    cout<<">>The tracing History:"<<endl;
    cout<<">> Starting at: T="<<t<<", [";
    for (size_t i = 0; i < x.size(); i++)
    {
        cout<<x[i]<<",";
    }
    cout<<"]"<<endl;
    #endif
    while (true)
    {
        // index++;
        #if VERBOSE == 1
        cout<<".";//<<endl;
        cout.flush();
        #endif
        // ! Get the values at the next step
        tnext = t+dt;
        xnext = fmin(x+dxdt*dt,tnext);
        tie(dxdt_next, negeig) = dxmindt(xnext,tnext);
        #if VERBOSE == 2
        index++;
        cout<<">>>>Step-"<<index<<endl;
        cout<<">>>>>> T="<<tnext<<" ,X=[";
        for (size_t itemp = 0; itemp < xnext.size(); itemp++)
        {
            cout<<xnext[itemp]<<",";
        }
        cout<<"], dXdT=[";
        for (size_t itemp = 0; itemp < xnext.size(); itemp++)
        {
            cout<<dxdt_next[itemp]<<",";
        }
        cout<<"], negeig="<<negeig<<endl;
        #endif
        if (negeig)
        {
        // ! If negeig is true from dxmindt, that means, (xnext, tnext) already become saddle/maximum point, so we need to shrink the step and try again. And at most, at tnext and xnext, this phase will disappear, becoming saddle/maximum point
            dt *= 0.5;
            overX = xnext;
            overT = tnext;
        }
        else
        {
            // ! The step might still be too big if it's outside of our error tolerance.
            xerr = sqrt(max((x+dxdt*dt - xnext)*(x+dxdt*dt - xnext),(xnext-dxdt_next*dt - x)*(xnext-dxdt_next*dt - x)));
            if (xerr < deltaX_tol)
            {
                // ! The step is good, push back the point we got.
                T.push_back(tnext);
                X.push_back(xnext);
                dXdT.push_back(dxdt_next);
                if (overT < -10)
                {
                    // ! If we haven't yet encounter a saddle point, we can change the step size to accelerate the trace
                    // ! The step size is changed according to the current xerr between guessed point and true (from fmin) point and the precision in X we want to achieve.
                    dt *= pre_control.deltaX_target/(xerr+1e-30);
                }
                x = xnext;
                t = tnext;
                dxdt = dxdt_next;
                overX.clear();
                overT = -11;  
            }
            else
            {
                // ! Either stepsize was too big, or we hit a transition. Just cut the step in half.
                dt *= 0.5;
                overX = xnext;
                overT = tnext;
            }   
        }
        if (abs(dt) < abs(dtmin))
        {
            break;
        }
        if ((dt > 0 && t >= tstop) || (dt < 0 && t <= tstop))
        {
            dt = tstop-t;
            x = fmin(x+dxdt*dt, tstop);
            tie(dxdt,negeig) = dxmindt(x,tstop);
            t = tstop;
            X[X.size()-1] = x;
            T[T.size()-1] = t;
            dXdT[dXdT.size()-1] = dxdt;
            overX.clear();
            overT = -11; //! As we hit the stop point, no over point, or you can also say the over point is the last point, which will be set after the loop in the following.
            break;
        }
        double dtmax = max(t*pre_control.dtfracMax, dtabsMax);
        if (abs(dt) > dtmax)
        {
            dt = abs(dt)/dt*dtmax;
        }
    }
    if (overT < -10)
    {
        overX = X[X.size()-1];
        overT = T[T.size()-1];
    }
    #if VERBOSE == 1
    cout<<endl;
    #endif
    return {X, T, dXdT, overX, overT}; 
}

/*
 * @description: This function tries to find all phases for one potential
 * @param
 *      f:
 *          The potential that we want to find all phases
 *      df_dx, d2f_dxdt, d2f_dx2:
 *          The derivatives of the potential
 *      points:
 *          Starting points (can be approximate local minima)
 *      tLow, tHigh:
 *          The lowest and highest temperature we want to trace
 *      deltaX_target:
 *          The tolerance in minimization
 *      dtstart: default 1e-3
 *          The starting step size;
 *      tjump: default 1e-3
 *          The jump step size in `t` from the end of one phase 
 *          to another initial tracing point
 *      fc: default AllPass
 *          The function determine whether we want to forbid a phase
 * @return: 
 *      MP: std::map<int,Phase>
 *          A map from key => Phase
 *          The key is unique for each phase.
 */
MP traceMultiMin(ScalarFunction f, dScalarFunction df_dx, dScalarFunction d2f_dxdt, HM d2f_dx2, VT points, double tLow, double tHigh, /*double deltaX_target, double dtstart, double tjump,*/ forbidCrit fc)
{
    MP phases;
    if (points.size()==0)
    {
        // ! No point provide, nothing to be started with
        return phases;
    }
    MP::iterator iter;
    int Ndim = get<0>(points[0]).size();
    const double xeps = pre_control.deltaX_target*1e-2;
    function<VD(VD, double)> fmin = [&](VD xin, double t){
        const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
        gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, Ndim);
        gsl_multimin_function_fdf minex_func;
        VD x_start = xin+xeps;
        // VD x_start = xin + 0.1;
        gsl_vector_view X = gsl_vector_view_array(x_start.data(),Ndim);
        minex_func.n = Ndim;
        minex_func.f = f_min_wrap;
        minex_func.df = df_min_wrap;
        minex_func.fdf = fdf_min_wrap;
        fdf_for_min par = {.f = f, .df_dx = df_dx, .t = t};
        minex_func.params = (void *)&par;
        gsl_multimin_fdfminimizer_set(s, &minex_func, &X.vector, xeps*0.1,xeps);
        int iter = 0;
        int status;
        do
        {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate(s);
            if (status)
            {
                break;
            }
            status = gsl_multimin_test_gradient(s->gradient, xeps*1e-3);
            
        } while (status == GSL_CONTINUE && iter < 100);
        
        VD res(Ndim,0);
        if (status == GSL_SUCCESS || status == GSL_ENOPROG)
        {
            for (int i = 0; i < Ndim; i++)
            {
                res[i] = gsl_vector_get(s->x,i);
            }
        }
        
        gsl_multimin_fdfminimizer_free(s);
        return res;
    };

    const double dtstart = pre_control.dtstart * (tHigh-tLow);
    const double tjump = pre_control.tjump * (tHigh-tLow);
    VnT nextPoint;
    VD x;
    double t;
    for (int i = 0; i < points.size(); i++)
    {
        tie(x,t) = points[i];
        nextPoint.push_back(make_tuple(t, dtstart, fmin(x,t),-1));
    }
    
    VD x1;
    double t1, dt1;
    int linkedFrom;
    int round = 0;
    while (nextPoint.size()!=0)
    {
        #if VERBOSE == 2
        round++;
        cout<<"At round: "<<round<<endl;
        for (size_t i = 0; i < nextPoint.size(); i++)
        {
            tie(t1,dt1,x1,linkedFrom) = nextPoint[i];
            cout<<"Point-"<<i<<" T="<<t1<<", X=[";
            for (size_t j = 0; j < x1.size(); j++)
            {
                cout<<x1[j]<<",";
            }
            cout<<"]"<<endl;
        }
        #endif

        tie(t1,dt1,x1,linkedFrom) = nextPoint.back();
        nextPoint.pop_back();
        
        x1 = fmin(x1,t1); // make sure we start as accurately as possible

        // * Check the bounds
        if ((t1 < tLow) || (t1 == tLow && dt1 < 0))
        {
            continue;
        }
        if ((t1 > tHigh) || (t1 == tHigh and dt1 > 0))
        {
            continue;
        }
        if (fc(x1))
        {
            continue;
        }
        
        // * Check if it is redundant with another phase
        int index;
        Phase *ph;
        bool covered = false;
        for (iter = phases.begin(); iter != phases.end(); ++iter)
        {
            index = iter->first;
            ph = &(iter->second);
            if ((t1 < ph->GetTmin()) || t1 > ph->GetTmax())
            {
                continue;
            }
            x = fmin(ph->valAt(t1),t1);
            if (sqrt((x-x1)*(x-x1)) < 2*pre_control.deltaX_target)
            {
                // * This point is already covered
                if (linkedFrom != index && linkedFrom != -1)
                {
                    ph->addLinkFrom(&(phases.at(linkedFrom)));
                }
                covered = true;
                break;
            }
        }
        _traceMinimum_rval down_trace, up_trace;
        VVD X_down, X_up, dXdT_down, dXdT_up, X, dXdT;
        VD T_down, T_up, nX, T;
        double nT;
        VD x2;
        double t2, dt2;
        VVD points;
        if (!covered)
        {
            // * Not covered, Trace the phase
            #if VERBOSE == 1
            cout<<"...Tracing phase starting at x = (";
            for (int ix = 0; ix < x1.size(); ix++)
            {
                cout<<" "<<x1[ix]<<" ";
            }
            cout<<"); t = "<<t1<<" ..."<<endl;
            #endif
            int phase_key = phases.size();
            int oldNumPoints = nextPoint.size();
            if (t1 > tLow)
            {
                #if VERBOSE == 1
                cout<<"......Tracing minimum down to "<<tLow<<" ......"<<endl;
                #endif
                down_trace = traceMinimum(f,df_dx,d2f_dxdt, d2f_dx2, x1, t1, tLow, -dt1);
                X_down = down_trace.X;
                T_down = down_trace.T;
                dXdT_down = down_trace.dXdT;
                nX = down_trace.overX;
                nT = down_trace.overT;
                t2 = nT - tjump;
                dt2 = 0.1*tjump;
                x2 = fmin(nX,t2);
                nextPoint.push_back(make_tuple(t2,dt2,x2,phase_key));
                if (sqrt((X_down.back()-x2)*(X_down.back()-x2))>pow(pre_control.deltaX_target,2))
                {
                    points = findApproxLocalMin(f,X_down.back(),x2,t2);
                    for (int ipoints = 0; ipoints < points.size(); ipoints++)
                    {
                        nextPoint.push_back(make_tuple(t2,dt2,fmin(points[ipoints],t2),phase_key));
                    }
                }
                reverse(X_down.begin(),X_down.end());
                reverse(T_down.begin(),T_down.end());
                reverse(dXdT_down.begin(),dXdT_down.end());
            }
            if (t1 < tHigh)
            {
                #if VERBOSE == 1
                cout<<"......Tracing minimum up to "<<tHigh<<" ......"<<endl;
                #endif
                up_trace = traceMinimum(f,df_dx,d2f_dxdt,d2f_dx2,x1,t1,tHigh,dt1);
                X_up = up_trace.X;
                T_up = up_trace.T;
                dXdT_up = up_trace.dXdT;
                nX = up_trace.overX;
                nT = up_trace.overT;
                t2 = nT+tjump;
                dt2 = 0.1*tjump;
                x2 = fmin(nX,t2);
                nextPoint.push_back(make_tuple(t2,dt2,x2,phase_key));
                if (sqrt((X_up.back()-x2)*(X_up.back()-x2))>pow(pre_control.deltaX_target,2))
                {
                    points = findApproxLocalMin(f,X_up.back(),x2,t2);
                    for (int ipoints = 0; ipoints < points.size(); ipoints++)
                    {
                        nextPoint.push_back(make_tuple(t2,dt2,fmin(points[ipoints],t2),phase_key));
                    }
                }
            }
            if (t1 <= tLow)
            {
                X = X_up;
                T = T_up;
                dXdT = dXdT_up;
            }
            else if (t1 >= tHigh)
            {
                X = X_down;
                T = T_down;
                dXdT = dXdT_down;
            }
            else
            {
                X = VVD(X_down);
                X.insert(X.end(),X_up.begin()+1,X_up.end());
                T = VD(T_down);
                T.insert(T.end(),T_up.begin()+1,T_up.end());
                dXdT = VVD(dXdT_down);
                dXdT.insert(dXdT.end(),dXdT_up.begin()+1,dXdT_up.end());
            }
            if (fc(X.front()) || fc(X.back()))
            {
                nextPoint.resize(oldNumPoints);
            }
            else if (X.size()>1)
            {
                Phase newphase(phase_key,X,T,dXdT);
                if (linkedFrom >= 0)
                {
                    newphase.addLinkFrom(&(phases.at(linkedFrom)));
                }
                phases.insert(pair<int, Phase>(phase_key,newphase));
            }
            else
            {
                nextPoint.resize(oldNumPoints);
            }
        }
    }
    return phases;
}

/*
 * @description: Find minima between two points
 * @param
 *      f:
 *          The function to minimize
 *      x0,x2:
 *          The points between which to find minima
 *      ti:
 *          The temperature suppling to f
 *      n:
 *          Number of points to test for local minima
 *      edge:
 *          Determine whether test for minima directly next to input
 *          points
 *          edge = 0, do such test
 *          edge = 0.5 means only test the center point
 * @return: 
 *      A list of points that are approximate minima
 */
VVD findApproxLocalMin(ScalarFunction f, VD x0, VD x1, double ti, int n, double edge)
{
    if (edge > 0.5)
    {
        edge = 0.5;
    }
    VVD res;
    VVD points;
    bool *test = new bool[n];
    double min = edge;
    double step = (1.0-2*edge)/(n-1);
    for (int i = 0; i < n; i++)
    {
        points.push_back(x0+(x1-x0)*(min+i*step));
        // printf("%d, (%.4f,%.4f), %.5f\n",i,points[i][0],points[i][1],f(points[i],ti));
        test[i] = false;
        if (i!=0 && i!=n-1)
        {
            if (f(points[i],ti)>f(points[i-1],ti))
            {
                test[i] = false;
                test[i-1] *= true;
            }
            else
            {
                test[i] = true;
                test[i-1] *= false;
            } 
            if (test[i-1])
            {
                res.push_back(points[i-1]);
            }
        }
    }
    delete []test;
    return res;
}

/*
 * @description: Remove redundant phases
 * @param 
 *      f:
 *          The potential function which was passed to `traceMultiMin`
 *      phases:
 *          The return value of `traceMultiMin`
 *      xeps: default 1e-5
 *          The Error tolerance in minimization
 *      diftol: default 1e-2
 *          Maximum separation between two phases before 
 *          they are considered to be coincident
 * @return: 
 *      directly modify phases
 */
void removeRedundantPhases(ScalarFunction f, dScalarFunction df_dx, MP &phases, double xeps, double diftol)
{
    // xeps = 0.1;
    function<VD(VD, double)> fmin = [&](VD x, double t){
        int Ndim = x.size();
        const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
        gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, Ndim);
        gsl_multimin_function_fdf minex_func;
        VD x_start = x+xeps;
        // printf("starting at: %.5f, %.5f\n",x_start[0],x_start[1]);
        gsl_vector_view X = gsl_vector_view_array(x_start.data(),Ndim);
        minex_func.n = Ndim;
        minex_func.f = f_min_wrap;
        minex_func.df = df_min_wrap;
        minex_func.fdf = fdf_min_wrap;
        fdf_for_min par = {.f = f, .df_dx = df_dx, .t = t};
        minex_func.params = (void *)&par;
        gsl_multimin_fdfminimizer_set(s, &minex_func, &X.vector, xeps*0.1,xeps);
        int iter = 0;
        int status;
        do
        {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate(s);
            if (status)
            {
                break;
            }
            status = gsl_multimin_test_gradient(s->gradient, xeps);
            // if (status == GSL_SUCCESS)
            //     printf ("Minimum found at:\n");

            // printf ("%5d %.5f %.5f %10.5f\n", iter,
            //   gsl_vector_get (s->x, 0),
            //   gsl_vector_get (s->x, 1),
            //   s->f);
        } while (status == GSL_CONTINUE && iter < 100);
        
        VD res(x);
        if (status == GSL_SUCCESS || status == GSL_ENOPROG)
        {
            for (int i = 0; i < Ndim; i++)
            {
                res[i] = gsl_vector_get(s->x,i);
            }
        }
        gsl_multimin_fdfminimizer_free(s);
        return res;
    };
    bool has_redundant_phase = true;
    MP::iterator iter1;
    MP::iterator iter2;
    int index1,index2;
    int index_low, index_high;
    int index_rej;
    int newkey;
    double tmax,tmin;
    double tmax1,tmax2,tmin1,tmin2;
    VD x1, x2;
    VVD X;
    VD T;
    VVD dXdT;
    double dif;
    bool same_at_tmax, same_at_tmin;
    while (has_redundant_phase)
    {
        has_redundant_phase = false;
        for (iter1 = phases.begin(); iter1 != phases.end(); iter1++)
        {
            index1 = iter1->first;
            for (iter2 = phases.begin(); iter2 != phases.end(); iter2++)
            {
                index2 = iter2->first;
                // cout<<"----"<<index1<<" "<<index2<<endl;
                if (index1 == index2)
                {
                    continue;
                }
                Phase phase1 = phases.at(index1);
                Phase phase2 = phases.at(index2);
                tmin1 = phase1.GetTmin();
                tmax1 = phase1.GetTmax();
                tmin2 = phase2.GetTmin();
                tmax2 = phase2.GetTmax();
                tmax = min(tmax1,tmax2);
                tmin = max(tmin1,tmin2);
                if (tmin > tmax)
                {
                    continue;
                }
                if (tmax == tmax1)
                {
                    x1 = phase1.GetXatTmax();
                }
                else
                {
                    x1 = fmin(phase1.valAt(tmax),tmax);
                }
                // cout<<"Phase 1 xmax:"<<x1[0]<<" "<<x1[1]<<endl;
                if (tmax == tmax2)
                {
                    x2 = phase2.GetXatTmax();
                }
                else
                {
                    x2 = fmin(phase2.valAt(tmax),tmax);
                }
                // cout<<"Phase 2 xmax:"<<x2[0]<<" "<<x2[1]<<endl;
                dif = sqrt((x2-x1)*(x2-x1));
                same_at_tmax = (dif < diftol);
                if (tmin == tmin1)
                {
                    x1 = phase1.GetXatTmin();
                }
                else
                {
                    x1 = fmin(phase1.valAt(tmin),tmin);
                }
                // cout<<"Phase 1 xmin:"<<x1[0]<<" "<<x1[1]<<endl;
                if (tmin == tmin2)
                {
                    x2 = phase2.GetXatTmin();
                }
                else
                {
                    x2 = fmin(phase2.valAt(tmin),tmin);
                }
                // cout<<"Phase 2 xmin:"<<x2[0]<<" "<<x2[1]<<endl;
                dif = sqrt((x2-x1)*(x2-x1));
                same_at_tmin = (dif < diftol);
                if (same_at_tmin && same_at_tmax)
                {
                    // * Phases are redundant, need to remove one of them
                    has_redundant_phase = true;
                    index_low = tmin1<tmin2?index1:index2;
                    index_high = tmax1>tmax2?index1:index2;
                    if (index_low == index_high)
                    {
                        index_rej = index_low == index1?index2:index1;
                        _removeRedundantPhase(phases,index_rej,index_low);
                    }
                    else
                    {
                        newkey = Phase::GetNPhases();
                        phases.insert(pair<int, Phase>(newkey,Phase(phase1,phase2)));
                        _removeRedundantPhase(phases,index_low,newkey);
                        _removeRedundantPhase(phases,index_high,newkey);
                    }
                    break;
                }
                else if (same_at_tmin || same_at_tmax)
                {
                    cout << "Not Implemented Yet" << endl;
                }
            }
            if (has_redundant_phase)
            {
                break;
            }            
        }
    }
}

/*
 * @description: Direct remove the phase, and reset the links
 * @param
 *      phases:
 *          The phase map
 *      index_removed:
 *          The key for the phase need to be removed
 *      index_phase:
 *          They key for the phase that is considered to be duplicated with the removed one
 * @return: 
 *      Directly modify `phases`
 */
void _removeRedundantPhase(MP &phases,int index_removed, int index_phase)
{
    SI low_trans_removed = phases.at(index_removed).GetLowTrans();
    SI high_trans_removed = phases.at(index_removed).GetHighTrans();
    SI::iterator iter;
    for (iter = low_trans_removed.begin(); iter != low_trans_removed.end(); iter++)
    {
        if (index_phase != (*iter))
        {
            phases.at(index_phase).addLowTrans(*iter);
            phases.at(*iter).eraseHighTrans(index_removed);
            phases.at(*iter).addHighTrans(index_phase);
        }
    }
    for (iter = high_trans_removed.begin(); iter != high_trans_removed.end(); iter++)
    {
        if (index_phase != (*iter))
        {
            phases.at(index_phase).addHighTrans(*iter);
            phases.at(*iter).eraseLowTrans(index_removed);
            phases.at(*iter).addLowTrans(index_phase);
        }
    }
    phases.erase(index_removed);
}
bool TransCompTc(TransCritical t1, TransCritical t2)
{
    return t1.Tcrit > t2.Tcrit;
}
struct minipara
{
    ScalarFunction f;
    Phase ph1;
    Phase ph2;
};
double Vcrit_wrap(double T, void *param)
{
    minipara *par = (minipara*)param;
    double v1 = (par->f)((par->ph1).valAt(T),T);
    double v2 = (par->f)((par->ph2).valAt(T),T);
    return v1-v2;
}

/*
 * @description: Find all possible critical temperature
 * @param
 *      phases:
 *          The phase map, from traceMultimin (and remove redundant phase)
 *      f:
 *          The potential
 * @return:
 *      A list of critical point
 *          Tcrit,Tnuc:
 *              The critical temperature
 *                  either Tcrit or Tnuc, the other one will be set to -1;
 *          low_vev, high_vev:
 *              The filed value at critical point for low and high temperature phase
 *          low_phase, high_phase:
 *              The key identifying the phase at low and high temperature
 *          trantype:
 *              1: first order
 *              2: second order
 */
VTC findCriticalTemperatures(MP phases, ScalarFunction f)
{
    VTC transitions;
    MP::iterator iter1,iter2;
    SI::iterator iterSI;
    int index1,index2;
    double Tcrit;
    for (iter1 = phases.begin(); iter1 != phases.end(); iter1++)
    {
        for (iter2 = phases.begin(); iter2 != phases.end(); iter2++)
        {
            if (iter1 == iter2)
            {
                continue;
            }
            index1 = iter1->first;
            index2 = iter2->first;
            Phase phase1 = phases.at(index1);
            Phase phase2 = phases.at(index2);
            double tmax = min(phase1.GetTmax(),phase2.GetTmax());
            double tmin = max(phase1.GetTmin(),phase2.GetTmin());
            if (tmin >= tmax)
            {
                // * No overlap, check whether it is second-order PT
                SI ph1SILow = phase1.GetLowTrans();
                if (ph1SILow.find(index2) != ph1SILow.end())
                {
                    transitions.push_back(secondOrderTrans(phase1,phase2,"Tcrit"));
                }
                continue;
            }

            // For 1-st order PT, just consider phase1->phase2
            if (f(phase1.valAt(tmin),tmin)-f(phase2.valAt(tmin),tmin) < 0)
            {
                continue;
            }
            if (f(phase1.valAt(tmax),tmax)-f(phase2.valAt(tmax),tmax) > 0)
            {
                continue;
            }
            function<double() > findTc = [&](){
                int status;
                int iter = 0;
                int max_iter = 200;
                const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
                gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
                minipara params = {f, phase1, phase2};
                gsl_function F;
                F.function = &Vcrit_wrap;
                F.params = &params;
                double t_low = tmin;
                double t_high = tmax;
                double t_crit;
                gsl_root_fsolver_set(s, &F, t_low, t_high);

                do
                {
                    iter++;
                    status = gsl_root_fsolver_iterate(s);
                    t_crit = gsl_root_fsolver_root(s);
                    t_low = gsl_root_fsolver_x_lower(s);
                    t_high = gsl_root_fsolver_x_upper(s);
                    status = gsl_root_test_interval (t_low,t_high,0.001,0.001);
                } while (status == GSL_CONTINUE && iter < max_iter);
                
                t_crit = gsl_root_fsolver_root(s);
                gsl_root_fsolver_free(s);
                return t_crit;
            };
            Tcrit = findTc();
            transitions.push_back({Tcrit,-1,phase2.valAt(Tcrit),phase1.valAt(Tcrit),index2,index1,1});
        }
    }
    sort(transitions.begin(),transitions.end(),TransCompTc);
    return transitions;
}

TransCritical secondOrderTrans(Phase high_phase, Phase low_phase, string str)
{
    TransCritical res;
    if (str=="Tcrit")
    {
        res.Tcrit = 0.5*(high_phase.GetTmin()+low_phase.GetTmax());
        res.Tnuc = -1;
    }
    else
    {
        res.Tnuc = 0.5*(high_phase.GetTmin()+low_phase.GetTmax());
        res.Tcrit = -1;
    }
    res.low_vev = low_phase.GetXatTmax();
    res.high_vev = high_phase.GetXatTmin();
    res.low_phase = low_phase.GetKey();
    res.high_phase = high_phase.GetKey();
    res.trantype = 2;
    return res;
}