#ifndef PATH_DEFORMATION_H
#define PATH_DEFORMATION_H

#include "VTypes.h"
#include "GSL_Wraper.h"

class SplinePath
{
private:
    int _NDim;
    double _L;

    GSL_Spline_Inter _V_Inter;

    GSL_Multi_Spline_Inter _path_Inter;


public:
    SplinePath(VVD pts, ScalarFunction V_, double T,int V_spline_samples = 100, bool extend_to_minima=false, bool re_eval_distances=true);
    ~SplinePath(){};

    VD pts_at_dist(double x, int deri=0);
    VVD pts_at_dist(VD X);
    double GetDistance(){return _L;};
    ScalarFunction V;
    dScalarFunction dV;
    HM d2V;
};
enum Deformation_Status {
    REACHMAXITER = -2,
    DIVERGENT = -1,
    FAILCALLBACK = 0,
    DF_CONVERGED = 1
};
class Deformation_Spline
{
private:
    double _L; // Total length of the path;
    VD _t; // List from 0 to 1 which parameterized the path.
    GSL_BSpline_Fit _fitter;

    VVD _phi;

    dScalarFunction _dV;
    VD _v2; // The `velocity square along the path`

    bool _fix_start;
    bool _fix_end;
    bool _save_all_steps;
    int _num_steps;
    std::vector<VVD> _F_list;
    std::vector<VVD> _phi_list;
    VVD _phi_last;
    VVD _F_last;

public:
    Deformation_Spline(VVD phi, VD dphidr, dScalarFunction dV, int ncoeffs = 10, int K = 4, double v2min = 0.0, bool fix_start = false, bool fix_end = false, bool save_all_steps = false);
    ~Deformation_Spline();

    std::tuple<VVD,VVD> forces();
    std::tuple<double,bool,double> step(double last_step, double maxstep=0.1, double minstep=1e-4, double reverseCheck=0.15, double stepIncrease=1.5, double stepDecrease=5.0, bool checkAfterFit=true);
    Deformation_Status deformPath(double start_step=2e-3,double fRatioConv=1e-2,double converge_0 = 5.0, double fRatioIncrease=5.0, int maxiter=500, std::function<bool(Deformation_Spline*)> callback=[](Deformation_Spline* cur){return true;});

    VVD GetPhi(){return _phi;};
    int GetSteps(){return _num_steps;};
};



#endif //PATH_DEFORMATION_H