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
    ScalarFunction V;
    dScalarFunction dV;
    HM d2V;
};




#endif //PATH_DEFORMATION_H