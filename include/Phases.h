/*
 * @Description  : Phase Class store the information of a phase in a temperature period
 * @Author       : Yongcheng Wu
 * @Date         : 2019-12-20 17:28:39
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-16 16:03:43
 */
#ifndef PHASES_H
#define PHASES_H
#include <string>
#include "VTypes.h"
#include <gsl/gsl_spline.h>

typedef std::vector<gsl_spline*> VSP;

class Phase
{
private:
    static int NPhases;
    int key; // Phase ID
    int DimX; // Dimension of X
    int DimData; // Dimension of T
    VVD X; // Store the Field information at every temperature
    VD T; // Temperature
    VVD dXdT; // The derivative of X respect to Temperature
    SI low_trans;
    SI high_trans;

    gsl_interp_accel *acc;
    VSP inters;

public:
    Phase();
    Phase(int _key, VVD _X, VD _T, VVD _dXdT);
    Phase(const Phase &Ph);
    Phase(const Phase &ph1, const Phase &ph2); //Merge two phases // ! When using this constructor, please make sure one of the two phases are redundant. We don't do such check
    ~Phase();

    void SetPhase(int _key, VVD _X, VD _T, VVD _dXdT);
    void SetUpInterpolation();
    void FreeInterpolation();
    VD valAt(double _T);
    int GetKey() const {return key;}
    double GetTmin() const {return T.front();}
    double GetTmax() const {return T.back();}
    VD GetXatTmin() const {VD Xmin; for(int i = 0; i < DimX; i++ ) {Xmin.push_back(X[i].front());} return Xmin;}
    VD GetXatTmax() const {VD Xmax; for(int i = 0; i < DimX; i++ ) {Xmax.push_back(X[i].back());} return Xmax;}
    void addLowTrans(int _key) {low_trans.insert(_key);}
    void eraseLowTrans(int _key) {low_trans.erase(_key);}
    void addHighTrans(int _key) {high_trans.insert(_key);}
    void eraseHighTrans(int _key) {high_trans.erase(_key);}
    void addLinkFrom(Phase *other_phase);
    std::string repr();

    VVD GetX() const {return X;}
    VVD GetdXdT() const {return dXdT;}
    VD GetT() const {return T;}
    SI GetLowTrans() const {return low_trans;}
    SI GetHighTrans() const {return high_trans;}
    static int GetNPhases() {return NPhases;}
};

#endif