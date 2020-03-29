/*
 * This is the examples from GSL Documents
 * https://www.gnu.org/software/gsl/doc/html/bspline.html#examples
*/

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "GSL_Wraper.h"
#include <iostream>
#include <fstream>

using namespace std;

/* Number of data points to fit */
#define N 200


int main(int argc, char const *argv[])
{
    const int n = N;

    // ! Generating the data that is needed to be fitted
    // double *x = new double[n];
    // double *y = new double[n];
    VD x(n);
    VD y(n);
    ofstream fdata("Spline_Data.dat");
    fdata<<"xi\tyi"<<endl;
    for (int i = 0; i < n; i++)
    {
        double sigma;
        double xi = (15.0/(N-1))*i; // * [0,15]
        double yi = cos(xi)*exp(-0.1*xi); 

        // ! Store the data
        fdata<<xi<<"\t"<<yi<<endl;
        x[i] = xi;
        y[i] = yi;
    }
    fdata.close();

    
    GSL_Spline_Inter spline;
    spline.SetData(&y,&x);
    // gsl_interp_accel *acc = gsl_interp_accel_alloc();
    // gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,n);
    // gsl_spline_init(spline,x,y,n);
    // // ! Output the Smoothed Curve
    // {
    //     double xi,yi,dyi,ddyi,yerr;
    //     ofstream fcurve("Spline_inter_curve.dat");
    //     fcurve<<"xi\tyi\tdyi\tddyi"<<endl;
    //     for ( xi = -3; xi < 18.0; xi += 0.1)
    //     {
    //         yi = gsl_spline_eval(spline,xi,acc);
    //         dyi = gsl_spline_eval_deriv(spline,xi,acc);
    //         ddyi = gsl_spline_eval_deriv2(spline,xi,acc);
    //         fcurve<<xi<<"\t"<<yi<<"\t"<<dyi<<"\t"<<ddyi<<endl;
    //     }
    //     fcurve.close();
    // }
    {
        double xi,yi,dyi,ddyi,yerr;
        ofstream fcurve("Spline_inter_curve.dat");
        fcurve<<"xi\tyi\tdyi\tddyi"<<endl;
        for ( xi = -3; xi < 18.0; xi += 0.1)
        {
            yi = spline.valAt(xi,0);
            dyi = spline.valAt(xi,1);
            ddyi = spline.valAt(xi,2);
            fcurve<<xi<<"\t"<<yi<<"\t"<<dyi<<"\t"<<ddyi<<endl;
        }
        fcurve.close();
    }



    return 0;
}

