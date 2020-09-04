/*
 * This is the examples from GSL Documents
 * https://www.gnu.org/software/gsl/doc/html/bspline.html#examples
*/

#include <cmath>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <iostream>
#include <fstream>

using namespace std;

/* Number of data points to fit */
#define N 500

/* Number of fit coefficients */
#define NCOEFFS 12

#define K 4 // Cubic order
/* NBREAK = NCOEFFS + 2 - K */
#define NBREAK (NCOEFFS + 2 - K)

int main(int argc, char const *argv[])
{
    const size_t n = N;
    const size_t ncoeffs = NCOEFFS;
    const size_t nbreak = NBREAK;
    size_t i, j;
    gsl_bspline_workspace *bw;
    gsl_vector *B;
    gsl_matrix *dB;
    double dy;
    gsl_rng *r;
    gsl_vector *c, *w;
    gsl_vector *x, *y;
    gsl_matrix *X, *cov;
    gsl_multifit_linear_workspace *mw;
    double chisq, Rsq, dof, tss;

    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);

    /* Allocate a cubic bspline workspace k=4 */
    bw = gsl_bspline_alloc(K,nbreak);
    B = gsl_vector_alloc(ncoeffs); // store the coefficient at each x

    x = gsl_vector_alloc(n); // Data x
    y = gsl_vector_alloc(n); // Data y
    w = gsl_vector_alloc(n); // Data 1/sigma^2;
    X = gsl_matrix_alloc(n,ncoeffs);
    c = gsl_vector_alloc(ncoeffs);
    cov = gsl_matrix_alloc(ncoeffs,ncoeffs);
    mw = gsl_multifit_linear_alloc(n,ncoeffs);

    // ! Generating the data that is needed to be fitted
    ofstream fdata("BSpline_Data.dat");
    fdata<<"xi\tyi"<<endl;
    for (i = 0; i < n; i++)
    {
        double sigma;
        double xi = (15.0/(N-1))*i; // * [0,15]
        double yi = cos(xi)*exp(-0.1*xi); 

        // ! Smearing the data, adding noise
        sigma = 0.1*yi;
        dy = gsl_ran_gaussian(r, sigma);
        yi += dy;

        // ! Store the data
        fdata<<xi<<"\t"<<yi<<endl;
        gsl_vector_set(x,i,xi);
        gsl_vector_set(y,i,yi);
        gsl_vector_set(w,i,1.0/(sigma*sigma));
    }
    fdata.close();

    // ! Using uniform breakpoints on [0,15]
    gsl_bspline_knots_uniform(0.0,15.0,bw);

    // ! Construct the fit matrix X
    for (i = 0; i < n; i++)
    {
        double xi = gsl_vector_get(x,i);

        // ! compute B_j(xi)
        gsl_bspline_eval(xi,B,bw);

        // ! fill X
        // cout<<xi;
        for (j = 0; j < ncoeffs; j++)
        {
            double Bj = gsl_vector_get(B,j);
            gsl_matrix_set(X,i,j,Bj);
            // cout<<"\t"<<Bj;
        }
        // cout<<endl;
    }

    // ! Do the fit
    // gsl_multifit_wlinear(X,w,y,c,cov,&chisq,mw);
    gsl_multifit_linear(X,y,c,cov,&chisq,mw);
    gsl_vector_set(c,0,gsl_vector_get(y,0));
    gsl_vector_set(c,ncoeffs-1,gsl_vector_get(y,n-1));
    dof = n - ncoeffs;
    // tss = gsl_stats_wtss(w->data,1,y->data,1,y->size);
    tss = gsl_stats_tss(y->data,1,y->size);
    Rsq = 1.0 - chisq/tss;
    cout<<"chisq/dof = "<<chisq/dof<<", Rsq = "<<Rsq<<endl;

    // ! Output the Smoothed Curve
    {
        double xi,yi,dyi,ddyi,yerr;
        ofstream fcurve("BSpline_curve.dat");
        fcurve<<"xi\tyi\tdyi\tddyi"<<endl;
        dB = gsl_matrix_alloc(ncoeffs,3);
        gsl_vector_view dBi;
        for ( xi = 0; xi < 15.0; xi += 0.1)
        {
            gsl_bspline_deriv_eval(xi,2,dB,bw);
            dBi=gsl_matrix_column(dB,0);
            gsl_multifit_linear_est(&dBi.vector,c,cov,&yi,&yerr);
            dBi=gsl_matrix_column(dB,1);
            gsl_multifit_linear_est(&dBi.vector,c,cov,&dyi,&yerr);
            dBi=gsl_matrix_column(dB,2);
            gsl_multifit_linear_est(&dBi.vector,c,cov,&ddyi,&yerr);
            fcurve<<xi<<"\t"<<yi<<"\t"<<dyi<<"\t"<<ddyi<<endl;
        }
        fcurve.close();
    }

    gsl_rng_free(r);
    gsl_bspline_free(bw);
    gsl_vector_free(B);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_vector_free(w);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(mw);

    return 0;
}

