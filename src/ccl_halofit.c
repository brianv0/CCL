
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ccl_core.h"
#include "ccl_power.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

#define R_LOW 0.1
#define R_HIGH 10

//Params for sigma(R) integrand
typedef struct {
  ccl_cosmology *cosmo;
  double R;
} SigmaR_pars;

static double sigmaR_gauss_integrand(double lk,void *params)
{
  SigmaR_pars *par=(SigmaR_pars *)params;
  double k=pow(10.,lk);
  double pk=ccl_linear_matter_power(par->cosmo,1.,k);
  double kR=k*par->R;
  double w;
  // gaussian window function
  w = exp(-kR*kR);
  
  return pk*k*k*k*w;
}

double ccl_sigmaR_gauss(ccl_cosmology* cosmo, double R)
{ // returns sigma_r with gaussian window function

  SigmaR_pars par;
  par.cosmo=cosmo;
  par.R=R;

  gsl_integration_cquad_workspace *workspace=gsl_integration_cquad_workspace_alloc(1000);
  gsl_function F;
  F.function=&sigmaR_gauss_integrand;
 // F.function=&sigmaR_integrand;
  F.params=&par;

  double sigma_R;
  gsl_integration_cquad(&F,log10(K_MIN_INT),log10(K_MAX_INT),0.0,1E-9,workspace,&sigma_R,NULL,NULL);
  //TODO: log10 could be taken already in the macros.
  //TODO: 1E-5 should be a macro
  //TODO: we should check for integration success
  gsl_integration_cquad_workspace_free(workspace);

  return sqrt(sigma_R*M_LN10/(2*M_PI*M_PI));
}


static double sigmaR_gauss_root_finding(double R, void *params)
{
    ccl_cosmology* cosmo = (ccl_cosmology*)params;
    return ccl_sigmaR_gauss(cosmo, R) - 1.;
}



double ccl_k_sigma_m(ccl_cosmology* cosmo){
// returns the nonlinear scale (k_sigma)^-1

    int status;
    int iter = 0;
    double R, R_low, R_high;

    gsl_function F;
    F.function = &sigmaR_gauss_root_finding;
    F.params = cosmo;

    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, 1/K_MAX_INT, 1/K_MIN_INT);

    do
    {
	iter++;
        status = gsl_root_fsolver_iterate (s);
        R_low = gsl_root_fsolver_x_lower (s);
        R_high = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (R_low, R_high, 0, 1E-5);
    }
    while (status == GSL_CONTINUE);
    printf("# Root found in %i iterations.\n", iter);
    R = gsl_root_fsolver_root (s);
    gsl_root_fsolver_free (s);

    return R;
}
