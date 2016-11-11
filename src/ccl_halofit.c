
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ccl_core.h"
#include <gsl/gsl_integration.h>

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

double ccl_sigmaR_gauss(ccl_cosmology* cosmo, double R){
	
  SigmaR_pars par;
  par.cosmo=cosmo;
  par.R=R;

  gsl_integration_cquad_workspace *workspace=gsl_integration_cquad_workspace_alloc(1000);
  gsl_function F;
  F.function=&sigmaR_gauss_integrand;
  F.params=&par;

  double sigma_R;
  gsl_integration_cquad(&F,log10(K_MIN_INT),log10(K_MAX_INT),0.0,1E-5,workspace,&sigma_R,NULL,NULL);
  //TODO: log10 could be taken already in the macros.
  //TODO: 1E-5 should be a macro
  //TODO: we should check for integration success
  gsl_integration_cquad_workspace_free(workspace);

  return sqrt(sigma_R*M_LN10/(2*M_PI*M_PI));
}
