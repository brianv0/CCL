
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ccl_core.h"
#include "ccl_power.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>

#define R_LOW 0.1
#define R_HIGH 10

//Params for sigma(R) integrand
typedef struct {
  ccl_cosmology *cosmo;
  double R;
} SigmaR_pars;

typedef struct {
	double a, b, c, gamma, alpha, beta, mu, nu, f1, f2, f3;
} halofit_param;


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

static double log_sigma(double R, void *params)
{
    ccl_cosmology* cosmo = (ccl_cosmology*)params;
    return log(ccl_sigmaR_gauss(cosmo, R));
}

static double sigmaR_deriv(double R, void *params)
{ // d ln(sigma) / d ln(r)
	gsl_function F;
	double result, abserr;

	F.function = &log_sigma;
	F.params = params;

	gsl_deriv_central (&F, R, 1e-6, &result, &abserr);
	result *= R; // d ln R = 1/R * dR
	
//	printf("R = %f, derivative d ln(sigma) / d ln(r) = %.10f (+/- %.10f)\n", R, result, abserr);
	return result;
}

double ccl_n_eff(ccl_cosmology* cosmo)
{ // returns effective spectral index
	return -sigmaR_deriv(ccl_k_sigma_m(cosmo), cosmo)-3;
}

static double sigmaR_gauss_2nd_deriv(double R, void *params)
{ // d^2 ln(sigma) / d ln(r)^2
	gsl_function F;
	double result, abserr;

	F.function = &sigmaR_deriv;
	F.params = params;
	
	gsl_deriv_central (&F, R, 1e-6, &result, &abserr);
	result *= R; // d ln R = 1/R * dR

//	printf("R = %f, derivative d^2 ln(sigma) / d ln(r)^2 = %.10f (+/- %.10f)\n", R, result, abserr);
	return result;
}

double ccl_curvature(ccl_cosmology* cosmo)
{ // returns curvature
	return -sigmaR_gauss_2nd_deriv(ccl_k_sigma_m(cosmo), cosmo);
}

void ccl_set_Takahashi_fit(ccl_cosmology* cosmo, halofit_param* param, double n_eff, double C, double a)
{
	double w = cosmo->params.w0 + (1-a)*cosmo->params.wa;
	
	param->a = 1.5222 + 2.8553*n_eff + 2.3706*pow(n_eff, 2.) + 0.9903*pow(n_eff, 3.) + 0.2250*pow(n_eff, 4.) \
				- 0.6038*C +0.1749*cosmo->params.Omega_l*(1.+w);
	param->a = pow(10., param->a);
	
	param->b = -0.5642 + 0.5864*n_eff + 0.5716*pow(n_eff, 2.) + 1.5474*C + 0.2279*cosmo->params.Omega_l*(1+w);
	param->b = pow(10., param->b);
	
	param->c = 0.3698 + 2.0404*n_eff + 0.8161*pow(n_eff, 2.) + 0.5869*C;
	param->c = pow(10., param->c);
	
	param->gamma = 0.1971 - 0.0843*n_eff + 0.8460*C;
	
	param->alpha = abs(6.0835 + 1.3373*n_eff - 0.1959*pow(n_eff, 2.) - 5.5274*C);
	
	param->beta = 2.0379 - 0.7354*n_eff + 0.3157*pow(n_eff, 2.) + 1.2490*pow(n_eff, 3.) + 0.3980*pow(n_eff, 4.) - 0.1682*C;
	
	param->mu = 0;
	
	param->nu = 5.2105 + 3.6902*n_eff;
	param->nu = pow(10., param->nu);
	
	param->f1 = pow(cosmo->params.Omega_m, -0.0307);
	
	param->f2 = pow(cosmo->params.Omega_m, -0.0585);
	
	param->f3 = pow(cosmo->params.Omega_m, 0.0743);
}