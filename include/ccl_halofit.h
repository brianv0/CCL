#pragma once

#include "ccl_core.h"

typedef struct {
	double a, b, c, gamma, alpha, beta, mu, nu, f1, f2, f3;
} halofit_param;

// ccl_cosmology_halofit :: to be incorporated into ccl_cosmology
typedef struct
{
	halofit_param param;
	double k_sigma_m, n_eff, C;
	bool computed_halo_fit;
	
} ccl_cosmology_halofit;

ccl_cosmology_halofit ccl_new_ccl_cosmology_halofit();
void ccl_cosmology_init_halofit(ccl_cosmology* cosmo, double a, ccl_cosmology_halofit* cosmo_halo_fit);
double ccl_nonlin_matter_power_halofit(ccl_cosmology* cosmo, double a, double k, ccl_cosmology_halofit* cosmo_halo_fit);