#pragma once

// sigmaR defined with gaussian window function
double ccl_sigmaR_gauss(ccl_cosmology* cosmo, double R);

// the nonlinear scale (k_sigma)^-1
double ccl_k_sigma_m(ccl_cosmology* cosmo);

// effective spectral index
double ccl_n_eff(ccl_cosmology* cosmo);

// curvature C
double ccl_curvature(ccl_cosmology* cosmo);