#pragma once

// sigmaR defined with gaussian window function
double ccl_sigmaR_gauss(ccl_cosmology* cosmo, double R);

// the nonlinear scale (k_sigma)^-1
double ccl_k_sigma_m(ccl_cosmology* cosmo);
