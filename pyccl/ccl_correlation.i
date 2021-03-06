%module ccl_correlation

%{
#define SWIG_FILE_WITH_INIT
#include "../include/ccl_correlation.h"
%}

// Automatically document arguments and output types of all functions
%feature("autodoc", "1");

// Strip the ccl_ prefix from function names
%rename("%(strip:[ccl_])s") "";

%include "../include/ccl_correlation.h"

// Enable vectorised arguments for arrays
%apply (int DIM1,double* IN_ARRAY1) {
                                     (int nlarr,double* larr),
                                     (int nclarr,double* clarr),
                                     (int nt,double *theta)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* output, int nout)};

%inline %{

void correlation_vec(ccl_cosmology *cosmo,
		     int nlarr,double *larr,
		     int nclarr,double *clarr,
		     int nt,double *theta,
		     int corr_type,int method,
		     double *output,int nout,
		     int *status)
{
  assert(nlarr==nclarr);
  assert(nt==nout);

  ccl_correlation(cosmo,nlarr,larr,clarr,nt,theta,output,corr_type,0,NULL,method,status);
}

%}
