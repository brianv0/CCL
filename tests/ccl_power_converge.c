#include "ccl.h"
#include <math.h>
#include <stdio.h>
#include "ccl_params.h"

int main(int argc, char * argv[])
{
  int status=0;
  double Omega_c = 0.25;
  double Omega_b = 0.05;
  double h = 0.7;
  double normp = 2.1e-9;
  double n_s = 0.96;

  ccl_configuration config = default_config;
  config.matter_power_spectrum_method=ccl_linear;

  ccl_parameters params = ccl_parameters_create_flat_lcdm(Omega_c, Omega_b, h, normp, n_s, &status);
  ccl_cosmology * cosmo = ccl_cosmology_create(params, config);

  FILE * f = fopen("./direct_class/ccl_linear_z0_pk.dat", "r");

  char str[1024];
  fgets(str, 1024, f);
  fgets(str, 1024, f);
  fgets(str, 1024, f);
  fgets(str, 1024, f);
  int i,nk=134;

  printf("# k [h/Mpc],P(k,z=0),P(k,z=3) [(Mpc/h)^3]\n");
  
  double k,p,p3=0;
  double a_at_z3=0.25;
  if(cosmo->config.matter_power_spectrum_method==ccl_linear){
    for (i=0;i<nk;i++){
      fscanf(f,"%le %*le\n",&k);
      p = ccl_linear_matter_power(cosmo, k*h,1.0, &status);
      p3 = ccl_linear_matter_power(cosmo,k*h, a_at_z3,&status);
      printf("%le %le %le \n", k, p*h*h*h,p3*h*h*h);
    }
  } else {
    if(cosmo->config.matter_power_spectrum_method==ccl_halofit){
      for (i=0;i<nk;i++){
	fscanf(f,"%le %*le\n",&k);
	p = ccl_nonlin_matter_power(cosmo, k*h,1.0,&status);
	p3 = ccl_nonlin_matter_power(cosmo,k*h, a_at_z3,&status);
	printf("%le %le %le \n", k, p*h*h*h,p3*h*h*h);
      }
    } else {
      printf("ccl_power_converge.c: Unknown power spectrum method.\n");
      return NAN;
    }
  }
  
  fclose(f);
  ccl_cosmology_free(cosmo);

  return 0;

}
