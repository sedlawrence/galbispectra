
/**
 * Workspace that contains the intermediate results for the integration of the galaxy
 * bispectrum.
 *
 */



 /** @file galbispectra2.h Documented header file for the galbispectra2.c module */

 #ifndef __GALBISPECTRA2__
 #define __GALBISPECTRA2__

 #include "common.h"
 #include "common2.h"
 #include "input.h"
 #include "perturbations.h"
 #include "perturbations2_macros.h"


 struct galbispectra2
 {

   double *** first_order_sources;

   double *** first_order_sources_integrand;

   double **** first_order_sources_integ;



   double ***** Dl;
   double ***** Dl2;

   double * tau_sampling_cls;

   double * tau_sampling_bessel;
   double ** tau_rp;
   double * test_array;

   double ***** Cl;

   ErrorMsg error_message;

   double * w_trapz_k; /* Corresponding to the array of trapezoidal weights w_trapz_k[index_k] for the k integration */

   double * w_trapz_r;

   double * w_trapz_alpha;

   double * w_trapz_lens;

   double * test_weights;

   int tau_size_selection;

   int tau_size_bessel;

   double * k_bessel;

   double * vecback; /* Pointer to array filled with background information */

   int r_size;

   int alpha_size;

   double * r; /* Parameter to parameterise time grid */

   double * alpha; /* Parameter to parameterise time grid */

   int k_size_bessel; /* Number of points in the k_bessel grid */

   int type_size; /* Number of number count contributions */

   int source_size; /* Number of distinct terms in the first_order_sources array */



   int index_type_density; /* The index type corresponding to the rsd term in the galaxy number overdensity. */

   int index_type_rsd; /* The index type corresponding to the density term in the galaxy number overdensity.*/

   int index_type_g1;

   int index_type_g2;

   int index_type_g3;

   int index_type_lens;

   int index_type_g4;

   int index_type_g5;

   int index_source_phi;

   int index_source_psi;

   int index_source_phi_plus_psi;

   int index_source_phi_plus_psi_prime;

   int index_source_phi_prime;

   int index_source_delta_cdm;

   int index_source_theta;  /* Velocity divergence = -k*k*v */

   int index_source_v; /* velocity */

   int index_type_d1; /* First Doppler term */

   int index_type_d2; /* Second Doppler term */

   double ** tau_sampling_selection;

   double ****** integral_over_single_window;

   double ********* asym_redgalbispectrum;

   double ********** redgalbispectrum;

 };

 struct galbispectra2_workspace{


 };


 #endif

 /**************************************************************/
