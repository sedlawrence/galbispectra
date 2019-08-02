
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



   double ***** Cl;

   double * tau_sampling_cls;

   double * tau_sampling_bessel;

   double ***** Cl_final;

   ErrorMsg error_message;

   double * w_trapz_k; /* Corresponding to the array of trapezoidal weights w_trapz_k[index_k] for the k integration */

   double * w_trapz_lens;

   double * test_weights;

   int tau_size_selection;

   int tau_size_bessel;

   double * k_bessel;



   int k_size_bessel; /* Number of points in the k_bessel grid */

   int type_size; /* Number of different terms in the first_order_sources array */

   int index_type_delta_cdm; /* The index type corresponding to the rsd term in the galaxy number overdensity. */

   int index_type_rsd; /* The index type corresponding to the rsd term in the galaxy number overdensity.*/

   int index_type_v; /* velocity */

   int index_type_theta; /* velocity divergence, theta = -k*k*v */

   int index_type_d1; /* Index type corresponding to the D1 (Doppler contribution) */

   int index_type_d2; /* Index type corresponding to the D2 (Doppler contribution) */

   int index_type_rsd_gr; /* Anisotropic stress is neglected such that Phi = Psi */

   int index_type_lens;

   int index_type_phi_plus_psi;

   int index_type_phi_prime;

   int index_type_phi;

   double ** tau_sampling_selection;

   double ****** integral_over_single_window;

   double ********* asym_redgalbispectrum;

   double ********** redgalbispectrum;

 };

 struct galbispectra2_workspace{


 };


 #endif

 /**************************************************************/
