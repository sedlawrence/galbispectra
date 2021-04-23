
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
 #include "wigxjpf.h"


 struct galbispectra2
 {

   double *** first_order_sources;

   double *** second_order_sources_eq;

   double *** first_order_sources_integrand;

   double **** first_order_sources_integ;

   int flag;

   double ******* Dl;
   double ******k_integrand;
   double ***** Dl2;
   double ******* Dl3;

   double * tau_sampling_cls;
   double * tau_sampling_bessel;
   double ** tau_sampling_bessel2;
   double *** tau_sampling_selection_hires;
   double ** tau_sampling_cls_hires;
   double *** tau_sampling_lens_bi;
   double ** tau_rp;
   double * test_array;

   double ***** Cl;
   double ***** Cl3;

   ErrorMsg error_message;

   double * w_trapz_k; /* Corresponding to the array of trapezoidal weights w_trapz_k[index_k] for the k integration */

   double * w_trapz_r;
   double ** w_trapz_r2;

   double * w_trapz_alpha;

   double * w_trapz_lens;

   double * test_weights;

   int tau_size_selection;

   int tau_size_cls;

   int test_integer;

   int  tau_size_bessel;

   double * k_bessel;

   double * vecback; /* Pointer to array filled with background information */

   int r_size;

   int alpha_size;

   int alpha_log_size;

   int alpha_window_size;

   double *** r; /* Parameter to parameterise time grid */

   double ** r2;

   double * alpha; /* Parameter to parameterise time grid */

   int k_size_bessel; /* Number of points in the k_bessel grid */

   int type_size; /* Number of number count contributions */

   int source_size; /* Number of distinct terms in the first_order_sources array */



   int index_type_density; /* The index type corresponding to the rsd term in the galaxy number overdensity. */

   int index_type_rsd; /* The index type corresponding to the density term in the galaxy number overdensity.*/

   int index_type_delta;

   int index_type_g1;

   int index_type_g2;

   int index_type_g3;

   int index_type_lens;

   int index_type_g4;

   int index_type_g5;

   int index_type_quad_density_p;

   int index_type_quad_v;

   int index_type_quad_v_p; /* Index to denote the term that appears quadratically at second order -2*v'/a/a/H/H */

   int index_type_quad_v_pp; /* Index to denote the term that appears quadratically at second order v'' */

   int index_source_phi;

   int index_source_psi;

   int index_source_phi_plus_psi;

   int index_source_phi_plus_psi_prime;

   int index_source_phi_prime;

   int index_source_delta_cdm;

   int index_source_delta_b;
   int index_source_delta_m;

   int index_source_theta;  /* Velocity divergence = -k*k*v */

   int index_source_v; /* velocity */
   int index_type_vp; /* velocity */

   int index_type_d1; /* First Doppler term */

   int index_type_d2; /* Second Doppler term */

   /* Bispectrum types */
   /* Newtonian */
   int index_bisp_dens_mono; //3.19
   int index_bisp_dens_di; //3.20
   int index_bisp_dens_quad; //3.21
   int index_bisp_v_vpp; //3.25
   int index_bisp_vp_squared; //3.26
   int index_bisp_v_densp; // 3.27
   int index_bisp_vp_dens; // 3.28
   int index_bisp_so_rsd; // 3.46

   /* Terms that include lensing */
   int index_bisp_lens_dens; //3.32
   int index_bisp_vp_lens; //3.33
   int index_bisp_lens_squared; //3.34
   int index_bisp_Ddelta_Dpsi; //3.35
   int index_bisp_Dvp_Dpsi; // 3.36
   int index_bisp_Dlens_Dpsi; // 3.37
   int index_bisp_int_Dlens_DPsi1; //3.38
   int index_bisp_int_nabla2_DPsi1_DPsi1; // 3.39
   int index_bisp_so_lens; // 3.57

   int bisp_type_size;

   double ** tau_sampling_selection;

   double ******* redbi;

   double **** fo_dens_integ_hires_in_l;
   double ******* densdens_nDl1l2;
   int l_exact;
   int l_minus_two;
   int l_minus_one;
   int l_plus_one;
   int l_plus_two;

   int n_is_minus_one;
   int n_is_plus_one;
   int n_is_minus_two;
   int n_is_plus_two;
   int n_is_zero;
   int * n;

   int ** l_dual;

   double ****** integral_over_single_window;

   double ********* asym_redgalbispectrum;

   double ********** redgalbispectrum;

 };

 struct galbispectra2_workspace{


 };


 #endif

 /**************************************************************/
