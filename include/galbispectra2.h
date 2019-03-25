
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

   double *** Cl;

   double * tau_sampling_cls;

   double *** Cl_final;

   ErrorMsg error_message;

   double * w_trapz_k; /* Corresponding to the array of trapezoidal weights w_trapz_k[index_k] for the k integration */

   int tau_size_selection;

   double ** tau_sampling_selection;

 };

 struct galbispectra2_workspace{


 };

/*
 int galbispectra2_init (
      struct precision * ppr,
      struct precision2 * ppr2,
      struct background * pba,
      struct thermo * pth,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct galbispectra2 * pgb2,
      struct bessels * pbs,
      struct transfers * ptr
    );
    */

 #endif

 /**************************************************************/
