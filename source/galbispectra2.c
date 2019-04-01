/* Module to compute galaxy bispectra */


#include "galbispectra2.h"
#include "perturbations2.h"




//herehere
/* Define integral over k to give angular power spectrum without window functions */
int integral(
  struct background * pba, /* input pointer to background */
  struct primordial * ppm,
  struct bessels * pbs, /* input, pointer to bessels structure */
  struct galbispectra2 * pgb2, /* input pointer to galaxy bispectra structure */
  struct perturbs * ppt,
  struct transfers * ptr,
  int index_type, /* index corresponding to type (e.g. delta_m) */
  int index_tau_first,  /* index corresponding to ppt->tau_sampling_cls[index_tau] */
  int index_tau_second, /* index corresponding to ppt->tau_sampling_cls[index_tau] */
  int index_l,
  double * result) /**< order of the Bessel functions from pbs->l[index_l] */
  {
  int index_md;
  double tmp = 0.;

  double j1;
  double j2;
  double Pk;


  //for(int index_k = 0; index_k < ppt->k_size[ppt->index_md_scalars]; index_k++){
  //for(int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++){
    //double x1 = ppt->k[ppt->index_md_scalars][index_k]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first]);
    //double x2 = ppt->k[ppt->index_md_scalars][index_k]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second]);
    //double x1 = pgb2->k_bessel[index_k_bessel];/*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first]);*/
    //double x2 = pgb2->k_bessel[index_k_bessel];/*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second]);*/
    //if(x1>pbs->x_max || x2 > pbs->x_max){/*printf("Skipping x=%.10e or %.10e \n",x1,x2);*/continue;}
    //if(x1>pbs->x_max){/*printf("Skipping x=%.10e or %.10e \n",x1,x2);*/continue;}
    //if(index_l > pbs->l_size){printf("Skipping l=%.10e \n", ptr->l[index_l]);continue;}
    //class_call(bessel_at_x(pbs,	x1, index_l, &j1), pbs->error_message, pgb2->error_message);
    //class_call(bessel_at_x(pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);
    //class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, ppt->k[ppt->index_md_scalars][index_k], &Pk), ppm->error_message, pgb2->error_message);
    //printf("j_%i (%e) = %e\n", ptr->l[index_l], x1, j1);
    //printf("P(%g)= %g\n", ppt->k[ppt->index_md_scalars][index_k], Pk);

/* Using area of trapezoid, we can sum a number of areas of trapezoids to approximate the integral */
    /*tmp+=pow(ppt->k[ppt->index_md_scalars][index_k],-1.) * 4. * _PI_ * pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau_first][index_k]
      * pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau_second][index_k] * Pk
        * j1 * j2  * pgb2->w_trapz_k[index_k];*/



    for(int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++){

      double x1 = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first]);
      double x2 = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second]);

      if(x1>pbs->x_max || x2 > pbs->x_max){
      /*printf("ALERT! x1= %g x2 = %g \n",x1,x2);*/continue;}

      class_call(bessel_at_x(pbs,	x1, index_l, &j1), pbs->error_message, pgb2->error_message);
      class_call(bessel_at_x(pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);
      class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, pgb2->k_bessel[index_k_bessel], &Pk), ppm->error_message, pgb2->error_message);


      /* Using area of trapezoid, we can sum a number of areas of trapezoids to approximate the integral */
      tmp += pow(pgb2->k_bessel[index_k_bessel],-1.) * 4. * _PI_ * pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau_first][index_k_bessel]
        * pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau_second][index_k_bessel] * Pk
          * j1 * j2  * pgb2->w_trapz_k[index_k_bessel];


      //printf("tmp in k loop outside integral function = %g\n", tmp);
      //printf("w_trapz_k = %g \n", pgb2->w_trapz_k[index_k_bessel]);

    }




      //printf("pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau_second][index_k] = %g\n", pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau_second][index_k] );
      //printf("ppt->k[ppt->index_md_scalars][index_k] = %g\n",ppt->k[ppt->index_md_scalars][index_k] );

    //printf("tmp in k loop inside integral function = %.10e\n", tmp );


  //printf("TEMP = %.10e\n", tmp);
  *result = tmp;
  return _SUCCESS_;
} //End of integral()



/***************************************************************************
====================    Define Search Functions ============================
****************************************************************************/

/* index_tau() - function to assign indices to tau values

* @param tau1                     Input: conformal time
* @param tau2                     Input: conformal time
* @param index_tau1               Output: index in the array tau_sampling_cls
* @param index_tau1               Output: index in the array tau_sampling_cls
* @param perturbs                 Input: perturbs structure
*/

int index_of_tau_sampling_cls(double tau1,
                 int * index_tau1,
                 struct galbispectra2 * pgb2){

                 * index_tau1 = 0;
                 double tau;

                  /* Scan through the tau_sampling_cls grid and assign the index which gives the best estimate of tau1 and tau2 */

                    for (int index = 0; index < pgb2->tau_size_selection ; index++) {
                      tau = pgb2->tau_sampling_cls[index];

                      if (tau > tau1 && *index_tau1 == 0) {
                        *index_tau1=index;
                      }



                      if (*index_tau1 != 0){
                        break;
                      }



                    }
                  }

int index_of_tau_sampling_quadsources(double tau1,
                 int * index_tau1,
                 struct perturbs * ppt){

                 * index_tau1 = 0;

                 double tau;

                  /* Scan through the tau_sampling_quadsources grid and assign the index which gives the best estimate of tau1 and tau2 */

                    for (int index = 0; index < ppt->tau_size_quadsources ; index++) {
                      tau = ppt->tau_sampling_quadsources[index];

                      if (tau > tau1 && *index_tau1 == 0) {
                        *index_tau1=index;
                      }

                      if (*index_tau1 != 0){
                        break;
                      }



                    }
                  }


/* index_of_k() - function to assign indices of ppt->k to a given k value

* @param k                    Input: k-value
* @param index_tau1           Output: index in the array ppt->k
* @param perturbs             Input: perturbs structure
*/



int index_of_k(double k,
               int * index_k,
               struct perturbs * ppt){

               * index_k = 0;
               double k_ppt;

                  /* Scan through the k grid pgb2->k_bessel[index] and assign the index which gives the closest estimate of k */

                    for (int index = 0; index < ppt->k_size[ppt->index_md_scalars]; index++) {
                      k_ppt = ppt->k[ppt->index_md_scalars][index];
 /* NOTE: the sign < >*/
                      if (k < k_ppt && *index_k == 0) {
                        *index_k=index;
                      }

                      if (*index_k != 0){
                        break;
                      }



                    }
                  }


/* index_of_k_bessel() - function to assign indices of ppt->k to a given k value

* @param k                    Input: k-value
* @param index_tau1           Output: index in the array pgb2->k_bessel
* @param perturbs             Input: perturbs structure
*/

int index_of_k_bessel(double k,
               int * index_k_bessel,
               struct galbispectra2 * pgb2){

               * index_k_bessel = 0;
               double k_bessel;

                  /* Scan through the k_bessel grid pgb2->k_bessel[index] and assign the index which gives the closest estimate of k */

                    for (int index = 0; index < pgb2->k_size_bessel; index++) {
                      k_bessel = pgb2->k_bessel[index];
 /* NOTE: the sign < >*/
                      if (k < k_bessel && *index_k_bessel == 0) {
                        *index_k_bessel=index;
                      }

                      if (*index_k_bessel != 0){
                        break;
                      }



                    }
                  }



/* Initialise the module, this function is called by the song.c main file. */

int galbispectra2_init (
     struct precision * ppr,
     struct precision2 * ppr2,
     struct background * pba,
     struct primordial * ppm,
     struct thermo * pth,
     struct perturbs * ppt,
     struct perturbs2 * ppt2,
     struct galbispectra2 * pgb2,
     struct bessels * pbs,
     struct transfers * ptr
     )
{
    /* Define local variables */

  int index_l;
  int index_k;
  int index_tau;
  int index_tau_first;
  int index_tau_second;
  int index_type;
  int index_md;
  int bin_first;
  int bin_second;
  double * j;
  double * pvecback;
  double * w_trapz_tau;
  double tau0 = pba->conformal_age;
  double ** w_trapz;
  double ** tau0_minus_tau;
  int i;
  int bin;
  pgb2->tau_size_selection = 50;
  pgb2->k_size_bessel = 3000;
  printf("Starting galaxy bispectra module...\n");

      /* Manually set selection function paramters. NOTE: this should be accounted for in the .ini file */

      //ppt->selection_mean[0] = 1.0;
      //ppt->selection_mean[1] = 1.5;
    //  ppt->selection_mean[2] = 2.0;
    //  ppt->selection_num = 3;

  /*
  if (ppt2->has_perturbations2 == _FALSE_) {
    if (ppt2->perturbations2_verbose > 0)
      printf("No second-order sources requested. Second-order perturbations module skipped.\n");
    return _SUCCESS_;
  }
  */
  class_alloc(pvecback,pba->bg_size*sizeof(double), pgb2->error_message);
  if (ppt2->perturbations2_verbose > 0)
    printf("Computing second-order perturbations\n");




  // ====================================================================================
  // =                              Indices and samplings                               =
  // ====================================================================================

  /* Determine which sources need to be computed and their k-sampling */

/*  class_call( perturb2_indices_of_perturbs(
                ppr,
                ppr2,
                pba,
                pth,
                ppt,
                ppt2),
  ppt2->error_message,
  ppt2->error_message); */


  /* Determine the time sampling for the sources */

/*  class_call (perturb2_timesampling_for_sources (
                ppr,
                ppr2,
                pba,
                pth,
                ppt,
                ppt2),
    ppt2->error_message,
    ppt2->error_message); */


  // ====================================================================================
  // =                            Solve first-order system                              =
  // ====================================================================================

  /* Run the first-order perturbations module in order to:
    - Compute and store in ppt->quadsources the perturbations needed to solve the
      2nd-order system.
    - Compute and store in ppt->sources the line-of-sight sources needed to compute
      the first-order transfer functions (not the C_l's).
    The k-sampling for the first-order perturbations (ppt->k) has been already determined
    in perturb2_get_k_lists() and matches the one for the second-order sources (ppt2->k).
    Similarly, their time sampling (ppt->tau_sampling_quadsources) has been computed
    in perturb2_timesampling_for_sources().  */

  /*class_call (perturb_init (
                ppr,
                pba,
                pth,
                ppt),
    ppt->error_message, ppt2->error_message);*/


  /* Stop here if the user asked to compute only the first-order perturbations */
  if (ppt2->stop_at_perturbations1 == _TRUE_) {

    ppt->has_perturbations = _FALSE_;
    ppt->has_cls = _FALSE_;
    ppt->has_cmb_bispectra = _FALSE_;
    ppt2->has_cmb_spectra = _FALSE_;
    ppt2->has_cmb_bispectra = _FALSE_;

    if (ppt2->perturbations2_verbose > 0)
      printf(" -> Exiting after computation of first-order perturbations\n");

    return _SUCCESS_;
  }

  /* Print some info to screen */
  if (ppt2->perturbations2_verbose > 0) {
    printf(" -> computing %s2nd-order sources ",
      ppt2->rescale_cmb_sources==_TRUE_ ? "RESCALED " : "");

    if (ppt->gauge == newtonian)
      printf("in Newtonian gauge ");
    if (ppt->gauge == synchronous)
      printf("in synchronous gauge ");

    printf ("for m=");
    for (int index_m=0; index_m < (ppr2->m_size-1); ++index_m)
      printf("%d,", ppr2->m[index_m]);
    printf("%d\n", ppr2->m[ppr2->m_size-1]);
  }


  /* Apart from ppt2->sources, all the arrays needed by the subsequent modules have been filled.
  If the user requested to load the line of sight sources from disk, we can stop the execution of
  this module now without regrets. */

  if (ppr2->load_sources_from_disk == _TRUE_) {

    if (ppt2->perturbations2_verbose > 0)
      printf(" -> leaving perturbs2 module; line-of-sight sources will be read from disk\n");

    /* Uncomment to produce the sources output files again */
    // if ((ppt2->k_out_size > 0) || (ppt2->tau_out_size > 0))
    //   class_call_parallel (perturb2_output (
    //                          ppr,
    //                          ppr2,
    //                          pba,
    //                          ppt,
    //                          ppt2),
    //      ppt2->error_message,
    //      ppt2->error_message);

    return _SUCCESS_;

  }


  /* Do not evaluate the subsequent modules if ppt2->stop_at_perturbations2 == _TRUE_ */
  if (ppt2->stop_at_perturbations2 == _TRUE_) {
    ppt->has_perturbations = _FALSE_;
    ppt->has_cls = _FALSE_;
    ppt->has_cmb_bispectra = _FALSE_;
    ppt2->has_cmb_spectra = _FALSE_;
    ppt2->has_cmb_bispectra = _FALSE_;
  }





  /* Define an array of values of first order transfer functions:
              pgb2->first_order_sources[index_type][index_tau][index_k_bessel] */

  class_alloc(pgb2->first_order_sources, ppt->qs_size[ppt->index_md_scalars] * sizeof(double **), ppt->error_message);
    for (int index_type = 0; index_type < ppt->qs_size[ppt->index_md_scalars]; index_type++) {
      /* Allocate memory for pgb2->first_order_sources[index_type][index_tau] */
      class_alloc(pgb2->first_order_sources[index_type],
                  ppt->tau_size_quadsources * sizeof(double *),
                  ppt->error_message);
      /* Allocate memory for pgb2->first_order_sources[index_type] */
      for (int index_tau = 0; index_tau < ppt->tau_size_quadsources; index_tau++) {
          /* Loop over type and tau. For each of them, allocate memory
           for pgb2->first_order_sources[index_type][index_tau][index_k_bessel]  */
           class_alloc(pgb2->first_order_sources[index_type][index_tau],
                      pgb2->k_size_bessel * sizeof(double),
                      ppt->error_message);
    }
  }

  printf("l_size is %d \n", ptr->l_size);

  /* Allocate array for Cl[index_l][index_tau_first][index_tau_second] */
  class_alloc(pgb2->Cl, ptr->l_size[ppt->index_md_scalars] * sizeof(double **), ppt->error_message);
  for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
    class_alloc(pgb2->Cl[index_l], ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
    for (int index_tau_first = 0; index_tau_first < ppt->tau_size_quadsources; index_tau_first++) {
      class_alloc(pgb2->Cl[index_l][index_tau_first], ppt->tau_size_quadsources * sizeof(double), ppt->error_message);
    }
  }

  /* Allocate array for Cl_final[index_l][bin_first][bin_second] */
  class_alloc(pgb2->Cl_final, ptr->l_size[ppt->index_md_scalars] * sizeof(double **), ppt->error_message);
  for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
    class_alloc(pgb2->Cl_final[index_l], ppt->selection_num * sizeof(double *), ppt->error_message);
    for (int bin_first = 0; bin_first < ppt->selection_num; bin_first++) {
      class_alloc(pgb2->Cl_final[index_l][bin_first], ppt->selection_num * sizeof(double), ppt->error_message);
    }
  }

  /* Allocate and fill array for the trapezoidal weights for line of sight integration */




  class_alloc(tau0_minus_tau,
              ppt->selection_num * sizeof(double*),
              ppt->error_message);

  for ( bin = 0; bin < ppt->selection_num; bin++) {

    class_alloc(tau0_minus_tau[bin],
              pgb2->tau_size_selection * sizeof(double),
              ppt->error_message);
  }


  class_alloc(pgb2->tau_sampling_selection,
              ppt->selection_num * sizeof(double*),
              pgb2->error_message);

  for (bin = 0; bin < ppt->selection_num; bin++) {

    class_alloc(pgb2->tau_sampling_selection[bin],
              pgb2->tau_size_selection * sizeof(double),
              pgb2->error_message);
  }

  class_alloc(pgb2->tau_sampling_cls,
              pgb2->tau_size_selection * sizeof(double),
              pgb2->error_message);


  double overall_tau_min = 160000.;
  for (bin = 0; bin < ppt->selection_num; bin++) {
    //finer sampling of bins
    double tau_min;
    double tau_max;
     double z_max, z_min;
     z_max =   ppt->selection_mean[bin]+ 5. * ppt->selection_width[bin];
     z_min =  ppt->selection_mean[bin]- 5. * ppt->selection_width[bin];
     if (z_max > 1000.) {
        z_max = 1000.;
        printf('Rubbish bin');

     }
     if (z_min < 0.){
       z_min = 0.;
     }
    class_call(background_tau_of_z(
                            pba,
                            z_max,
                            &tau_min
                          ),ppt->error_message,ppt->error_message);
    class_call(background_tau_of_z(
                            pba,
                            z_min,
                            &tau_max
                          ),ppt->error_message,pgb2->error_message);

    if (tau_max > ppt->tau_sampling_quadsources[ppt->tau_size_quadsources-1]) {
      tau_max = ppt->tau_sampling_quadsources[ppt->tau_size_quadsources-1];
    }


    if (tau_min < 400.) {
      printf("Rubbish, Window function extends to before recombination\n");
      tau_min = 400.;
    }

    overall_tau_min = MIN(tau_min,overall_tau_min);

    for (index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
      pgb2->tau_sampling_selection[bin][index_tau] = tau_min + index_tau*(tau_max-tau_min)/(pgb2->tau_size_selection-1);
//      printf("pgb2->tau_sampling_selection[%d][%g] = ??\n", bin, pgb2->tau_sampling_selection[bin][index_tau] );
    }
  }

  for (index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
    pgb2->tau_sampling_cls[index_tau] = overall_tau_min + index_tau*(pba->conformal_age-overall_tau_min)/(pgb2->tau_size_selection-1);
  }

  /* New Bessel k-sampling to capture features of the oscillations */

  class_alloc(pgb2->k_bessel,
              pgb2->k_size_bessel * sizeof(double),
              ppt->error_message);

  // NOTE: this part may need to be more general (i.e. not fixed to scalars)

  double k_min = ppt->k[ppt->index_md_scalars][0];

  double k_max =  ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1];

  for (int i = 0; i < pgb2->k_size_bessel; i++) {

    pgb2->k_bessel[i] = k_min + i*(k_max-k_min)/pgb2->k_size_bessel;

    //printf("k_bessel[%d] =  %g\n", i, pgb2->k_bessel[i]);
  }

  /* Interpolate the Spherical Bessel Function */






  /* Allocate and fill array for the trapezoidal weights for chi integration */
  class_alloc(w_trapz,
              ppt->selection_num * sizeof(double*),
              ppt->error_message);


  for ( bin = 0; bin < ppt->selection_num; bin++) {
    class_alloc(w_trapz[bin],
                pgb2->tau_size_selection * sizeof(double),
                ppt->error_message);
  }


  /* Fill the array ptw->tau0_minus_tau[index_tau] */
  for (bin = 0; bin < ppt->selection_num ; bin++) {

    for(int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
      tau0_minus_tau[bin][index_tau] = tau0 - pgb2->tau_sampling_selection[bin][index_tau];
  //    printf("reversed time  = %g = ??\n", tau0_minus_tau[bin][index_tau]  );
    }

    class_call(array_trapezoidal_mweights(tau0_minus_tau[bin],
                                          pgb2->tau_size_selection,
                                          w_trapz[bin],
                                          pgb2->error_message),
                                          ppt2->error_message,
                                          ppt2->error_message);

  }


  /* Declaration of temporary pointer */
  double ** selection;


  /* Allocation of first dimension selection[bin] */
  class_alloc(selection,
              ppt->selection_num * sizeof(double*),
              ppt->error_message);


  /* Allocation of second dimension selection[bin][index_tau] */
  for(int bin = 0; bin < ppt->selection_num; bin++){
    class_alloc(selection[bin],
                /*ppt->tau_size_quadsources*/ pgb2->tau_size_selection * sizeof(double),
                ppt->error_message);

  /* transfer_selection_compute writes in to selection[bin] */
  class_call(transfer_selection_compute(ppr, pba, ppt, ptr, selection[bin], tau0_minus_tau[bin], w_trapz[bin], pgb2->tau_size_selection, pvecback, tau0, bin),
             pgb2->error_message,
             pgb2->error_message);
  //  for (i = 0; i < pgb2->tau_size_selection; i++) {
  //    printf("selection[bin = %d][pgb2->tau_sampling_selection[bin][index_tau] = %g] %g\n", bin, pgb2->tau_sampling_selection[bin][i], selection[bin][i]);
  //  }

  }

  double t1,t2,k1;
  int index_k1, index_t1, index_t2;
  int index_k_bessel;

  printf("ppt->tau_sampling_quadsources has %d points sampled between (%g,%g)\n", ppt->tau_size_quadsources, ppt->tau_sampling_quadsources[0], ppt->tau_sampling_quadsources[ppt->tau_size_quadsources-1] );
  printf("pgb2->tau_sampling_cls has %d points sampled between (%g,%g)\n", pgb2->tau_size_selection, pgb2->tau_sampling_cls[0], pgb2->tau_sampling_cls[pgb2->tau_size_selection-1] );
  printf("pgb2->tau_sampling_selection has %d points sampled between (%g,%g)\n", pgb2->tau_size_selection, pgb2->tau_sampling_selection[0], pgb2->tau_sampling_selection[pgb2->tau_size_selection-1]);

/*
  for (double k = 0; k < 10; k += 0.25) {

    class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);
    class_call(index_of_k_bessel(k, &index_k_bessel, pgb2), pgb2->error_message, pgb2->error_message);

    printf("k = %g, k_original = %g, k_bessel = %g\n",k, ppt->k[ppt->index_md_scalars][index_k],pgb2->k_bessel[index_k_bessel]);
  }







  k1 = 5.156;

  class_call(index_of_k(k1, &index_k_bessel, pgb2), pgb2->error_message, pgb2->error_message);

  printf("k_original = %g, k_bessel = %g\n", k1, pgb2->k_bessel[index_k_bessel]);*/

  /*for (int index_k = 0; index_k < ppt->k_size[ppt->index_md_scalars]; index_k++) {
    k1 = ppt->k[ppt->index_md_scalars][index_k];

    class_call(index_of_k(k1, &index_k_bessel, pgb2), pgb2->error_message, pgb2->error_message);

    printf("k_original = %g, k_bessel = %g\n", k1, pgb2->k_bessel[index_k_bessel]);
  }*/

  /*for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
    printf("pgb2->k_bessel[%d] = %g\n", index_k_bessel, pgb2->k_bessel[index_k_bessel]);
  }*/

  printf("ppt->k_size[ppt->index_md_scalars] = %d\n",ppt->k_size[ppt->index_md_scalars] );
  printf("ppt->tau_size_quadsources = %d\n", ppt->tau_size_quadsources);
  printf("pgb2->tau_size_selection = %d\n", pgb2->tau_size_selection );
  printf("pgb2->k_bessel[0] = %g\n", pgb2->k_bessel[0] );




  //printf("%g\n",ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][100 * ppt->k_size[ppt->index_md_scalars] + 20]);

  /* Now fill and interpolate this new pointer-array pgb2->first_order_sources[index_type][index_tau][index_k]
    with information from the pre computed ppt->quadsources[index_md][index_ic*ppt->tp_size[index_md]+index_type]
    [index_tau * ppt->k_size[index_md] + index_k]. */
  int dump = 0;
  double f,g;
  int i2;
  int index;
  double tau;

  double intermediate,intermediate_plus;
  //for (int index_type = 0; index_type < ppt->tp_size[index_md]; index_type++) {
    for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){


      for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
        index = 0;
        index_k = 0;

        double tau = pgb2->tau_sampling_cls[index_tau];

        class_call(index_of_tau_sampling_quadsources(tau, &index, ppt),ppt->error_message,pgb2->error_message);

        double k = pgb2->k_bessel[index_k_bessel];

        class_call(index_of_k(k, &index_k, ppt), ppt->error_message, pgb2->error_message);


        f = (pgb2->tau_sampling_cls[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

        intermediate  = (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k]+
            (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k]);

        intermediate_plus =  (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k+1]+
            (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k+1]);

        //printf("quadsources = %g\n");
        //printf("f = %g\n",f);


        g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

        pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau][index_k_bessel] = g*intermediate_plus +(1-g)*intermediate;




         //printf("pgb2->k_bessel[%d], ppt->k[][%d] = (%g, %g)\n",index_k1, index_k, pgb2->k_bessel[index_k1], ppt->k[ppt->index_md_scalars][index_k]);
         //printf("pgb2->tau_sampling_cls[%d], ppt->tau_sampling_quadsources[%d] = (%g, %g)\n", index_t1, index, pgb2->tau_sampling_cls[index_t1], ppt->tau_sampling_quadsources[index]  );

      }


        //printf("g = %g\n",g);




  /*    printf("first_order_sources(k = %g, tau = %g), quadsources(k = %d, tau = %d) = (%g, %g)\n",
        pgb2->k_bessel[index_k1],
        pgb2->tau_sampling_cls[index_t1],
        ppt->k[ppt->index_md_scalars][25],
        ppt->tau_sampling_quadsources[50],
        pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_t1][index_k1],
        ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][50 * ppt->k_size[ppt->index_md_scalars] + 25]);

    }*/



    }




  printf("First order source array filled.\n");



  printf("before k-int\n");

  /* Integrate over our integrand w.r.t. k*/
  //printf("the integral of k^2 between [%g,%g] is %g\n", ppt->k[ppt->index_md_scalars][0],ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1], integral(pba, pbs, pgb2, ppt, ppt->index_qs_delta_cdm, index_tau_first, index_tau_second, index_l));
  /* Write the result of the angular power spectrum integral into an array */


  /*class_alloc(pgb2->w_trapz_k,
              ppt->k_size[ppt->index_md_scalars] * sizeof(double),
              ppt->error_message);*/

  class_alloc(pgb2->w_trapz_k,
              pgb2->k_size_bessel * sizeof(double),
              ppt->error_message);


  /*class_call(array_trapezoidal_weights(ppt->k[ppt->index_md_scalars],
                                       ppt->k_size[ppt->index_md_scalars],
                                       pgb2->w_trapz_k,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);*/

   class_call(array_trapezoidal_weights(pgb2->k_bessel,
                                        pgb2->k_size_bessel,
                                        pgb2->w_trapz_k,
                                        pgb2->error_message),
                                        pgb2->error_message,
                                        pgb2->error_message);


  /*for (int index_k = 0; index_k < ppt->k_size[ppt->index_md_scalars]; index_k++) {
    printf("pgb2->w_trapz_k[index_k] = %g\n", pgb2->w_trapz_k[index_k]);
  }
*/

   /*class_call(background_at_tau(pba,
                                pgb2->tau_sampling_cls[10],
                                pba->long_info,
                                pba->inter_normal,
                                &dump,
                                pvecback),
              pba->error_message,
              pgb2->error_message);*/
              //herehere
  double integ;
  double result2;


  printf("k lower bound = %g\n",ppt->k[ppt->index_md_scalars][0]);
  printf("k upper bound = %g\n",ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1] );

  /*class_call(integral(pba,
                      ppm,
                      pbs,
                      pgb2,
                      ppt,
                      ptr,
                      ppt->index_tp_delta_m,
                      0,
                      0,
                      0,
                      &result2),
                      pgb2->error_message,
                      pgb2->error_message);*/



  double p1,p2,f1,f2,j1,j2;
  int k_test_index;

  printf("BEFORE LOOP\n");
  printf("pbs->x_max = %g\n",pbs->x_max );

  //class_call(index_of_k(0.001,&k_test_index,pgb2),pgb2->error_message,pgb2->error_message);
  double tmp = 0;

  for(int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++){

    double x1 = 1000.*pgb2->k_bessel[index_k_bessel];
    double x2 = pgb2->k_bessel[index_k_bessel];

    if(x1>pbs->x_max || x2 > pbs->x_max){
    /*printf("ALERT! x1= %g x2 = %g \n",x1,x2);*/continue;}
    //printf(" x1= %g x2 = %g \n",x1,x2);
    class_call(bessel_at_x(pbs,	x1, 6, &j1), pbs->error_message, pgb2->error_message);
    class_call(bessel_at_x(pbs, x2 , 0, &j2), pbs->error_message, pgb2->error_message);


    tmp +=  (1/(pgb2->k_bessel[index_k_bessel]*pgb2->k_bessel[index_k_bessel])) * j1 * j2 * pgb2->w_trapz_k[index_k_bessel];


    //printf("tmp in k loop outside integral function = %g\n", tmp);
    //printf("w_trapz_k = %g \n", pgb2->w_trapz_k[index_k_bessel]);

  }


  double m1 = 0.001;
//  double m2 = ppt->k[ppt->index_md_scalars][0];



  //class_call(bessel_at_x(pbs,	m2, 0, &f2), pbs->error_message, pgb2->error_message);

/*  for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){

    for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++){

      for(int index_k = 0; index_k < ppt->k_size[ppt->index_md_scalars]; index_k++){

        double x1 = ppt->k[ppt->index_md_scalars][index_k]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first]);

        if(x1>pbs->x_max){continue;}

        class_call(bessel_at_x(pbs,	x1, index_l, &f1), pbs->error_message, pgb2->error_message);

        printf("The result of j_%d(%g)*j_%d(%g) is %g\n",ptr->l[index_l], x1, ptr->l[index_l],x1, f1*f1);
      }
    }
  }*/
/*
  for(int index_k = 0; index_k < ppt->k_size[ppt->index_md_scalars]; index_k++){
    printf("k[%d] = %g\n", index_k, ppt->k[ppt->index_md_scalars][index_k]);
  }
*/






  //printf("The result of integral of j_%d(k)*j_%d(k) from %g to %g  is  %g (call of integral)\n", ptr->l[0], ptr->l[0], ppt->k[ppt->index_md_scalars][0], ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1],result2);
  printf("The result of integral of j_%d(1000*k)*j_%d(k)*(1/k*k) from %g to %g  is  %g (loop outside of function)\n", ptr->l[6], ptr->l[0], pgb2->k_bessel[k_test_index], ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1],tmp);





  for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
    //pbs->l_size=ptr->l_size[ppt->index_md_scalars];
    for (int k = 0; k< index_l+1; k++){
      printf("X");
    }
    printf("\n");

    for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++){
      for(int index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++){

        class_call(integral(pba,ppm, pbs, pgb2, ppt, ptr, ppt->index_qs_delta_cdm, index_tau_first, index_tau_second, index_l, &integ),
                    pgb2->error_message,
                    pgb2->error_message);
        pgb2->Cl[index_l][index_tau_first][index_tau_second] = integ;


      //  printf("selection[0][%d] = %g\n", index_tau_first, selection[0][index_tau_first]);
      //  printf("Cl[index_l = %d][index_tau_first = %d][index_tau_second = %d] = %g\n",index_l, index_tau_first, index_tau_second, pgb2->Cl[index_l][index_tau_first][index_tau_second]);
        //printf("%g\n", integrand(pba,pbs,pgb2, ppt, index_type, index_tau_first, index_tau_second, index_k, index_l));
        //printf("integral(index_l = %d, index_tau_first = %d, index_tau_second = %d) = %g\n",index_l,index_tau_first, index_tau_second, integral(pba, pbs, pgb2, ppt, ppt->index_qs_delta_cdm, index_tau_first, index_tau_second, index_l));

      }
    }

    //printf("Cl[index_l = %d][index_tau_first = %d][index_tau_second = %d] = %g\n",index_l, index_tau_first, index_tau_second, pgb2->Cl[index_l][index_tau_first][index_tau_second]);
  }
  free(pgb2->w_trapz_k);
  printf("Done k-integration.\n");


  int  bin1 = 0;
  int  bin2 = 0;
  double tmp2 = 0;
  int index_tau1;
  int index_tau2;
  double temp, temp_minus, temp_plus;

  for(index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
  //  for (bin1 = 0; bin1 < ppt->selection_num; bin1++){
    //  for(bin2 = 0; bin2 < ppt->selection_num; bin2++){
    pgb2->Cl_final[index_l][0][0] = 0.;


    for(index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++){
      //printf("index_tau2 = %d of %d\n", index_tau_second,  pgb2->tau_size_selection );

      double temp123 = 0.;
      double tau2 = pgb2->tau_sampling_selection[0][index_tau_second];


      for(index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++){

        double tau1 = pgb2->tau_sampling_selection[0][index_tau_first];

        index_of_tau_sampling_cls(tau1, &index_tau1, pgb2);
        index_of_tau_sampling_cls(tau2, &index_tau2, pgb2);

        temp_minus = pgb2->Cl[index_l][index_tau1-1][index_tau2-1]*(pgb2->tau_sampling_cls[index_tau1]-tau1)
                      + pgb2->Cl[index_l][index_tau1][index_tau2-1]*(tau1-pgb2->tau_sampling_cls[index_tau1-1]);
        temp_minus /= (pgb2->tau_sampling_cls[index_tau1] - pgb2->tau_sampling_cls[index_tau1-1]);

        temp_plus = pgb2->Cl[index_l][index_tau1-1][index_tau2]*(pgb2->tau_sampling_cls[index_tau1]-tau1)
                      + pgb2->Cl[index_l][index_tau1][index_tau2]*(tau1-pgb2->tau_sampling_cls[index_tau1-1]);
        temp_plus /= (pgb2->tau_sampling_cls[index_tau1] - pgb2->tau_sampling_cls[index_tau1-1]);

        temp = temp_minus * (pgb2->tau_sampling_cls[index_tau2]-tau2)
                      + temp_plus * (tau2-pgb2->tau_sampling_cls[index_tau2-1]);

        temp /= (pgb2->tau_sampling_cls[index_tau2] - pgb2->tau_sampling_cls[index_tau2-1]);

        temp123 += temp * w_trapz[0][index_tau_first] * selection[0][index_tau_first];

      }


      //printf("selection[0][%g] = %g\n",pgb2->tau_sampling_selection[0][index_tau_second], selection[0][index_tau_second] );
    //  printf("temp %g temp123 %g \n",temp,temp123 );

      pgb2->Cl_final[index_l][0][0] += temp123 * w_trapz[0][index_tau_second]
          * selection[0][index_tau_second];

    }
    //  }
  //  }
  }

  for  (index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
    printf("l*(l+1)*Cl_final[%d][0][0] = %g\n",ptr->l[index_l], pgb2->Cl_final[index_l][0][0]*ptr->l[index_l]*(ptr->l[index_l]+1.));
  }


  printf("We are here 20!\n");

  // ================================================================================
  // =                          Create sources directory                            =
  // ================================================================================


  /* Create the files to store the source functions in */
  if ((ppr2->store_sources_to_disk == _TRUE_) || (ppr2->load_sources_from_disk == _TRUE_)) {

    /* We are going to store the sources in n=k_size files, one for each requested k1 */
    class_alloc (ppt2->sources_files, ppt2->k_size*sizeof(FILE *), ppt2->error_message);
    class_alloc (ppt2->sources_paths, ppt2->k_size*sizeof(char *), ppt2->error_message);

    for (int index_k1=0; index_k1<ppt2->k_size; ++index_k1) {

      /* The name of each sources file will have the k1 index in it */
      class_alloc (ppt2->sources_paths[index_k1], _FILENAMESIZE_*sizeof(char), ppt2->error_message);
      sprintf (ppt2->sources_paths[index_k1], "%s/sources_%03d.dat", ppt2->sources_dir, index_k1);

    } // end of loop on index_k1

    if (ppr2->store_sources_to_disk == _TRUE_)
      if (ppt2->perturbations2_verbose > 2)
        printf ("     * will create %d files to store the sources\n", ppt2->k_size);

  } // end of if(ppr2->store_sources_to_disk)

printf("We are here 21!\n");


  // ================================================================================
  // =                          Print result to terminal                            =
  // ================================================================================
  int l_max = 200;

/*
  for(index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
    for (index_tau_first = 0; index_tau_first < ppt->tau_size; index_tau_first++){
      for(index_tau_second = 0; index_tau_second < ppt->tau_size; index_tau_second++){
        printf("Test dens-dens angular power spectrum C_%d(%g)(%g) \n", ptr->l[index_l] , ppt->tau_sampling[index_tau_first], ppt->tau_sampling[index_tau_second],
          pgb2->Cl_final[index_l][index_tau_first][index_tau_second]);
      }
    }
  }
*/

free(pgb2->tau_sampling_cls);
printf("End of galbispectra2!\n");
  return _SUCCESS_;

} // end of galbispectra2_init()
