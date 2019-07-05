/* Module to compute galaxy bispectra */


#include "galbispectra2.h"
#include "perturbations2.h"


/* bessel_at_x_second_deriv() - given the redshift and the index_l of the spherical Bessel function, this function
  will compute its second derivative using a recursion relation, storing it as * result.


* @param galbispectra2          Input: galaxy bispectra structre
* @param bessels                Input: Bessel structure
* @param index_l                Input: l
* @param result                 Output: j''_l(x) (second derivative of the spherical Bessel function)
*/

double bessel_at_x_second_deriv(struct galbispectra2 * pgb2,
                       struct bessels * pbs,
                       double z,
                       int index_l,
                       double * result){

    double out;
    double j1;
    double j2;

    if (pbs->l[index_l]==0){
      out = 0.0;
    }
    else if(pbs->l[index_l]==1){
      out = 0.0;
    }
    /*using recurrence relation */
    else
    class_call(bessel_at_x(pbs, z, index_l, &j1), pbs->error_message, pgb2->error_message);
    class_call(bessel_j(pbs, pbs->l[index_l], z, &j2), pbs->error_message, pgb2->error_message);

    out = ((pbs->l[index_l]*pbs->l[index_l]-pbs->l[index_l]-z*z)*j1+2*z*j2)/(z*z);

    *result = out;
    return _SUCCESS_;
}

/* integral() - given two types of number count perturbations, this function will integrate over k
  the integral is over Pk, two Bessel functions, two transfer functions and a prefactor, ready to be
  integrated over time for subsequent angular power/bi spectra.


* @param background             Input: background structure
* @param primordial             Input: primordial structure
* @param bessels                Input: Bessel structure
* @param galbispectra2          Input: galaxy bispectra structre
* @param perturbs               Input: perturbations structure
* @param transfers              Input: transfer structure
* @param index_type_first       Input: perturbation type
* @param index_type_second      Input: perturbation type
* @param index_tau_first        Input: conformal time
* @param index_tau_second       Input: conformal time
* @param index_l                Input: l
* @param pvecback               Output: temporary background vector
* @param result                 Output: result of integration
*/

int integral(
  struct background * pba,
  struct primordial * ppm,
  struct bessels * pbs,
  struct galbispectra2 * pgb2,
  struct perturbs * ppt,
  struct transfers * ptr,
  int index_type_first,
  int index_type_second,
  int index_tau_first,
  int index_tau_second,
  int index_l,
  double * pvecback,
  double * result)
  {
  int index_md;
  double tmp = 0.;
  double lensing_result = 0.;

  double j1;
  double j2;
  double Pk;
  double type1;
  double type2;

  /* Pair input type with corresponding source type. Each type is made up of, the transfer function of the perturbation, a prefactor
    and a Bessel function (or a derivative of). For integrated terms inside pgb2->first_order_sources_integ, the spherical Bessel function
    may have already been integrated over before this stage, thus is will not appear explicitly here. */
  int last_index1 = 0;
  int last_index2 = 0;
  int index_tau_lens;
  double prefactor1;
  double prefactor2;
  double lensing_result1;
  double lensing_result2;

  for(int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++){

    double x1 = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first]);
    double x2 = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second]);

    double r1 = pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first];
    double r2 = pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second];

    if(x1>pbs->x_max || x2 > pbs->x_max){
      printf("ALERT! x1= %g x2 = %g \n",x1,x2);continue;}

    class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, pgb2->k_bessel[index_k_bessel], &Pk), ppm->error_message, pgb2->error_message);

    /* First Type: Density */

    if( index_type_first == pgb2->index_type_delta_cdm){

      class_call(bessel_at_x(pbs, x1 , index_l, &j1), pbs->error_message, pgb2->error_message);

      type1 = pgb2->first_order_sources[pgb2->index_type_delta_cdm][index_tau_first][index_k_bessel] * j1;

    }

    /* First Type: RSD */

    else if( index_type_first == pgb2->index_type_rsd ){

      class_call(bessel_at_x_second_deriv(pgb2, pbs, x1 , index_l, &j1), pbs->error_message, pgb2->error_message);

      class_call(background_at_tau(pba,
                                   pgb2->tau_sampling_cls[index_tau_first],
                                   pba->long_info,
                                   pba->inter_normal,
                                   &last_index1,
                                   pvecback),
                                   pba->error_message,
                                   pgb2->error_message);

      prefactor1 = -1.0/pvecback[pba->index_bg_H];

      type1 = prefactor1 * pgb2->k_bessel[index_k_bessel] * pgb2->k_bessel[index_k_bessel]
                * pgb2->first_order_sources[pgb2->index_type_rsd][index_tau_first][index_k_bessel] * j1;
    }

    /* First Type: Lensing convergence */

    else if( index_type_first == pgb2->index_type_lens ){
      int index_tau_source1;
      int last_index1;
      int index_tau_lens1;
      double r_lens1;

      prefactor1 = -ptr->l[index_l] * (ptr->l[index_l] + 1.);

      type1 = prefactor1 * pgb2->first_order_sources_integ[pgb2->index_type_lens][index_l][index_tau_first][index_k_bessel];
    }

    else {
      type1 = 0.;
    }

    /* Establish terms that are called by index_type_second */

    /* Second Type: Density */

    if( index_type_second == pgb2->index_type_delta_cdm){

      class_call(bessel_at_x(pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);

      type2 = pgb2->first_order_sources[pgb2->index_type_delta_cdm][index_tau_first][index_k_bessel] * j2;
    }

    /* Second Type: RSD */

    else if( index_type_second == pgb2->index_type_rsd ){

      class_call(bessel_at_x_second_deriv(pgb2, pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);

      class_call(background_at_tau(pba,
                                   pgb2->tau_sampling_cls[index_tau_second],
                                   pba->long_info,
                                   pba->inter_normal,
                                   &last_index2,
                                   pvecback),
                                   pba->error_message,
                                   pgb2->error_message);

      prefactor2 = -1/pvecback[pba->index_bg_H];

      type2 = prefactor2 * pgb2->k_bessel[index_k_bessel] * pgb2->k_bessel[index_k_bessel]
                * pgb2->first_order_sources[pgb2->index_type_rsd][index_tau_second][index_k_bessel] * j2;
    }

    /* Second Type: Lensing convergence */

    else if( index_type_second == pgb2->index_type_lens ){

      prefactor2 = -ptr->l[index_l] * (ptr->l[index_l] + 1.);

      type2 = prefactor2 * pgb2->first_order_sources_integ[pgb2->index_type_lens][index_l][index_tau_second][index_k_bessel];
    }

    else {
      type2 = 0.;
    }

    // NOTE: Should have some CLASS test here to check whether the types called are valid.



    /* Using area of trapezoid, we can sum a number of areas of trapezoids to approximate the integral */
    tmp += pow(pgb2->k_bessel[index_k_bessel],-1.0) * 4. * _PI_ *  Pk * type1 * type2  * pgb2->w_trapz_k[index_k_bessel];
    }


  *result = tmp;
  return _SUCCESS_;
}



/***************************************************************************
====================    Define Search Functions ============================
****************************************************************************/

/* index_of_tau_sampling_cls() - given some conformal time tau1, this function will output the (closest)
  corresponding index in the pgb2->tau_sampling_cls[index_tau1] array.

* @param tau1                     Input: conformal time
* @param index_tau1               Output: index in the array tau_sampling_cls
* @param galbispectra2            Input: galbispectra2 structure
*/


int index_of_tau_sampling_cls(double tau1,
                 int * index_tau1,
                 int * last_index,
                 struct galbispectra2 * pgb2){



                 double output = (tau1 - pgb2->tau_sampling_cls[0])/(pgb2->tau_sampling_cls[pgb2->tau_size_selection-1] - pgb2->tau_sampling_cls[0]) * (pgb2->tau_size_selection-1);
                 * index_tau1 = (int) ceil(output);
                 if (* index_tau1 < 1) {
                   * index_tau1 = 1;
                 }
                 else if(* index_tau1 > pgb2->tau_size_selection -1){
                   * index_tau1 = pgb2->tau_size_selection -1;

                 }

               }

/* index_of_tau_sampling_ppt2() - given some conformal time tau1, this function will output the (closest)
 corresponding index in the ppt2->tau_sampling[index_tau1] array.

* @param tau1                     Input: conformal time
* @param index_tau1               Output: index in the array tau_sampling_cls
* @param perturbs                 Input: perturbs structure
*/


int index_of_tau_sampling_ppt2(double tau1,
                 int * index_tau1,
                 int * last_index,
                 struct perturbs2 * ppt2){

                 * index_tau1 = 0;
                 * last_index = 0;
                 double tau;
/* Scan through the ppt2->tau_sampling grid and assign the index which gives the best estimate of tau1 and tau2 */

                    for (int index = *last_index; index < ppt2->tau_size ; index++) {
                      tau = ppt2->tau_sampling[index];

                      if (tau > tau1 && *index_tau1 == 0) {
                        *index_tau1=index;
                      }



                      if (*index_tau1 != 0){
                        *last_index = *index_tau1;
                        break;
                      }
                    }
                  }

/* index_of_k() - function to output indices of ppt->k to a given input k value

* @param k                    Input: k-value
* @param index_tau1           Output: index in the array ppt->k
* @param perturbs             Input: perturbs structure
*/

int index_of_k_old(double k,
               int * index_k,
               int * last_index_k,
               struct perturbs * ppt){
  double k_start;
  double k_end;

  /* Initialise index_k to a value that does not exist within the grid, this will circumvent any chance of the
    wrong index being assigned */
  *index_k = -1;
  if (*last_index_k > 0){

    k_start = ppt->k[ppt->index_md_scalars][*last_index_k - 1];
    k_end = ppt->k[ppt->index_md_scalars][*last_index_k];

    if(k_end == k){
      *index_k=*last_index_k;
    }
  }

  else{
    k_start = ppt->k[ppt->index_md_scalars][* last_index_k];
  }

  double k_ppt;

  /*if (k == k_start){
    printf("k == k3_start\n");
    if (*last_index_k > 0){
      *index_k=*last_index_k-1;
    }

    if (*last_index_k == 0){
      *index_k=*last_index_k;
    }

  }*/

  /* If k< k_start, search grid from index = 0 */
  if (k < k_start){
    for (int index = 0; index < *last_index_k; index++) {
      k_ppt = ppt->k[ppt->index_md_scalars][index];

      if (k < k_ppt && *index_k == -1) {
        *index_k=index;
        *last_index_k=index;
        break;
      }

    }
  }

  else {

    for (int index = *last_index_k; index < ppt->k_size[ppt->index_md_scalars]; index++) {
     k_ppt = ppt->k[ppt->index_md_scalars][index];

     if (k <= k_ppt && *index_k == -1) {
       *index_k=index;
       *last_index_k=index;
       break;
     }
    }
  }

  if (*index_k == -1){
    printf("%s\n",'Search out of bounds in index_of_k()' );
  }


  return _SUCCESS_;
}

int index_of_k3_old(double k,
  int index_k1,
  int index_k2,
  int * index_k3,
  int * last_index_k3,
  struct perturbs2 * ppt2){
  double k3_start;
  double k_max;
  double k3_ppt_last;
  double k3_ppt;


  *index_k3 = -1;
  if (*last_index_k3 > 0){

    k3_start = ppt2->k3[index_k1][index_k2][*last_index_k3 - 1];
  }

  if (*last_index_k3 == 0){

    k3_start = ppt2->k3[index_k1][index_k2][* last_index_k3];
  }



  /*if (k == k3_start){
    printf("k == k3_start\n");
    if (*last_index_k3 > 0){
      *index_k3=*last_index_k3-1;
    }

    if (*last_index_k3 == 0){
      *index_k3=*last_index_k3;
    }

  }*/



  if (k < k3_start){ // need to search from beginning
    //printf("search for backward  k = %g, k3_start = %g\n",k, k3_start);

    for (int index = 0; index < *last_index_k3; index++) {
      k3_ppt = ppt2->k3[index_k1][index_k2][index];
      //printf("%g\n", k_ppt);
      if (k < k3_ppt && *index_k3 == -1) {
        *index_k3=index;
        *last_index_k3=index;
        break;
      }
    }
  }

  //k_max = ppt2->k3[index_k1][index_k2][ppt2->k3_size[index_k1][index_k2]-1];

  /*if (k >= k_max){
    /* then take *index_k3 to be the index that corresponds to the largest value in ppt2->k3[index_k1][index_k2][index_k3]*/

  /*  for(int index = 1; index < ppt2->k3_size[index_k1][index_k2]; index++){
      k3_ppt = ppt2->k3[index_k1][index_k2][index];
      k3_ppt_last = ppt2->k3[index_k1][index_k2][index-1];

      if (abs(k3_ppt) > abs(k3_ppt_last)) {
        *index_k3=index;
      }
    }
    printf("input k is larger than highest k on grid, index set to maximum.\n");
  }*/



  else {
    for (int index = *last_index_k3; index < ppt2->k3_size[index_k1][index_k2]; index++) {
     k3_ppt = ppt2->k3[index_k1][index_k2][index];

     //printf("%g\n", k_ppt);
     if (k <= k3_ppt && *index_k3 == -1){
       *index_k3=index;
       *last_index_k3=index;
       //printf("found %d \n",*index_k3);
       break;
     }
    }
  }

  if (*index_k3 == -1){
    //printf("%s\n",'Search out of bounds in index_of_k()' );
  }

/*  if (*index_k3 != -1){
    printf("k3 = %g, k3_ppt = %g \n",k, k3_ppt = ppt2->k3[index_k1][index_k2][*index_k3]);
  }*/



  return _SUCCESS_;

}

int index_of_k3_mod(double k,
  int index_k1,
  int index_k2,
  int * index_k3,
  int * last_index_k3,
  struct perturbs2 * ppt2){

  *index_k3 = 0;

  double a;
  double b;

  for (int index = 1; index < ppt2->k3[index_k1][index_k2]; index++) {

    a = abs(k-ppt2->k3[index_k1][index_k2][index]);

    b = abs(k-ppt2->k3[index_k1][index_k2][index-1]);

    if (a < b) {
      *index_k3=index;
    }
  }


  return _SUCCESS_;
}



int index_of_k(double k,
               int * index_k,
               struct perturbs * ppt){
  double k_start;
  double k_ppt;

  k_start = ppt->k[ppt->index_md_scalars][*index_k];
  if(k<k_start){printf("THING THAT SHOULD NOT HAPPEN HAPPENED\n");exit(2);}
  for (int index = *index_k+1; index < ppt->k_size[ppt->index_md_scalars]; index++) {
   k_ppt = ppt->k[ppt->index_md_scalars][index];
   if (k_ppt>k) {
     *index_k=index-1;
     break;
   }
  }

  return _SUCCESS_;
}

/* For a given value of k, this function finds the index in the ppt2->k[index] grid. Both k1 and k2 live on this grid. */
int index_of_k1(double k,
    int * index_k,
    struct perturbs2 * ppt2){
    double k_start;
    double k_ppt;

  k_start = ppt2->k[*index_k];
  if(k<k_start){printf("THING THAT SHOULD NOT HAPPEN HAPPENED\n");exit(2);}
  for (int index = *index_k+1; index < ppt2->k_size; index++) {
   k_ppt = ppt2->k[index];
   if (k_ppt>k) {
     *index_k=index-1;
     break;
   }
  }

  return _SUCCESS_;
}


/* For a given value of k, index_k1, index_k2, this function finds the index_k3 in the ppt2->k2[index_k1][index_k2][index_k3] grid. */
int index_of_k3(double k,
    int index_k1,
    int index_k2,
    int * index_k3,
    struct perturbs2 * ppt2){
    double k_start;
    double k_ppt;

  k_start = ppt2->k3[index_k1][index_k2][*index_k3];
  if(k<k_start){printf("THING THAT SHOULD NOT HAPPEN HAPPENED\n");exit(2);}
  for (int index = *index_k3+1; index < ppt2->k3_size[index_k1][index_k2]; index++) {
    k_ppt = ppt2->k3[index_k1][index_k2][index];
    if (k_ppt>k) {
     *index_k3=index-1;
     break;
   }
  }

  return _SUCCESS_;


}




/* index_of_k_bessel() - function to output indices of ppt->k to a given input k value

* @param k                    Input: k-value
* @param index_tau1           Output: index in the array pgb2->k_bessel
* @param perturbs             Input: perturbs structure
*/

int index_of_k_bessel(double k,
               int * index_k_bessel,
               int * last_index,
               struct galbispectra2 * pgb2){
               double k_start;

               * index_k_bessel = -1;

               if (*last_index > 0){
                 k_start = pgb2->k_bessel[*last_index - 1];
               }
               else
                 k_start = pgb2->k_bessel[* last_index];


               double k_bessel;

               /* Scan through the k_bessel grid pgb2->k_bessel[index] and assign the index which gives the closest estimate of k */
               if (k < k_start){ // need to search from beginning
                  for (int index = *last_index; index < pgb2->k_size_bessel; index++) {
                    k_bessel = pgb2->k_bessel[index];

                    if (k < k_bessel && *index_k_bessel == -1) {
                      *index_k_bessel=index;
                      break;
                    }
                  }
              }

              else {
                for (int index = *last_index; index < pgb2->k_size_bessel; index++) {

                  k_bessel = pgb2->k_bessel[index];

                  if (k < k_bessel && *index_k_bessel == -1) {
                    *index_k_bessel=index;
                    break;
                  }
                }
             }



              if (*index_k_bessel = -1){
                printf("%s\n",'Search out of bounds in index_of_k_bessel()' );
              }
              return _SUCCESS_;
    }

int k3_search_mod(double k,
                  int index_k1,
                  int index_k2,
                  int * index_mark,
                  double * a,
                  struct perturbs2 * ppt2){

  for (int index = 0; index < ppt2->k3_size[index_k1][index_k2]; index++) {

    double k3,k_diff;

    k3 = ppt2->k3[index_k1][index_k2][index];
    k_diff = k3-k;

    a[index] = k_diff;

    if (k_diff < 0) {
      a[index]=-k_diff;
    }

    if(a[index]<a[index-1]){

      *index_mark = index;
    }
  }

  return _SUCCESS_;
}

int index_of_tau_sampling_quadsources(double tau,
                 int * index_tau,
                 struct perturbs * ppt){

                   double tau_start;
                   double tau_ppt;

                   tau_start = ppt->tau_sampling_quadsources[*index_tau];
                   if(tau<tau_start){printf("THING THAT SHOULD NOT HAPPEN HAPPENED\n");exit(2);}
                   for (int index = *index_tau+1; index < ppt->tau_size_quadsources; index++) {
                    tau_ppt = ppt->tau_sampling_quadsources[index];
                    if (tau_ppt>tau) {
                      *index_tau=index-1;
                      break;
                    }
                   }

                   return _SUCCESS_;

                 }


/* galbispectra2_init - the main function of this module. Here we can compute all types of galaxy bispectra. This function
  is called in the song.c main file. Note that all of the input structures should be filled before galbispectra2_init is
  called i.e. all of the associated modules should be ran beforehand. The output of this function is stored in the galbispectra2
  structure which is poined to via the pointer pgb2.

* @param precision              Input: precision structure
* @param precision2             Input: second-order precision structure
* @param background             Input: background structure
* @param primordial             Input: primordial structure
* @param thermodynamics         Input: thermodynamics structure
* @param perturbations          Input: perturbations structure
* @param perturbations2         Input: perturbations2 structure
* @param galbispectra2          Output: galaxy bispectra structre
* @param bessels                Input: Bessel structure
*/

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
    /* Define local (stack) variables */
  int index_l;
  int index_k;
  int index_tau;
  int index_tau_ppt2;
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
  int index_k2;
  int index_k3;
  int index_tau_qs;
  int index_k1_fo;
  int index_k2_fo;
  double T1, T2, T3;
  double kernel;
  int * last_index_k3;
  double k1_fo;
  double k2_fo;
  double p,q;
  // NOTE: Alert
  /*the bessel grid is boosted by an integer factor to account for rapid oscillations. Perhaps the k grid should have a similar amount of points.*/
  int bessel_boost = 100;
  pgb2->tau_size_selection = 50;
  /* We wish to double the number of slots so that points in the two grids tau_sampling_cls and tau_sampling_bessel align. */
  pgb2->tau_size_bessel = bessel_boost * (pgb2->tau_size_selection-1) + 1;
  pgb2->k_size_bessel = 9000;  /*tau_size_bessel*/   /*k_size * boost*/
  printf("Starting galaxy bispectra module...\n");


  double * result;




  /* Fix k2 = 10e-5, find what this value is in both the ppt and ppt2 grids respectively */
  /*class_call(index_of_k1(0.00001,
                 &index_k2,
                 ppt2),
                 ppt2->error_message,
                 pgb2->error_message);

  double *** a;
  double * b;


  class_alloc(a,
        ppt2->k_size * sizeof(double**),
        pgb2->error_message);


  //printf("k = 7e-06, result = %g, ratio = %g\n", ppt2->k3[0][1][index_k3], 7.0e-06/ppt2->k3[0][1][index_k3]);

  k2_fo = ppt2->k[index_k2];

  class_call(index_of_k(k2_fo, &index_k2_fo, ppt), ppt->error_message, pgb2->error_message);

  index_tau_qs = 0;
  double tau_kernel;
  class_call(background_tau_of_z(pba, 500.0, &tau_kernel),ppt->error_message,pgb2->error_message);

  /* Choose some time slice, for example tau = ppt2->tau_sampling[5], find the corresponding index in the tau_sampling_quadsources grid */
  /*class_call(index_of_tau_sampling_quadsources(tau_kernel, &index_tau_qs, ppt), ppt->error_message, pgb2->error_message);
  int last_tau_kernel = 0;

  class_call(index_of_tau_sampling_ppt2(tau_kernel,
                 &index_tau_ppt2,
                 &last_tau_kernel,
                 pgb2), pgb2->error_message, pgb2->error_message);

  printf("ppt2->tau_sampling[5] = %g, ppt->tau_sampling_quadsources[index_tau_qs] = %g \n", ppt2->tau_sampling[5], ppt->tau_sampling_quadsources[index_tau_qs]);
  printf("ppt2->k_size = %d\n", ppt2->k_size);
  printf("ppt2->tau_size = %d\n", ppt2->tau_size);

  /*Initialise some variables */
  /*index_k3 = 0;
  index_k1_fo = 0;
  last_index_k3 = 0;
  int index_k1p = 0;


  /* Check that the tau values on both time grid are equal */
  //printf("index_tau_qs = %d, ppt->tau_sampling_quadsources[%d] = %g\n", index_tau_qs, index_tau_qs, ppt->tau_sampling_quadsources[index_tau_qs]);

  /*int * last_i_k;
  int zero = 0;
  last_i_k = &zero;
  //printf("last_i_k = %d\n", last_i_k);
  //printf("*last_i_k = %d\n", *last_i_k);

  for (int index_k1 = ppt2->k_size-1; index_k1 > 0; --index_k1) {

    class_alloc(a,
           ppt2->k3_size[index_k1][index_k2] * sizeof(double),
           ppt->error_message);
  /* Need to set k1 = k3 */
  /*  k1_fo = ppt2->k[index_k1];

    class_call(index_of_k_old(ppt2->k[index_k1], &index_k1_fo, last_i_k, ppt), pgb2->error_message, pgb2->error_message);

    class_call(k3_search_mod(ppt2->k[index_k1], index_k1, index_k2, &index_k3, a, ppt2), pgb2->error_message, pgb2->error_message);

    /* NOTE: index_k1 is looped over, index_k2 is fixed by k2=1e-05, index_k3 is searched (also fixed by k1=k2), index_k1_fo is also searched, index_k2_fo is fixed by k2=1e-15, both time indices are also fixed. */

  /*  T1 = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index_tau_qs * ppt->k_size[ppt->index_md_scalars] + index_k1_fo];

    T2 = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index_tau_qs * ppt->k_size[ppt->index_md_scalars] + index_k2_fo];

    T3 = ppt2->sources[ppt2->index_tp2_delta_cdm][index_k1][index_k2][index_tau_ppt2*ppt2->k3_size[index_k1][index_k2]+index_k3];


    kernel = T3/(T1*T2);
    //printf("%g  %g\n", ppt2->k[index_k1], kernel);

  }
  //printf("tau_kernel = %g\n",tau_kernel);
  //printf("ppt2->tau_sampling[%d] = %g, ppt->tau_sampling_quadsources[%d] = %g, ratio = %g\n",index_tau, ppt2->tau_sampling[index_tau], index_tau_qs, ppt->tau_sampling_quadsources[index_tau_qs],  ppt->tau_sampling_quadsources[index_tau_qs]/ppt2->tau_sampling[index_tau]);



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

  class_alloc(pgb2->tau_sampling_bessel,
              pgb2->tau_size_bessel * sizeof(double),
              pgb2->error_message);


  double overall_tau_min = 160000.;


  /* Define new tau_sampling_selection. This sampling is bin-dependent */


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
    }
  }

  for (index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
    pgb2->tau_sampling_cls[index_tau] = overall_tau_min + index_tau*(pba->conformal_age-overall_tau_min)/(pgb2->tau_size_selection-1);
  }

  for (int index_tau_bessel = 0; index_tau_bessel < pgb2->tau_size_bessel; index_tau_bessel++) {
    pgb2->tau_sampling_bessel[index_tau_bessel] = overall_tau_min + index_tau_bessel*(pba->conformal_age-overall_tau_min)/(pgb2->tau_size_bessel-1);
  }




  /* New Bessel k-sampling to capture features of the Bessel oscillations */

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







  /* Allocate and fill array for the trapezoidal weights for Chi integration w_trapz[bin][index_tau] */
  class_alloc(w_trapz,
              ppt->selection_num * sizeof(double*),
              ppt->error_message);


  for ( bin = 0; bin < ppt->selection_num; bin++) {
    class_alloc(w_trapz[bin],
                pgb2->tau_size_selection * sizeof(double),
                ppt->error_message);
  }

  /* Allocate and fill array for the trapezoidal weights for chi integration in the lensing term w_trapz_lens[index_tau] */
  class_alloc(pgb2->w_trapz_lens,
              pgb2->tau_size_bessel * sizeof(double*),
              ppt->error_message);

  class_call(array_trapezoidal_weights(pgb2->tau_sampling_bessel,
                                       pgb2->tau_size_bessel,
                                       pgb2->w_trapz_lens,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);

/*for (int index_tau = 0; index_tau < ppt->tau_size_quadsources; index_tau++) {
  printf("pgb2->w_trapz_lens[%d] = %g\n", index_tau, pgb2->w_trapz_lens[index_tau]);
}*/



  /* Fill the array ptw->tau0_minus_tau[index_tau] */
  for (bin = 0; bin < ppt->selection_num ; bin++) {

    for(int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
      tau0_minus_tau[bin][index_tau] = tau0 - pgb2->tau_sampling_selection[bin][index_tau];
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
  double * pvecback2;


  /* Allocation of first dimension selection[bin] */
  class_alloc(selection,
              ppt->selection_num * sizeof(double*),
              ppt->error_message);

  printf("selection_num = %d\n",ppt->selection_num);
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              ptr->error_message);

  /* Allocation of second dimension selection[bin][index_tau] */
  for(int bin = 0; bin < ppt->selection_num; bin++){
    printf("bin = %d\n", bin );
    class_alloc(selection[bin],
                pgb2->tau_size_selection * sizeof(double),
                ppt->error_message);

    printf("tau0_minus_tau[%d]= %g\n",bin, tau0_minus_tau[bin]);
    printf("w_trapz[%d] = %g\n", bin, w_trapz[bin]);
    printf("tau0 = %g\n", tau0 );
    printf("pvecback = %g\n",pvecback);
    printf("pvecback = %g\n",*pvecback);
    /* transfer_selection_compute writes in to selection[bin] */
    class_call(transfer_selection_compute(ppr,
                                          pba,
                                          ppt,
                                          ptr,
                                          selection[bin],
                                          tau0_minus_tau[bin],
                                          w_trapz[bin],
                                          pgb2->tau_size_selection,
                                          pvecback,
                                          tau0,
                                          bin),
               pgb2->error_message,
               pgb2->error_message);

  }

  double t1,t2,k1;
  int index_k1, index_t1, index_t2;
  int index_k_bessel;

  printf("ppt->tau_sampling_quadsources has %d points sampled between (%g,%g)\n", ppt->tau_size_quadsources, ppt->tau_sampling_quadsources[0], ppt->tau_sampling_quadsources[ppt->tau_size_quadsources-1] );
  printf("pgb2->tau_sampling_cls has %d points sampled between (%g,%g)\n", pgb2->tau_size_selection, pgb2->tau_sampling_cls[0], pgb2->tau_sampling_cls[pgb2->tau_size_selection-1] );
  printf("pgb2->tau_sampling_selection has %d points sampled between (%g,%g)\n", pgb2->tau_size_selection, pgb2->tau_sampling_selection[0], pgb2->tau_sampling_selection[pgb2->tau_size_selection-1]);



  printf("ppt->k_size[ppt->index_md_scalars] = %d\n",ppt->k_size[ppt->index_md_scalars] );
  printf("ppt->tau_size_quadsources = %d\n", ppt->tau_size_quadsources);
  printf("pgb2->tau_size_selection = %d\n", pgb2->tau_size_selection );
  printf("pgb2->k_bessel[0] = %g\n", pgb2->k_bessel[0]);

  int dump = 0;
  double f,g;
  int i2;
  int index; index_type;
  double tau;
  int last_index;
  int last_index_k;

  //NOTE: Fix this
  index_type = 0;
  double k5 = 5.0;

/* NOTE: The following if statements should be dependent on the unser input */
// if want delta_cdm
/*  if (k5 == 5.0){
    pgb2->index_type_delta_cdm = index_type;
    index_type++;
  }
//want rsd

  if (k5 == 5.0){
    pgb2->index_type_rsd = index_type;
    index_type++;
  }*/

  if (k5 == 5.0){
    pgb2->index_type_phi_plus_psi = index_type;
    index_type++;
  }

  if (k5 == 5.0){
    pgb2->index_type_lens = index_type;
    index_type++;
  }

  pgb2->type_size = index_type;
  printf("type size = %d\n", pgb2->type_size );

  /* Define an array of values of first order transfer functions:
              pgb2->first_order_sources[index_type][index_tau][index_k_bessel] */



  class_alloc(pgb2->first_order_sources, pgb2->type_size * sizeof(double **), ppt->error_message);
    for (int index_type = 0; index_type < pgb2->type_size; index_type++) {
      /* Allocate memory for pgb2->first_order_sources[index_type][index_tau] */
      class_alloc(pgb2->first_order_sources[index_type],
                  pgb2->tau_size_selection * sizeof(double *),
                  ppt->error_message);
      /* Allocate memory for pgb2->first_order_sources[index_type] */
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
          /* Loop over type and tau. For each of them, allocate memory
           for pgb2->first_order_sources[index_type][index_tau][index_k_bessel]  */
         class_alloc(pgb2->first_order_sources[index_type][index_tau],
                      pgb2->k_size_bessel * sizeof(double),
                      ppt->error_message);
    }
  }

  class_alloc(pgb2->first_order_sources_integrand, pgb2->type_size * sizeof(double **), ppt->error_message);
    for (int index_type = 0; index_type < pgb2->type_size; index_type++) {
      /* Allocate memory for pgb2->first_order_sources_integrand[index_type][index_tau] */
      class_alloc(pgb2->first_order_sources_integrand[index_type],
                  pgb2->tau_size_bessel * sizeof(double *),
                  ppt->error_message);
      /* Allocate memory for pgb2->first_order_sources_integrand[index_type] */
      for (int index_tau = 0; index_tau < pgb2->tau_size_bessel; index_tau++) {
          /* Loop over type and tau. For each of them, allocate memory
           for pgb2->first_order_sources_integrand[index_type][index_tau][index_k_bessel]  */
         class_alloc(pgb2->first_order_sources_integrand[index_type][index_tau],
                      pgb2->k_size_bessel * sizeof(double),
                      ppt->error_message);
    }
  }

  /* Define an array of values of first order transfer functions with integrals (lensing etc.), these source terms have an extra
    index (index_l):
              pgb2->first_order_sources_integ[index_type][index_l][index_tau][index_k_bessel] */


  class_alloc(pgb2->first_order_sources_integ, pgb2->type_size * sizeof(double ***), ppt->error_message);
    for (int index_type = 0; index_type < pgb2->type_size; index_type++) {

      /* Allocate memory for pgb2->first_order_sources_integ[index_type][index_tau] */
      class_alloc(pgb2->first_order_sources_integ[index_type],
                  ptr->l_size[ppt->index_md_scalars] * sizeof(double **),
                  ppt->error_message);

      for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {

        class_alloc(pgb2->first_order_sources_integ[index_type][index_l],
                    pgb2->tau_size_selection * sizeof(double *),
                    ppt->error_message);
      /* Allocate memory for pgb2->first_order_sources_integ[index_type] */
        for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {

            /* Loop over type, l and tau. For each of them, allocate memory
             for pgb2->first_order_sources_integ[index_type][index_l][index_tau][index_k_bessel]  */
           class_alloc(pgb2->first_order_sources_integ[index_type][index_l][index_tau],
                        pgb2->k_size_bessel * sizeof(double),
                        ppt->error_message);
        }
      }
    }



  //printf("%g\n",ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][100 * ppt->k_size[ppt->index_md_scalars] + 20]);


  /***************************************************************************
  ==============   Interpolate on first_order_sources array =================
  ****************************************************************************/

/* Now fill and interpolate this new pointer-array pgb2->first_order_sources[index_type][index_tau][index_k]
  with information from the pre computed ppt->quadsources[index_md][index_ic*ppt->tp_size[index_md]+index_type]
  [index_tau * ppt->k_size[index_md] + index_k]. */



  printf("Source array allocation completed.\n" );
  // set up sources for delta_cdm

  /* NOTE: The following first_order_sources interpolation is only necressary for each of the perturbations delta, v, phi etc. ,
    terms in the galaxy number over-density which share these terms, can just reuse the associated */

  double intermediate;
  double intermediate_plus;
  //for (int index_type = 0; index_type < ppt->tp_size[index_md]; index_type++) {
    for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
      //last_index_k = 0;
      index_k = 0;

      for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
        index = 0;

        double tau = pgb2->tau_sampling_cls[index_tau];

        class_call(index_of_tau_sampling_quadsources(tau, &index, ppt),ppt->error_message,pgb2->error_message);

        double k = pgb2->k_bessel[index_k_bessel];
        //printf("index_k = %d\n", index_k);
        class_call(index_of_k(k, &index_k, ppt), ppt->error_message, pgb2->error_message);
        //class_call(index_of_k(k, &index_k, &last_index_k, ppt), ppt->error_message, pgb2->error_message);
        //last_index_k = index_k;


        f = (pgb2->tau_sampling_cls[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

        intermediate  = (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k]+
            (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k]);

        intermediate_plus =  (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k+1]+
            (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k+1]);



        g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);



        pgb2->first_order_sources[pgb2->index_type_delta_cdm][index_tau][index_k_bessel] = g*intermediate_plus +(1-g)*intermediate;
        //printf("pgb2->first_order_sources[%d][%d][%d] = %g\n", pgb2->index_type_delta_cdm, index_tau, index_k_bessel);
      }
    }

    printf("Reach here1\n" );

    // set up sources for rsd


    intermediate = 0;
    intermediate_plus = 0;
    int last_index_rsd;
    int last_index_k_rsd;
    int index_k_rsd;


  //for (int index_type = 0; index_type < ppt->tp_size[index_md]; index_type++) {
    for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
      // NOTE: is ppt->index_qs_v_cdm velocity or velocity divergence?
      index_k_rsd = 0;


      for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
        index = 0;
        //last_index_k_rsd = 0;

        double tau = pgb2->tau_sampling_cls[index_tau];

        class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);

        double k = pgb2->k_bessel[index_k_bessel];

        class_call(index_of_k(k, &index_k_rsd, ppt), pgb2->error_message, pgb2->error_message);
        //class_call(index_of_k(k, &index_k_rsd, &last_index_k_rsd, ppt), pgb2->error_message, pgb2->error_message);
        //last_index_k = index_k;


        f = (pgb2->tau_sampling_cls[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

        intermediate  = (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k]+
            (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k]);

        intermediate_plus =  (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k+1]+
            (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k+1]);



        g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);



        pgb2->first_order_sources[pgb2->index_type_rsd][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);

      }
    }

      // set up sources for lensing


    intermediate = 0;
    intermediate_plus = 0;
    int last_index_lens;
    int last_index_k_lens;
    int index_k_lens;
    double ** phi_plus_psi;

    /* Allocate a temporary array phi_plus_psi[index_tau][index_k] that stores the quadsources transfer function for phi+psi
        purely for brevity*/

    class_alloc(phi_plus_psi, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
      for (int index = 0; index < ppt->tau_size_quadsources; index++) {
        /* Allocate memory for pgb2->first_order_sources[index_type][index_tau] */
        class_alloc(phi_plus_psi[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);
    }

    for (int index = 0; index < ppt->tau_size_quadsources; index++) {
      for (int index_k = 0; index_k < ppt->k_size[ppt->index_md_scalars]; index_k++) {
        phi_plus_psi[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_phi][index * ppt->k_size[ppt->index_md_scalars] + index_k]
          + ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_psi][index * ppt->k_size[ppt->index_md_scalars] + index_k];

      //  printf("phi_plus_psi[%d][%d] = %g\n",index, index_k, phi_plus_psi[index][index_k]);
      }
    }
    //herehere

    printf("starting pgb2->index_type_lens\n" );
    for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){

        last_index_k_lens = 0;
        index_k_lens = 0;
        int index_tau_source;
        int last_index;
        int index_tau_lens;
        double r_lens;
        double r = pba->conformal_age - pgb2->tau_sampling_cls[index_tau];
        double lensing_result;
        double x;

        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
          index = 0;

          //printf("index_tau = %d, index_k_bessel = %d\n", index_tau, index_k_bessel);
          //last_index_k_rsd = 0;

          double tau = pgb2->tau_sampling_cls[index_tau];

          class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);

          double k = pgb2->k_bessel[index_k_bessel];

          class_call(index_of_k(k, &index_k_lens, ppt), pgb2->error_message, pgb2->error_message);
          //class_call(index_of_k(k, &index_k_lens, &last_index_k_lens, ppt), pgb2->error_message, pgb2->error_message);
          //last_index_k = index_k;

          f = (pgb2->tau_sampling_cls[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

          intermediate  = f*phi_plus_psi[index+1][index_k]+(1-f)*phi_plus_psi[index][index_k];

          intermediate_plus =  f*phi_plus_psi[index+1][index_k+1]+(1-f)*phi_plus_psi[index][index_k+1];

          g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

          pgb2->first_order_sources[pgb2->index_type_phi_plus_psi][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);

          //printf("pgb2->first_order_sources[pgb2->index_phi_plus_psi][%d][%d] = %g\n", index_tau, index_k_bessel,pgb2->first_order_sources[pgb2->index_type_phi_plus_psi][index_tau][index_k_bessel]);

        }
    }

    for (int index_tau = 0; index_tau < pgb2->tau_size_bessel; index_tau++){

        last_index_k = 0;
        index_k_lens = 0;
        int last_index;
        int index_tau_lens;
        double lensing_result;
        f = 0.;
        g = 0;

        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
          index = 0;

          double tau = pgb2->tau_sampling_bessel[index_tau];

          class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);

          double k = pgb2->k_bessel[index_k_bessel];

          //class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);
          class_call(index_of_k_old(k,
                         &index_k,
                         &last_index_k,
                         ppt),
                         pgb2->error_message,
                         pgb2->error_message);
          //last_index_k = index_k;

          f = (pgb2->tau_sampling_bessel[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

          intermediate  = f*phi_plus_psi[index+1][index_k]+(1-f)*phi_plus_psi[index][index_k];

          intermediate_plus =  f*phi_plus_psi[index+1][index_k+1]+(1-f)*phi_plus_psi[index][index_k+1];

          g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);



          /*pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau][index_k_bessel] = intermediate_plus*(pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k])
             + intermediate*(1-(pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]));*/
          pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);
          /*if (index_k_bessel == 0) {
            printf("tau_bessel = %g, tau_qs = %g\n", tau, ppt->tau_sampling_quadsources[index]);
          }
          if (index_tau == pgb2->tau_size_bessel -2) {
            printf("%g    %g    %g    %g\n", k, pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau][index_k_bessel], ppt->k[ppt->index_md_scalars][index_k], phi_plus_psi[index][index_k]);
          }*/



        //if (index_tau == 4900 ) {
            //printf("pgb2->first_order_sources_integrand[pgb2->index_phi_plus_psi][%d][%d] = %g\n", index_tau, index_k_bessel, pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau][index_k_bessel]);
            //printf("f = %g\n",f);
            //printf("g = %g\n",g);
            //printf("intermediate = %g\n",intermediate );
            //printf("intermediate_plus = %g\n",intermediate_plus );
          //}
          //printf(" \n");
          //printf("phi_plus_psi[%d][%d] = *%g*, k = %g, tau_qs = %g\n", index, index_k, phi_plus_psi[index][index_k],ppt->k[ppt->index_md_scalars][index_k],ppt->tau_sampling_quadsources[index]);
          //printf("pgb2->first_order_sources_integrand[pgb2->index_phi_plus_psi][%d][%d] = *%g*, k = %g, tau_bessel = %g\n", index_tau, index_k_bessel,pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau][index_k_bessel], pgb2->k_bessel[index_k_bessel], pgb2->tau_sampling_bessel[index_tau]);
          //printf("phi_plus_psi[%d][%d] = *%g*, k =%g, tau_qs = %g\n", index, index_k+1, phi_plus_psi[index][index_k+1], ppt->k[ppt->index_md_scalars][index_k+1],ppt->tau_sampling_quadsources[index]);
        }
    }
    int index_kg = 5000;
    int index_kg_ppt;
    double kg = pgb2->k_bessel[index_kg];

    class_call(index_of_k_old(kg,
                   &index_kg_ppt,
                   &last_index_k,
                   ppt),
                   pgb2->error_message,
                   pgb2->error_message);

    /*for (int index_tau = 0; index_tau < pgb2->tau_size_bessel; index_tau++){

      tau = pgb2->tau_sampling_bessel[index_tau];

      class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);


      printf("%s\n", printf("%g    %g    %g    %g\n", tau, pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau][index_kg], ppt->tau_sampling_quadsources[index], phi_plus_psi[index][index_kg_ppt]));
    }

    exit(0);*/

    int index_k_f = 5000;
    double k_f;
    index = 0;
    last_index_k = 0;
    // INTERPOLATION INTERPOLATE
    for (int index_tau = 0; index_tau < pgb2->tau_size_bessel; index_tau++) {
      last_index_k_lens = 0;
      index_k_lens = 0;
      int last_index;
      int index_tau_lens;

      double tau = pgb2->tau_sampling_bessel[index_tau];

      class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);
      //class_call(index_of_k(k_f, &index_k, ppt), pgb2->error_message, pgb2->error_message);
      k_f = pgb2->k_bessel[index_k_f];

      //class_call(index_of_k(k_f, &index_k, ppt), pgb2->error_message, pgb2->error_message);
      class_call(index_of_k_old(k_f,
                     &index_k,
                     &last_index_k,
                     ppt),
                     pgb2->error_message,
                     pgb2->error_message);
      //last_index_k = index_k;

      /*pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau][index_k_f] = intermediate_plus*(pgb2->k_bessel[index_k_f]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k])
         + intermediate*(1-(pgb2->k_bessel[index_k_f]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]));*/
      pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);
      //printf(" \n");
      //printf("k_qs = %g\n", ppt->k[ppt->index_md_scalars][index_k]);
      //printf("k_bessel = %g\n", pgb2->k_bessel[index_k_f] );
      //printf("k_qs+1 = %g\n", ppt->k[ppt->index_md_scalars][index_k+1]);
      //printf("phi_plus_psi[%d][%d] = *%g*, k = %g, tau_qs = %g\n", index, index_k, phi_plus_psi[index][index_k],ppt->k[ppt->index_md_scalars][index_k],ppt->tau_sampling_quadsources[index]);
      //printf("pgb2->first_order_sources_integrand[pgb2->index_phi_plus_psi][%d][%d] = *%g*, k = %g, tau_bessel = %g\n", index_tau, index_k_f,pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau][index_k_f], pgb2->k_bessel[index_k_f], pgb2->tau_sampling_bessel[index_tau]);
      //printf("phi_plus_psi[%d][%d] = *%g*, k =%g, tau_qs = %g\n", index+1, index_k, phi_plus_psi[index+1][index_k], ppt->k[ppt->index_md_scalars][index_k],ppt->tau_sampling_quadsources[index+1]);


    }




    double *** lensing_test_array;

    class_alloc(lensing_test_array,
                ptr->l_size[ppt->index_md_scalars] * sizeof(double**),
                pgb2->error_message);
    for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      class_alloc(lensing_test_array[index_l],
                  pgb2->tau_size_selection * sizeof(double*),
                  pgb2->error_message);
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
        class_alloc(lensing_test_array[index_l][index_tau],
                    pgb2->k_size_bessel * sizeof(double),
                    pgb2->error_message);
      }
    }

    /*int index_k_bessel2 = 4000;
    for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
      printf(" \n");
      printf("ORIGINAL[ppgb2->index_type_phi_plus_psi][tau[%d] = %g][k[%d] = %g] = %g\n", index_tau, pgb2->tau_sampling_cls[index_tau], index_k_bessel2, pgb2->k_bessel[index_k_bessel2], pgb2->k_bessel[index_k_bessel2], pgb2->first_order_sources[pgb2->index_type_phi_plus_psi][index_tau][index_k_bessel2]);
      for (int index_tau_lens = (index_tau * bessel_boost); index_tau_lens < pgb2->tau_size_bessel-1; index_tau_lens++) {

        if (index_tau_lens == index_tau * bessel_boost){
          printf("NEW[pgb2->index_type_phi_plus_psi][tau_lens[%d] = %g][k[%d] = %g] = %g\n",index_tau_lens, pgb2->tau_sampling_bessel[index_tau_lens], index_k_bessel2, pgb2->k_bessel[index_k_bessel2], pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau_lens][index_k_bessel2]);
        }

        if (index_tau_lens == ceil(((index_tau * bessel_boost)-pgb2->tau_size_bessel+2)/2)) {
          printf("NEW[pgb2->index_type_phi_plus_psi][tau_lens[%d] = %g][k[%d] = %g] = %g\n",index_tau_lens, pgb2->tau_sampling_bessel[index_tau_lens], index_k_bessel2, pgb2->k_bessel[index_k_bessel2], pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau_lens][index_k_bessel2]);
        }

        if (index_tau_lens == (pgb2->tau_size_bessel-2)) {
          printf("NEW[pgb2->index_type_phi_plus_psi][tau_lens[%d] = %g][k[%d] = %g] = %g\n",index_tau_lens, pgb2->tau_sampling_bessel[index_tau_lens], index_k_bessel2, pgb2->k_bessel[index_k_bessel2], pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau_lens][index_k_bessel2]);
        }
      }
    }
    exit(0);*/

    /* This is the ratio between the two equally spaced time grids. */
    // NOTE: ALERT l has been changed to run up to 4.
    for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){

      last_index_k_lens = 0;
      index_k_lens = 0;
      int index_tau_cls;
      int last_index;
      int index_tau_lens;
      int first_index_tau_in_lens;
      double r_lens;
      double r = pba->conformal_age - pgb2->tau_sampling_cls[index_tau];
      double lensing_result;
      double f_minus, f_plus, w_minus, w_plus, f_interpolate, tau_lens;
      double weight, end_weights;
      double x;
      double j,j1;

      for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

          lensing_result = 0.;
          first_index_tau_in_lens = index_tau * bessel_boost;
          //printf("pgb2->first_order_sources[pgb2->index_type_phi_plus_psi][%d][%d] = %g\n",index_tau, index_k_bessel,pgb2->first_order_sources[pgb2->index_type_phi_plus_psi][index_tau][index_k_bessel]);

          for (index_tau_lens = (index_tau * bessel_boost); index_tau_lens < pgb2->tau_size_bessel-1; index_tau_lens++) {
            // NOTE that we skip the last index "pgb2->tau_size_bessel-1" to avoid the singularity. This would analytically drop out anyways
            // for l not 1.

            //printf("index_tau * bessel_boost = %d . %d\n", index_tau * bessel_boost, index_tau_lens );
            tau_lens = pgb2->tau_sampling_bessel[index_tau_lens];


            if(x>pbs->x_max ){
            printf("ALERT! x= %g \n",x);continue;}



            class_call(bessel_at_x(pbs,pgb2->k_bessel[index_k_bessel]*(pba->conformal_age -tau_lens), index_l, &j), pbs->error_message, pgb2->error_message);

            index_of_tau_sampling_cls(tau_lens, &index_tau_cls, &last_index, pgb2);

            first_index_tau_in_lens = index_tau * bessel_boost;

            /* If we are taking the first or last trapezoidal weight then take the first weight, else take the other weight */
            if (index_tau_lens == first_index_tau_in_lens || index_tau_lens == (pgb2->tau_size_bessel-1) ){
              weight = (pgb2->tau_sampling_bessel[first_index_tau_in_lens]-pgb2->tau_sampling_bessel[pgb2->tau_size_bessel-1])/(2.0*((pgb2->tau_size_bessel-1)-first_index_tau_in_lens));
            }

            else{
              weight = (pgb2->tau_sampling_bessel[first_index_tau_in_lens]-pgb2->tau_sampling_bessel[pgb2->tau_size_bessel-1])/((pgb2->tau_size_bessel-1)-first_index_tau_in_lens);
            }

            lensing_result += ((pba->conformal_age - pgb2->tau_sampling_bessel[index_tau_lens])-(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]))
              * j * weight * pgb2->first_order_sources_integrand[pgb2->index_type_phi_plus_psi][index_tau_lens][index_k_bessel]/((pba->conformal_age - pgb2->tau_sampling_bessel[index_tau_lens])*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]));
          }

          pgb2->first_order_sources_integ[pgb2->index_type_lens][index_l][index_tau][index_k_bessel] = lensing_result;
        }
      }
    }

  printf("First order source array filled.\n");

  /***************************************************************************
  ==========================   k-Integration =================================
  ****************************************************************************/

  printf("Starting k-integration\n");

  /* Allocate array for Cl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] */
  class_alloc(pgb2->Cl, pgb2->type_size * sizeof(double ****), pgb2->error_message);
  for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    class_alloc(pgb2->Cl[index_type_first], pgb2->type_size * sizeof(double ***), pgb2->error_message);
    for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
      class_alloc(pgb2->Cl[index_type_first][index_type_second], ptr->l_size[ppt->index_md_scalars] * sizeof(double **), pgb2->error_message);
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
        class_alloc(pgb2->Cl[index_type_first][index_type_second][index_l], pgb2->tau_size_selection * sizeof(double *), pgb2->error_message);
        for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++) {
          class_alloc(pgb2->Cl[index_type_first][index_type_second][index_l][index_tau_first], pgb2->tau_size_selection * sizeof(double), pgb2->error_message);
        }
      }
    }
  }
  printf("Allocating size %ix%ix%ix%ix%i bytes \n", pgb2->type_size, pgb2->type_size, ptr->l_size[ppt->index_md_scalars], pgb2->tau_size_selection, pgb2->tau_size_selection);

  /* Allocate array for Cl_final[index_type_first][index_type_second][index_l][bin_first][bin_second] */
  class_alloc(pgb2->Cl_final, pgb2->type_size * sizeof(double ****), pgb2->error_message);
  for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    class_alloc(pgb2->Cl_final[index_type_first], pgb2->type_size * sizeof(double ***), pgb2->error_message);
    for (int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
      class_alloc(pgb2->Cl_final[index_type_first][index_type_second], ptr->l_size[ppt->index_md_scalars] * sizeof(double **), ppt->error_message);
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
        class_alloc(pgb2->Cl_final[index_type_first][index_type_second][index_l], ppt->selection_num * sizeof(double *), ppt->error_message);
        for (int bin_first = 0; bin_first < ppt->selection_num; bin_first++) {
          class_alloc(pgb2->Cl_final[index_type_first][index_type_second][index_l][bin_first], ppt->selection_num * sizeof(double), ppt->error_message);
        }
      }
    }
  }


/* Allocate and fill the trapezoidal weights for the k-integration */
  class_alloc(pgb2->w_trapz_k,
              pgb2->k_size_bessel * sizeof(double),
              ppt->error_message);


  class_call(array_trapezoidal_weights(pgb2->k_bessel,
                                       pgb2->k_size_bessel,
                                       pgb2->w_trapz_k,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);



  double integ, integ_dens_rsd;
  double result2;


  printf("integrating k between %g and %g\n",ppt->k[ppt->index_md_scalars][0],ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1] );

  printf("pgb2->type_size = %d\n", pgb2->type_size);



  double p1,p2,f1,f2,j1,j2;
  int k_test_index;


  /* Fill the array pgb2->Cl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] by calling the integral() function */
  for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    for (int k = 0; k < index_type_first+1; k++){
      printf("X");
    }
    printf("\n");

    for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
        for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++){
          for(int index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++){

            class_call(integral(pba, ppm, pbs, pgb2, ppt, ptr,
                                index_type_first,
                                index_type_second,
                                index_tau_first,
                                index_tau_second,
                                index_l,
                                pvecback,
                                &integ),
                       pgb2->error_message,
                       pgb2->error_message);

            pgb2->Cl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] = integ;

            if (integ != integ) {
              printf("pgb2->Cl[%d][%d][%d][%d][%d] = %g\n", index_type_first, index_type_second, index_l, index_tau_first, index_tau_second, pgb2->Cl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second]);
            }
          }
        }
      }
    }
  }

  printf("Done k-integration.\n");








  /***************************************************************************
  ==========================   time-Integration ==============================
  ****************************************************************************/

  /* Allocate the array pgb2->integral_over_single_window[index_type_first][index_type_second][index_l][bin1][bin2][index_tau_second].
    Technically this quantity does not depend on bin2 (third argument), however it is necessary for the next integration
    over index_tau_sampling ,to give the cl_sampling_selection context of the index_tau_second. */

  printf("Allocating size %ix%ix%ix%ix%ix%i bytes \n", pgb2->type_size, pgb2->type_size, ptr->l_size[ppt->index_md_scalars], ppt->selection_num, ppt->selection_num,pgb2->tau_size_selection);
  class_alloc(pgb2->integral_over_single_window,
              pgb2->type_size * sizeof(double*****),
              pgb2->error_message);

  for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    class_alloc(pgb2->integral_over_single_window[index_type_first],
                pgb2->type_size * sizeof(double****),
                pgb2->error_message);

    for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
      class_alloc(pgb2->integral_over_single_window[index_type_first][index_type_second],
                  ptr->l_size[ppt->index_md_scalars] * sizeof(double***),
                  pgb2->error_message);

      for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
        class_alloc(pgb2->integral_over_single_window[index_type_first][index_type_second][index_l],
                    ppt->selection_num * sizeof(double**),
                    pgb2->error_message);

        for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
          class_alloc(pgb2->integral_over_single_window[index_type_first][index_type_second][index_l][bin1],
                      ppt->selection_num * sizeof(double*),
                      pgb2->error_message);

          for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
            class_alloc(pgb2->integral_over_single_window[index_type_first][index_type_second][index_l][bin1][bin2],
                        pgb2->tau_size_selection * sizeof(double),
                        pgb2->error_message);

          }
        }
      }
    }
  }


  printf("Starting time-integration.\n");


  int  bin1 = 0;
  int  bin2 = 0;
  int  bin3 = 0;
  double tmp2 = 0;
  int index_tau1;
  int index_tau2;
  int last_index1;
  int last_index2;
  double temp, temp_minus, temp_plus;

  for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
      for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
        for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
          for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){

            pgb2->Cl_final[index_type_first][index_type_second][index_l][bin1][bin2] = 0.;


    /* Final result Cl_final[index_type_first][index_type_second][index_l][bin1][bin2] is defined on the pgb2->tau_sampling_cls, so an interpolation is required
          from Cl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] which is defined on pgb2->tau_size_selection[bin][index_tau]*/
            for(index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++){


              double temp123 = 0.;
              double tau2 = pgb2->tau_sampling_selection[bin2][index_tau_second];
              ;

              for(index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++){
                //printf("indexing %dx%dx%dx%dx%dx%dx%d\n",index_type_first, index_type_second, bin1, bin2, index_l,index_tau_second,index_tau_first);

                double tau1 = pgb2->tau_sampling_selection[bin1][index_tau_first];

                index_of_tau_sampling_cls(tau1, &index_tau1, &last_index1, pgb2);
                index_of_tau_sampling_cls(tau2, &index_tau2, &last_index2, pgb2);


                temp_minus = pgb2->Cl[index_type_first][index_type_second][index_l][index_tau1-1][index_tau2-1]*(pgb2->tau_sampling_cls[index_tau1]-tau1)
                              + pgb2->Cl[index_type_first][index_type_second][index_l][index_tau1][index_tau2-1]*(tau1-pgb2->tau_sampling_cls[index_tau1-1]);
                temp_minus /= (pgb2->tau_sampling_cls[index_tau1] - pgb2->tau_sampling_cls[index_tau1-1]);

                temp_plus = pgb2->Cl[index_type_first][index_type_second][index_l][index_tau1-1][index_tau2]*(pgb2->tau_sampling_cls[index_tau1]-tau1)
                              + pgb2->Cl[index_type_first][index_type_second][index_l][index_tau1][index_tau2]*(tau1-pgb2->tau_sampling_cls[index_tau1-1]);
                temp_plus /= (pgb2->tau_sampling_cls[index_tau1] - pgb2->tau_sampling_cls[index_tau1-1]);

                temp = temp_minus * (pgb2->tau_sampling_cls[index_tau2]-tau2)
                              + temp_plus * (tau2-pgb2->tau_sampling_cls[index_tau2-1]);

                temp /= (pgb2->tau_sampling_cls[index_tau2] - pgb2->tau_sampling_cls[index_tau2-1]);

                temp123 += temp * w_trapz[bin1][index_tau_first] * selection[bin1][index_tau_first];

              }


              pgb2->Cl_final[index_type_first][index_type_second][index_l][bin1][bin2] += temp123 * w_trapz[bin2][index_tau_second]
                  * selection[bin2][index_tau_second];

              pgb2->integral_over_single_window[index_type_first][index_type_second][index_l][bin1][bin2][index_tau_second] = temp123;
            }
          }
        }
      }
    }
  }

printf("integrating between k =%g and %g\n",pgb2->k_bessel[index_k_bessel],pgb2->k_bessel[pgb2->k_size_bessel-1]);
for  (index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
    //printf("l*(l+1)*Cl_final[%d][0][0] = %g\n",ptr->l[index_l], pgb2->Cl_final[index_type_first][index_type_second][index_l][0][0]*ptr->l[index_l]*(ptr->l[index_l]+1.));
    printf("%d    %g\n",ptr->l[index_l], pgb2->Cl_final[pgb2->index_type_lens][pgb2->index_type_lens][index_l][0][0]*ptr->l[index_l]*(ptr->l[index_l]+1.)/(2*_PI_));
  }
//exit(0);


  printf("Starting bispectrum computation.\n");


/***************************************************************************
=====================   Bispectrum Computation =============================
****************************************************************************/

  int index_l1,index_l2,index_l3;

//NOTE: Should each dimension in the asym_redgalbispectrum array be allocated equally?
/* Allocation for the array pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l2][index_l3][bin1][bin2][bin3]. This quantity is then summed
  with different permutations to yield the full reduced galaxy bispectrum.*/
  class_alloc(pgb2->asym_redgalbispectrum,
              pgb2->type_size * sizeof(double ********),
              pgb2->error_message);

  for (int index_type_SO1 = 0; index_type_SO1 < pgb2->type_size; index_type_SO1++){
    class_alloc(pgb2->asym_redgalbispectrum[index_type_SO1],
                pgb2->type_size * sizeof(double *******),
                pgb2->error_message);

    for (int index_type_SO2 = 0; index_type_SO2 < pgb2->type_size; index_type_SO2++){
      class_alloc(pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2],
                  pgb2->type_size * sizeof(double ******),
                  pgb2->error_message);

      for (int index_type_FO1 = 0; index_type_FO1 < pgb2->type_size; index_type_FO1++){
        class_alloc(pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1],
                    pgb2->type_size * sizeof(double *****),
                    pgb2->error_message);

        for (int index_type_FO2 = 0; index_type_FO2 < pgb2->type_size; index_type_FO2++){
          class_alloc(pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2],
                      ptr->l_size[ppt->index_md_scalars] * sizeof(double ****),
                      pgb2->error_message);

          for(int index_l1 = 0; index_l1 < ptr->l_size[ppt->index_md_scalars]; index_l1++){

            class_alloc(pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1],
                        ptr->l_size[ppt->index_md_scalars] * sizeof(double ***),
                        pgb2->error_message);

            for (int index_l2 = 0; index_l2 < ptr->l_size[ppt->index_md_scalars]; index_l2++) {

              class_alloc(pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2],
                          ppt->selection_num * sizeof(double**),
                          pgb2->error_message);

              for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {

              class_alloc(pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][bin1],
                          ppt->selection_num * sizeof(double*),
                          pgb2->error_message);

                for ( int bin2 = 0; bin2 < ppt->selection_num; bin2++) {

                  class_alloc(pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][bin1][bin2],
                              ppt->selection_num * sizeof(double),
                              pgb2->error_message);
                }
              }
            }
          }
        }
      }
    }
  }


/* Allocation for the array pgb2->redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][index_l3][bin1][bin2][bin3]
  This is the reduced galaxy bispectrum that sums the perumations of pgb2->redgalbispectrum */

  printf("Allocating size %ix%ix%ix%ix%ix%ix%ix%ix%ix%i bytes \n", pgb2->type_size, pgb2->type_size, pgb2->type_size, ptr->l_size[ppt->index_md_scalars], ptr->l_size[ppt->index_md_scalars], ptr->l_size[ppt->index_md_scalars], ppt->selection_num, ppt->selection_num, ppt->selection_num);

  class_alloc(pgb2->redgalbispectrum,
              pgb2->type_size * sizeof(double *********),
              pgb2->error_message);

  for (int index_type_SO1 = 0; index_type_SO1 < pgb2->type_size; index_type_SO1++){
    class_alloc(pgb2->redgalbispectrum[index_type_SO1],
                pgb2->type_size * sizeof(double ********),
                pgb2->error_message);

    for (int index_type_SO2 = 0; index_type_SO2 < pgb2->type_size; index_type_SO2++){
      class_alloc(pgb2->redgalbispectrum[index_type_SO1][index_type_SO2],
                  pgb2->type_size * sizeof(double *******),
                  pgb2->error_message);

      for (int index_type_FO1 = 0; index_type_FO1 < pgb2->type_size; index_type_FO1++){
        class_alloc(pgb2->redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1],
                    pgb2->type_size * sizeof(double ******),
                    pgb2->error_message);

        for (int index_type_FO2 = 0; index_type_FO2 < pgb2->type_size; index_type_FO2++){
          class_alloc(pgb2->redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2],
                      ptr->l_size[ppt->index_md_scalars] * sizeof(double *****),
                      pgb2->error_message);

          for(int index_l1 = 0; index_l1 < ptr->l_size[ppt->index_md_scalars]; index_l1++){
            class_alloc(pgb2->redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1],
                        ptr->l_size[ppt->index_md_scalars] * sizeof(double ****),
                        pgb2->error_message);

            for (int index_l2 = 0; index_l2 < ptr->l_size[ppt->index_md_scalars]; index_l2++) {
              class_alloc(pgb2->redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2],
                          ptr->l_size[ppt->index_md_scalars] * sizeof(double***),
                          pgb2->error_message);

              for (int index_l3 = 0; index_l3 < ptr->l_size[ppt->index_md_scalars]; index_l3++) {
                class_alloc(pgb2->redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][index_l3],
                            ppt->selection_num * sizeof(double**),
                            pgb2->error_message);

                for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
                  class_alloc(pgb2->redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][index_l3][bin1],
                            ppt->selection_num * sizeof(double*),
                            pgb2->error_message);

                  for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
                    class_alloc(pgb2->redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][index_l3][bin1][bin2],
                                ppt->selection_num * sizeof(double),
                                pgb2->error_message);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  printf("Allocated size %ix%ix%ix%ix%ix%ix%ix%ix%ix%i bytes \n", pgb2->type_size, pgb2->type_size, pgb2->type_size, ptr->l_size[ppt->index_md_scalars], ptr->l_size[ppt->index_md_scalars], ptr->l_size[ppt->index_md_scalars], ppt->selection_num, ppt->selection_num, ppt->selection_num);

/* It is now necessary to introduce new indices <A^(1)B^(1)C^(1)D^(1)>, let A and B form the second order term, then
    index_type_SO1 and index_type_SO2 denote types A and B, then index_type_FO1 and index_type_FO2 denote the types two
     first-order terms C and D. */

/* NOTE: Need to loop correctly over each type such that the pairs of types that feature in integral_over_single_window
  are valid, terms A and B which together make up the second order term cannot be paired together. */


  for(int index_type_SO1 = 0; index_type_SO1 < pgb2->type_size; index_type_SO1++){
    for(int index_type_SO2 = 0; index_type_SO2 < pgb2->type_size; index_type_SO2++){
      for(int index_type_FO1 = 0; index_type_FO1 < pgb2->type_size; index_type_FO1++){
        for(int index_type_FO2 = 0; index_type_FO2 < pgb2->type_size; index_type_FO2++){
          for(int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
            for(int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
              for(int bin3 = 0; bin3 < bin2+1; bin3++) {
                for(int index_l2 = 0; index_l2 < ptr->l_size[ppt->index_md_scalars]; index_l2++){
                  for(int index_l3 = 0; index_l3 < index_l2+1; index_l3++){
                    pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l2][index_l3][bin1][bin2][bin3] = 0.;
  /* The bin on the left of pgb2->asym_redgalbispectrum is taken to be the outer integration (the one corresponding to the second order part).
  The outer integration time-variable is the one that corresponds to the second-order
    perturbaton. Since this choice is arbitary, we must permute the bins. */

                    for(int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
                      pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l2][index_l3][bin1][bin2][bin3] +=
                        pgb2->integral_over_single_window[index_type_SO1][index_type_FO1][index_l2][bin2][bin1][index_tau]
                        * pgb2->integral_over_single_window[index_type_SO2][index_type_FO2][index_l3][bin3][bin1][index_tau] * w_trapz[bin1][index_tau]
                            * selection[bin1][index_tau];

                      //printf("pgb2->integral_over_single_window[%d][%d][%d][%d] = %g\n",index_l2, bin1, bin2, index_tau, pgb2->integral_over_single_window[index_l2][bin1][bin2][index_tau]);
                      //printf("selection[%d][%d] = %g, w_trapz[%d][%d] = %g\n",bin2, index_tau,selection[bin2][index_tau], bin2, index_tau, w_trapz[bin2][index_tau]);
                    }
                    printf("pgb2->asym_redgalbispectrum[%d][%d][%d][%d][%d][%d][%d][%d][%d] = %g\n", index_type_SO1, index_type_SO2, index_type_FO1, index_type_FO2, index_l2, index_l3, bin1, bin2, bin3);
                    //printf("pgb2->asym_redgalbispectrum[%d][%d][%d][%d][%d] = %g\n", index_l2, index_l3, bin1, bin2, bin3,pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l2][index_l3][bin1][bin2][bin3]);

                  /*  printf("pgb2->asym_redgalbispectrum[%d][%d][%d] = %g\n",
                            index_l1,
                            index_l2,
                            index_l3,
                            pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][index_l3]); */


                  }
                }
              }
            }
          }
        }
      }
    }
  }


  for(int index_type_SO1 = 0; index_type_SO1 < pgb2->type_size; index_type_SO1++){
    for(int index_type_SO2 = 0; index_type_SO2 < index_type_SO1; index_type_SO2++){
      for(int index_type_FO1 = 0; index_type_FO1 < pgb2->type_size; index_type_FO1++){
        for(int index_type_FO2 = 0; index_type_FO2 < pgb2->type_size; index_type_FO2++){
          for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
            for (int bin2 = 0; bin2 < bin1+1; bin2++) {
              for (int bin3 = 0; bin3 < bin2+1; bin3++) {
                for(int index_l1 = 0; index_l1 < ptr->l_size[ppt->index_md_scalars]; index_l1++){
                  for(int index_l2 = 0; index_l2 < index_l1+1; index_l2++){

        // NOTE: the relevance of the left bin and the right bin need to be worked out, atm it is redundant because of symmetry, but in general this is not so.
                    for (int index_l3 = 0; index_l3 < index_l2+1; index_l3++) {

                    /*  printf("[%d][%d][%d][%d][%d][%d][%d][%d][%d][%d]\n",
                              index_type_SO1,
                              index_type_SO2,
                              index_type_FO1,
                              index_type_FO2,
                              index_l1,
                              index_l2,
                              index_l3,
                              bin1,
                              bin2,
                              bin3);*/

                      pgb2->redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][index_l3][bin1][bin2][bin3] =
                        pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l2][index_l3][bin1][bin2][bin3]
                          + pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l3][bin2][bin1][bin3]
                            + pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][bin3][bin1][bin2];

                        /*  double result  =  pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l2][index_l3][bin1][bin2][bin3]
                                + pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l3][bin2][bin1][bin3]
                                  + pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][bin3][bin1][bin2]; */

                          //printf("result = %g\n",result );
                      printf("pgb2->redgalbispectrum[%d][%d][%d][%d][%d][%d][%d][%d][%d][%d] = %g\n",
                              index_type_SO1,
                              index_type_SO2,
                              index_type_FO1,
                              index_type_FO2,
                              index_l1,
                              index_l2,
                              index_l3,
                              bin1,
                              bin2,
                              bin3,
                              pgb2->redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l1][index_l2][index_l3][bin1][bin2][bin3]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

/*  for (int index_l1 = 0; index_l1 < ptr->l_size[ppt->index_md_scalars]; index_l1++) {
    for (int index_l2 = 0; index_l2 < index_l1+1; index_l2++) {
      for (int index_l3 = 0; index_l3 < index_l2+1; index_l3++) {
        printf("pgb2->redgalbispectrum[%d][%d][%d][0][0][0] = %g\n", index_l1,index_l2,index_l3);
      }
    }
  }
*/


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

  } // end of(ppr2->store_sources_to_disk)



free(pgb2->tau_sampling_cls);
printf("End of galbispectra2!\n");
  return _SUCCESS_;

} // end of galbispectra2_init()
