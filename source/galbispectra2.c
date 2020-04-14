/* Module to compute galaxy bispectra */

#include <math.h>
#include "galbispectra2.h"
#include "perturbations2.h"

double bessel_at_x_first_deriv_old(struct galbispectra2 * pgb2,
                       struct bessels * pbs,
                       double z,
                       int index_l,
                       double * result){

    double out;
    double j1;
    double j2;

    if (z == 0.0) {
      *result = 0.0;
    }

    else{
     /*using recurrence relation https://dlmf.nist.gov/10.51 eq. 10.51.2 */
      class_call(bessel_at_x(pbs, z, index_l, &j1), pbs->error_message, pgb2->error_message);
      //class_call(bessel_j(pbs, pbs->l[index_l], z, &j1), pbs->error_message, pgb2->error_message);
      class_call(bessel_j(pbs, pbs->l[index_l]+1, z, &j2), pbs->error_message, pgb2->error_message);

      out = -j2 + (pbs->l[index_l]*j1)/z;

      *result = out;
    }

    return _SUCCESS_;
}

double bessel_at_x_first_deriv(struct galbispectra2 * pgb2,
                       struct bessels * pbs,
                       double z,
                       int index_l,
                       double * result){

    double out;
    double j1;
    double j2;

    if (z == 0.0) {
      *result = 0.0;
    }

    else{
     /*using recurrence relation https://dlmf.nist.gov/10.51 eq. 10.51.2 */
      //class_call(bessel_at_x(pbs, z, index_l, &j1), pbs->error_message, pgb2->error_message);
      class_call(bessel_j(pbs, pbs->l[index_l], z, &j1), pbs->error_message, pgb2->error_message);
      class_call(bessel_j(pbs, pbs->l[index_l]+1, z, &j2), pbs->error_message, pgb2->error_message);

      out = -j2 + (pbs->l[index_l]*j1)/z;

      *result = out;
    }

    return _SUCCESS_;
}

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
    double A,B,C,l;

    if (pbs->l[index_l] == 0){
      out = 0.0;
    }
    else if(pbs->l[index_l]==1){
      out = 0.0;
    }

    else if(z == 0.0){
      *result = 1.0;
      out = 0.0;
    }
    /*using recurrence relation */
    else{

      class_call(bessel_at_x(pbs, z, index_l, &j1), pbs->error_message, pgb2->error_message);

      class_call(bessel_at_x_first_deriv(pgb2,
                             pbs,
                             z,
                             index_l,
                             &j2),
                             pbs->error_message,
                             pgb2->error_message);

      A = -2.0*j2;
      B = -(z*z-pbs->l[index_l]*(pbs->l[index_l]+1.0))*j1;
      C = z*z;
      l = pbs->l[index_l];
      out = (-2*z*j2-(z*z-l*(l+1))*j1)/(z*z);


      *result = out;
    }

    return _SUCCESS_;
}

double bessel_at_x_second_deriv_old(struct galbispectra2 * pgb2,
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
    class_call(bessel_j(pbs, pbs->l[index_l], z, &j1), pbs->error_message, pgb2->error_message);
    class_call(bessel_j(pbs, pbs->l[index_l]+1, z, &j2), pbs->error_message, pgb2->error_message);

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
  double * pvecback1,
  double * pvecback2,
  double * result)
  {
  int index_md;
  double tmp = 0.;
  double lensing_result = 0.;


  double j1, j2;
  double f_evo1, f_evo2;
  double Pk;
  double j1l;
  double j1l_plus;
  double j2l;
  double j2l_plus;
  double type1, type2;

  /* Pair input type with corresponding source type. Each type is made up of, the transfer function of the perturbation, a prefactor, some k-factor,
    and a Bessel function (or a derivative of). For integrated terms inside pgb2->first_order_sources_integ, the spherical Bessel function
    may have already been integrated over before this stage, thus is will not appear explicitly here. */
  int last_index1 = 0;
  int last_index2 = 0;
  int index_tau_lens;
  double prefactor1, prefactor2;
  double lensing_result1, lensing_result2;
  double velocity1, velocity2;

  class_call(background_at_tau(pba,
                               pgb2->tau_sampling_cls[index_tau_first],
                               pba->long_info,
                               pba->inter_normal,
                               &last_index1,
                               pvecback1),
                               pba->error_message,
                               pgb2->error_message);

  class_call(background_at_tau(pba,
                              pgb2->tau_sampling_cls[index_tau_second],
                              pba->long_info,
                              pba->inter_normal,
                              &last_index2,
                              pvecback2),
                              pba->error_message,
                              pgb2->error_message);

  for(int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++){

    double x1 = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first]);
    double x2 = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second]);

    if(x1>pbs->x_max || x2 > pbs->x_max){
      printf("ALERT! x1= %g x2 = %g \n",x1,x2);continue;}

    /* Call primordial power spectrum (Fourier space)*/
    class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, pgb2->k_bessel[index_k_bessel], &Pk), ppm->error_message, pgb2->error_message);

    /* First Type: Density */

    if( index_type_first == pgb2->index_type_density){

      class_call(bessel_at_x(pbs, x1 , index_l, &j1), pbs->error_message, pgb2->error_message);

      type1 = pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau_first][index_k_bessel] * j1;

    }

    /* First Type: RSD */

    else if( index_type_first == pgb2->index_type_rsd ){

      class_call(bessel_at_x_second_deriv(pgb2, pbs, x1 , index_l, &j1), pbs->error_message, pgb2->error_message);

      prefactor1 = 1.0
                   /pvecback1[pba->index_bg_H]
                   /pvecback1[pba->index_bg_a];

      type1 = prefactor1
              *pgb2->first_order_sources[pgb2->index_source_theta][index_tau_first][index_k_bessel]
              *j1;
    }

    /* First Type: Doppler1 */

    else if( index_type_first == pgb2->index_type_d1 ){

      class_call(bessel_at_x_first_deriv(pgb2, pbs, x1 , index_l, &j1), pbs->error_message, pgb2->error_message);

      f_evo1 = 2.
               /pvecback1[pba->index_bg_H]
               /pvecback1[pba->index_bg_a]
               /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first])
               +pvecback1[pba->index_bg_H_prime]
               /pvecback1[pba->index_bg_H]
               /pvecback1[pba->index_bg_H]
               /pvecback1[pba->index_bg_a];

      prefactor1 = (1.
                    +pvecback1[pba->index_bg_H_prime]
                    /pvecback1[pba->index_bg_a]
                    /pvecback1[pba->index_bg_H]
                    /pvecback1[pba->index_bg_H]
                    +(2.-5.*ptr->s_bias)
                    /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first])
                    /pvecback1[pba->index_bg_a]
                    /pvecback1[pba->index_bg_H]
                    +5.*ptr->s_bias
                    -f_evo1);

      type1 = prefactor1
              *pgb2->first_order_sources[pgb2->index_source_theta][index_tau_first][index_k_bessel]
              *j1
              /pgb2->k_bessel[index_k_bessel];
    }

    /* First Type: Doppler2 */

    else if( index_type_first == pgb2->index_type_d2 ){

      class_call(bessel_at_x(pbs, x1 , index_l, &j1), pbs->error_message, pgb2->error_message);

      f_evo1 = 2.
               /pvecback1[pba->index_bg_H]
               /pvecback1[pba->index_bg_a]
               /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first])
               +pvecback1[pba->index_bg_H_prime]
               /pvecback1[pba->index_bg_H]
               /pvecback1[pba->index_bg_H]
               /pvecback1[pba->index_bg_a];
               //alert f_evo1 skipped
      //prefactor1 = -3.0*pvecback1[pba->index_bg_a]*pvecback1[pba->index_bg_H];

      type1 = (f_evo1-3.0)
              *pvecback1[pba->index_bg_a]
              *pvecback1[pba->index_bg_H]
              *pgb2->first_order_sources[pgb2->index_source_theta][index_tau_first][index_k_bessel]
              *j1
              /pgb2->k_bessel[index_k_bessel]
              /pgb2->k_bessel[index_k_bessel];
    }

    /* First Type: g1 (first of the GR terms) */

    else if( index_type_first == pgb2->index_type_g1){

      class_call(bessel_at_x(pbs, x1 , index_l, &j1), pbs->error_message, pgb2->error_message);


      prefactor1 = -2.0+5.0*ptr->s_bias;

      type1 = prefactor1
              *pgb2->first_order_sources[pgb2->index_source_phi][index_tau_first][index_k_bessel]
              *j1;
    }

    /* First Type: g2 */

    else if( index_type_first == pgb2->index_type_g2 ){

      class_call(bessel_at_x(pbs, x1 , index_l, &j1), pbs->error_message, pgb2->error_message);

      f_evo1 = 2.
               /pvecback1[pba->index_bg_H]
               /pvecback1[pba->index_bg_a]
               /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first])
               +pvecback1[pba->index_bg_H_prime]
               /pvecback1[pba->index_bg_H]
               /pvecback1[pba->index_bg_H]
               /pvecback1[pba->index_bg_a];

      prefactor1 = (2.0
                   +pvecback1[pba->index_bg_H_prime]
                   /pvecback1[pba->index_bg_H]
                   /pvecback1[pba->index_bg_H]
                   /pvecback1[pba->index_bg_a]
                   +(2.0-5.0*ptr->s_bias)
                   /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_first])
                   /pvecback1[pba->index_bg_H]
                   /pvecback1[pba->index_bg_a]
                   +5*ptr->s_bias
                   -f_evo1);

      type1 = prefactor1
              *pgb2->first_order_sources[pgb2->index_source_psi][index_tau_first][index_k_bessel]
              *j1;

    }

    /* First Type: g3 */
    else if( index_type_first == pgb2->index_type_g3 ){

      class_call(bessel_at_x(pbs, x1 , index_l, &j1), pbs->error_message, pgb2->error_message);


      prefactor1 = 1.0
                   /pvecback1[pba->index_bg_H]
                   /pvecback1[pba->index_bg_a];

      type1 = prefactor1
              *pgb2->first_order_sources[pgb2->index_source_phi_prime][index_tau_first][index_k_bessel]
              *j1;
    }

    /* First Type: g4 */

    else if( index_type_first == pgb2->index_type_g4 ){

      type2 = pgb2->first_order_sources_integ[pgb2->index_type_g4][index_l][index_tau_first][index_k_bessel];
    }

    /* First Type: g5 */

    else if( index_type_first == pgb2->index_type_g5 ){

      type2 = pgb2->first_order_sources_integ[pgb2->index_type_g5][index_l][index_tau_first][index_k_bessel];
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

    if( index_type_second == pgb2->index_type_density){

      class_call(bessel_at_x(pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);

      type2 = pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau_second][index_k_bessel] * j2;
      f_evo2 = 2.
              /pvecback2[pba->index_bg_H]
              /pvecback2[pba->index_bg_a]
              /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second])
              +pvecback2[pba->index_bg_H_prime]
              /pvecback2[pba->index_bg_H]
              /pvecback2[pba->index_bg_H]
              /pvecback2[pba->index_bg_a];

      /*type2 = j2
              *(ptr->s_bias*
              pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau_second][index_k_bessel]
              +(f_evo2-3.)
              *pvecback2[pba->index_bg_H]
              *pvecback2[pba->index_bg_a]
              *pgb2->first_order_sources[pgb2->index_source_v][index_tau_second][index_k_bessel]
              /pgb2->k_bessel[index_k_bessel]
              -(3.*pgb2->first_order_sources[pgb2->index_source_phi][index_tau_second][index_k_bessel]));*/


      /*type2 = j2
              *(pvecback2[pba->index_bg_Omega_m]
              +3.
              *pvecback2[pba->index_bg_H]
              *pvecback2[pba->index_bg_a]
              *pgb2->first_order_sources[pgb2->index_source_theta][index_tau_second][index_k_bessel]
              /pgb2->k_bessel[index_k_bessel]
              /pgb2->k_bessel[index_k_bessel]);*/

     // Using Poisson eq. A.8 in [1307.1459]
      /*type2 = j2
              *-2.
              *pgb2->k_bessel[index_k_bessel]
              *pgb2->k_bessel[index_k_bessel]
              *pgb2->first_order_sources[pgb2->index_source_phi][index_tau_second][index_k_bessel]
              /3.
              /pvecback2[pba->index_bg_H]
              /pvecback2[pba->index_bg_H]
              /pvecback2[pba->index_bg_a]
              /pvecback2[pba->index_bg_a]
              /pvecback2[pba->index_bg_Omega_m];*/
    }

    /* Second Type: RSD */

    else if( index_type_second == pgb2->index_type_rsd ){

      class_call(bessel_at_x_second_deriv(pgb2, pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);

      prefactor2 = 1.0
                   /pvecback2[pba->index_bg_H]
                   /pvecback2[pba->index_bg_a];

      type2 = prefactor2
              *pgb2->first_order_sources[pgb2->index_source_theta][index_tau_second][index_k_bessel]
              *j2;
    }

    /* Second Type: Doppler1 */

    else if( index_type_second == pgb2->index_type_d1 ){

      class_call(bessel_at_x_first_deriv(pgb2, pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);

      f_evo2 = 2.
              /pvecback2[pba->index_bg_H]
              /pvecback2[pba->index_bg_a]
              /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second])
              +pvecback2[pba->index_bg_H_prime]
              /pvecback2[pba->index_bg_H]
              /pvecback2[pba->index_bg_H]
              /pvecback2[pba->index_bg_a];

      prefactor2 = (1.
                    +pvecback2[pba->index_bg_H_prime]
                    /pvecback2[pba->index_bg_a]
                    /pvecback2[pba->index_bg_H]
                    /pvecback2[pba->index_bg_H]
                    +(2.-5.*ptr->s_bias)
                    /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second])
                    /pvecback2[pba->index_bg_a]
                    /pvecback2[pba->index_bg_H]
                    +5.*ptr->s_bias
                    -f_evo2);

      type2 = prefactor2
              *pgb2->first_order_sources[pgb2->index_source_theta][index_tau_second][index_k_bessel]
              *j2
              /pgb2->k_bessel[index_k_bessel];
    }

    /* Second Type: Doppler2 */

    else if( index_type_second == pgb2->index_type_d2 ){

      class_call(bessel_at_x(pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);

      f_evo2 = 2.
               /pvecback2[pba->index_bg_H]
               /pvecback2[pba->index_bg_a]
               /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second])
               +pvecback2[pba->index_bg_H_prime]
               /pvecback2[pba->index_bg_H]
               /pvecback2[pba->index_bg_H]
               /pvecback2[pba->index_bg_a];
               //alert f_evo skipped
      //prefactor2 = -3.0*pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];



      type2 = (f_evo2-3.0)
              *pvecback2[pba->index_bg_a]
              *pvecback2[pba->index_bg_H]
              *pgb2->first_order_sources[pgb2->index_source_theta][index_tau_second][index_k_bessel]
              *j2
              /pgb2->k_bessel[index_k_bessel]
              /pgb2->k_bessel[index_k_bessel];
    }

    /* Second Type: g1 (first of the GR terms) */

    else if( index_type_second == pgb2->index_type_g1 ){

      class_call(bessel_at_x(pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);


      prefactor2 = -2.0+5.0*ptr->s_bias;

      type2 = prefactor2
              *pgb2->first_order_sources[pgb2->index_source_phi][index_tau_second][index_k_bessel]
              *j2;
    }

    /* Second Type: g2 */

    else if( index_type_second == pgb2->index_type_g2 ){

      class_call(bessel_at_x(pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);

      f_evo2 = 2.
               /pvecback2[pba->index_bg_H]
               /pvecback2[pba->index_bg_a]
               /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second])
               +pvecback2[pba->index_bg_H_prime]
               /pvecback2[pba->index_bg_H]
               /pvecback2[pba->index_bg_H]
               /pvecback2[pba->index_bg_a];

      prefactor2 = (2.0
                   +pvecback2[pba->index_bg_H_prime]
                   /pvecback2[pba->index_bg_H]
                   /pvecback2[pba->index_bg_H]
                   /pvecback2[pba->index_bg_a]
                   +(2.0-5.0*ptr->s_bias)
                   /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau_second])
                   /pvecback2[pba->index_bg_H]
                   /pvecback2[pba->index_bg_a]
                   +5*ptr->s_bias
                   -f_evo2);

      type2 = prefactor2
              *pgb2->first_order_sources[pgb2->index_source_psi][index_tau_second][index_k_bessel]
              *j2;

    }

    /* Second Type: g3 */
    else if( index_type_second == pgb2->index_type_g3 ){

      class_call(bessel_at_x(pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);


      prefactor2 = 1.0
                   /pvecback2[pba->index_bg_H]
                   /pvecback2[pba->index_bg_a];

      type2 = prefactor2
              *pgb2->first_order_sources[pgb2->index_source_phi_prime][index_tau_second][index_k_bessel]
              *j2;
    }

    /* Second Type: g4 */

    else if( index_type_second == pgb2->index_type_g4 ){

      type2 = pgb2->first_order_sources_integ[pgb2->index_type_g4][index_l][index_tau_second][index_k_bessel];
    }

    /* Second Type: g5 */

    else if( index_type_second == pgb2->index_type_g5 ){

      type2 = pgb2->first_order_sources_integ[pgb2->index_type_g5][index_l][index_tau_second][index_k_bessel];
    }



    /* Second Type: Lensing convergence */

    else if( index_type_second == pgb2->index_type_lens ){

      prefactor2 = -ptr->l[index_l] * (ptr->l[index_l] + 1.);

      type2 = prefactor2 * pgb2->first_order_sources_integ[pgb2->index_type_lens][index_l][index_tau_second][index_k_bessel];
    }

    else {
      type2 = 0.;
    }


    // DEBUG
    //printf("index_l x index_tau_first x index_tau_second x index_k = %dx%dx%dx%d\n", index_l, index_tau_first, index_tau_second, index_k_bessel);
    //printf("j1 = %g\n", j1);
    //printf("j2 = %g\n", j2);
    //printf("type1 = %g\n", type1);
    //printf("pgb2->first_order_sources_integ[pgb2->index_source_theta][%d][%d][%d] = %g\n", index_l, index_tau_first, index_k_bessel, pgb2->first_order_sources_integ[pgb2->index_source_theta][index_l][index_tau_first][index_k_bessel]);
    //printf("type2 = %g\n", type2);
    //printf("pgb2->first_order_sources_integ[pgb2->index_source_theta][%d][%d][%d] = %g\n", index_l, index_tau_second, index_k_bessel, pgb2->first_order_sources_integ[pgb2->index_source_theta][index_l][index_tau_second][index_k_bessel]);
    //printf(" \n");



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


int index_of_tau_sampling_cls(double tau,
                 int * index_tau,
                 struct galbispectra2 * pgb2){

                 double output = (tau - pgb2->tau_sampling_cls[0])/(pgb2->tau_sampling_cls[pgb2->tau_size_cls-1] - pgb2->tau_sampling_cls[0]) * (pgb2->tau_size_cls-1);
                 //printf("output = %g\n",output );
                 * index_tau = (int) ceil(output);

                 if (* index_tau < 1) {
                   * index_tau = 1;
                 }

                 else if(* index_tau > pgb2->tau_size_cls -1){
                   * index_tau = pgb2->tau_size_cls -1;

                 }
                 //printf("index_given = %d\n",*index_tau );
               }
/* index_of_tau_sampling_cls() - given some conformal time tau1, this function will output the (closest)
 corresponding index in the pgb2->tau_sampling_selection[bin][index_tau1] array.

* @param tau1                     Input: conformal time
* @param index_tau1               Output: index in the array tau_sampling_cls
* @param galbispectra2            Input: galbispectra2 structure
*/
int index_of_tau_sampling_selection(double tau,
                int bin,
                int * index_tau,
                struct galbispectra2 * pgb2){

                double output = (tau - pgb2->tau_sampling_selection[bin][0])/(pgb2->tau_sampling_selection[bin][pgb2->tau_size_selection-1] - pgb2->tau_sampling_selection[bin][0]) * (pgb2->tau_size_selection-1);
                //Warning changed to floor from ceil
                * index_tau = (int) ceil (output);

                if (* index_tau < 0) {
                  * index_tau = 0;
                }

                else if(* index_tau > pgb2->tau_size_selection -1){
                  * index_tau = pgb2->tau_size_selection -1;

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
int index_of_alpha(double alpha,
                 int * index_alpha,
                 struct galbispectra2 * pgb2){

                 double output = (alpha - pgb2->alpha[0])/(pgb2->alpha[pgb2->alpha_size-1] - pgb2->alpha[0]) * (pgb2->alpha_size-1);
                 //printf("output = %g\n",output );
                 * index_alpha = (int) ceil(output);

                 if (* index_alpha < 1) {
                   * index_alpha = 1;
                 }

                 else if(* index_alpha > pgb2->alpha_size -1){
                   * index_alpha = pgb2->alpha_size -1;

                 }
               }

 int index_of_r(double r,
                int index_alpha,
                int * index_r,
                struct galbispectra2 * pgb2){

                  double output = (r - pgb2->r2[index_alpha][0])/(pgb2->r2[index_alpha][pgb2->r_size-1] - pgb2->r2[index_alpha][0]) * (pgb2->r_size-1);
                  //printf("output = %g\n",output );
                  * index_r = (int) ceil(output);

                  if (* index_r < 1) {
                    * index_r = 1;
                  }

                  else if(* index_r > pgb2->r_size -1){
                    * index_r = pgb2->r_size -1;

                  }
                }

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


int compute_simps_weights_selection(int bin,
                double * weights,
                struct galbispectra2 * pgb2){

                int size = pgb2->tau_size_selection;
                double delta_x;


                delta_x = (pgb2->tau_sampling_selection[bin][size-1]-pgb2->tau_sampling_selection[bin][0])/(size-1);



                /* End weights */
                weights[0] = delta_x/3.0;
                weights[size-1] = delta_x/3.0;
                /* Even sums */
                for (int index = 1; index < (size-3)/2+1; index++) {
                  weights[2*index] = 2.0*delta_x/3.0;
                  //printf("even inside function weights[%d]=%g\n",index, weights[index]);
                }
                /* Odd sums */
                for (int index = 1; index < (size-3)/2+2; index++) {
                  weights[2*index-1] = 4.0*delta_x/3.0;
                  //printf("odd inside function weights[%d]=%g\n",index, weights[index]);
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
  //int index_tau_first;
  //int index_tau_second;
  int index_type;
  int index_md;
  int bin_first;
  int bin_second;
  double * pvecback;
  double * w_trapz_tau;
  double tau0;
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
  /*the bessel grid is boosted by an integer factor to account for rapid oscillations.*/
  int bessel_boost = 50;
  tau0 = pba->conformal_age;
  /* SIZES */


  /* Please ensure the pgb2->tau_size_selection grid is of a higher resolution than pgb2->r_size * pgb2->alpha_size = 13*/
  pgb2->tau_size_cls = 750; //prev on 750
  pgb2->tau_size_selection = 601;
  pgb2->r_size = 150; //previously on 75

  /* alpha_size must be an odd positive-integer in order to correctly fill the pgb2->r and pgb2->alpha grids using symmetries */
  pgb2->alpha_size = 65; //previously on 13

  if (pgb2->alpha_size % 2 == 0) {
    printf("ERROR! alpha_size must be an ODD positive-integer\n" );
    exit(0);
  }

  /* We wish to double the number of slots such that points in the two grids tau_sampling_cls and tau_sampling_bessel align. */
  pgb2->tau_size_bessel = bessel_boost * (pgb2->tau_size_cls-1) + 1;
  pgb2->k_size_bessel = 25000;  /*tau_size_bessel*/   /*k_size * boost*/
  printf("Starting galaxy bispectra module...\n");
  if (ppt->selection==gaussian) {
    printf("We have a Gaussian Window Function on our hands.\n");
  }
  if (ppt->selection==dirac) {
    printf("We have a Dirac Window Function on our hands.\n");
  }

  double * result;

  clock_t begin = clock();

/* here, do your time-consuming job */




/************************************************************
 ********* Second-order matter kernel computation *********
************************************************************/
/* Uncomment the following lines to print the second-order CDM kernel */
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


/* Time-Sampling */
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
              pgb2->tau_size_cls * sizeof(double),
              pgb2->error_message);

  class_alloc(pgb2->tau_sampling_bessel,
              pgb2->tau_size_bessel * sizeof(double),
              pgb2->error_message);

  class_alloc(pgb2->r,
              pgb2->r_size * sizeof(double),
              pgb2->error_message);

  class_alloc(pgb2->alpha,
              pgb2->alpha_size * sizeof(double),
              pgb2->error_message);

  class_alloc(pgb2->r2,
              pgb2->alpha_size * sizeof(double*),
              pgb2->error_message);

  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {

    class_alloc(pgb2->r2[index_alpha],
              pgb2->r_size * sizeof(double),
              pgb2->error_message);
  }



  double overall_tau_min = 160000.;
  double overall_tau_max = -1.;


  /* Define new tau_sampling_selection. This sampling is bin-dependent */
  double selection_mean;
  for (bin = 0; bin < ppt->selection_num; bin++) {
    //finer sampling of bins
    double tau_min;
    double tau_max;
    double z_max, z_min;
     //printf("ppt->selection_mean[%d]=%g, width = %g\n",bin,ppt->selection_mean[bin], ppt->selection_width[bin] );
    z_max =  ppt->selection_mean[bin]+ 5. * ppt->selection_width[bin];
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
                            &tau_min),
                            ppt->error_message,
                            ppt->error_message);

    class_call(background_tau_of_z(
                            pba,
                            z_min,
                            &tau_max),
                            ppt->error_message,
                            pgb2->error_message);

    class_call(background_tau_of_z(
                            pba,
                            ppt->selection_mean[bin],
                            &selection_mean),
                            ppt->error_message,
                            pgb2->error_message);

    if (tau_max > ppt->tau_sampling_quadsources[ppt->tau_size_quadsources-1]) {
      tau_max = ppt->tau_sampling_quadsources[ppt->tau_size_quadsources-1];
    }


    if (tau_min < 400.) {
      printf("Rubbish, Window function extends to before recombination\n");
      tau_min = 400.;
    }

    overall_tau_min = MIN(tau_min,overall_tau_min);
    overall_tau_max = MAX(tau_max,overall_tau_max);
    // ALERT SHOULD THE TIME GRIDS RUN TO conformal_age or tau_max
    for (index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
      pgb2->tau_sampling_selection[bin][index_tau] = tau_min + index_tau*(tau_max-tau_min)/(pgb2->tau_size_selection-1);
      //printf("tau[%d][%d] = %g\n", bin, index_tau, pgb2->tau_sampling_selection[bin][index_tau]);
    }
  }


  // ALERT HARDCODED!!
  for (index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {
    //pgb2->tau_sampling_cls[index_tau] = overall_tau_min + index_tau*(pba->conformal_age-overall_tau_min)/(pgb2->tau_size_cls-1);
    pgb2->tau_sampling_cls[index_tau] = overall_tau_min + index_tau*(overall_tau_max-overall_tau_min)/(pgb2->tau_size_cls-1);

  }

  /* Find the ranges of alpha to span such that we sample within tau_sampling_cls[0] to tau_sampling_cls[pgb2->tau_size_selection-1]*/
  double tau_max_cls = pgb2->tau_sampling_cls[pgb2->tau_size_cls-1];
  double tau_min_cls = pgb2->tau_sampling_cls[0];
  double r_max;
  double r_min = sqrt((1.-(11750./pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))*(1.-(11750./pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))+(1.-(11750./pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))*(1.-(11750./pgb2->tau_sampling_cls[pgb2->tau_size_cls-1])));
  double r_min2 = 0.2411267469;

  printf("r_min2 = %g\n", r_min2);





  /* Under new parameterisation tau1 = tau0(1-rsin(alpha)). tau2 = tau0(1-rcos(alpha)). Max r value can be found by
  rearranging tau_min = tau1(r_max,pi/2)=tau2(r_max,0). This is assuming alpha runs over 0 to pi/2. */
  // NOTE this has been edited to be hardcoded



  printf("r spans (%g,%g) with %d points\n",r_min, r_max, pgb2->r_size);


  printf("tau_sampling_cls has %d points that span (%g, %g)\n", pgb2->tau_size_cls, pgb2->tau_sampling_cls[0], pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]);
  printf("conformal_age = %g\n",pba->conformal_age);

  // Time-grid values for tau1 = tau0(1-r*sin(alpha)), tau2 = tau0(1-r*cos(alpha)). The motivation for this parameterisation is to
  // effectively sample highly oscillatory regions of upcoming integrands. We sample one dimension with fine precision and the other
  // with coarse precision. This way when we integrate over alpha and r we have fewer iterations since alpha_size * r_size < tau_size_selection^2.



  // NOTE: THIS HAS BEEN CHANGED
  double alpha_min = 0.6071841396;
  //double alpha_min = 0.0;

  double alpha_max = _PI_/2.0-alpha_min;
  //double alpha_max = _PI_/2.0;

  int middle_index_alpha = (pgb2->alpha_size-1)/2;
  printf("middle_index_alpha = %d\n", middle_index_alpha );

  printf("alpha_min =%g,  %g\n", alpha_min, (tau0-tau_max_cls)/(r_min2*tau_max_cls));
  printf("alpha_max =%g\n", alpha_max );



  /*for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    // Linear spacing
    pgb2->alpha[index_alpha] =  alpha_min+index_alpha*(alpha_max-alpha_min)/(pgb2->alpha_size-1);


  }*/

  double base = 1.3;
  double alpha_middle = alpha_min + (alpha_max-alpha_min)/2.;

  for (int index_alpha = middle_index_alpha-1; index_alpha > -1; index_alpha--) {

    pgb2->alpha[index_alpha] =  alpha_middle-pow(base,middle_index_alpha-index_alpha)*(alpha_middle-alpha_min)/(pow(base,middle_index_alpha));


  }
  pgb2->alpha[middle_index_alpha] = alpha_middle;
  // log spacing in z1=z2=1 case


  for (int index_alpha = middle_index_alpha+1; index_alpha < pgb2->alpha_size; index_alpha++) {
    //printf("index_alpha = %d\n", index_alpha );
    pgb2->alpha[index_alpha] =  alpha_middle+pow(base,index_alpha)*(alpha_max-alpha_middle)/(pow(base,pgb2->alpha_size-1));
    //printf("upper alpha[%d] = %g\n", index_alpha, pgb2->alpha[index_alpha]);

  }
  double alpha_sum = 0.;
  printf("#cumulative sum     alpha\n");
  //printf("#base = %g, alpha_size =%d\n", base, pgb2->alpha_size);
  /*for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    alpha_sum += pgb2->alpha[index_alpha];
    printf("%g      %g\n", pgb2->alpha[index_alpha], alpha_sum);
  }
  exit(0);
  */










  printf("tau_max_cls = %g\n", pgb2->tau_sampling_cls[pgb2->tau_size_selection-1] );



  for (int index_tau_bessel = 0; index_tau_bessel < pgb2->tau_size_bessel; index_tau_bessel++) {
    pgb2->tau_sampling_bessel[index_tau_bessel] = overall_tau_min + index_tau_bessel*(pba->conformal_age-overall_tau_min)/(pgb2->tau_size_bessel-1);
  }

  /* New Bessel k-sampling to capture features of the Bessel oscillations */

  class_alloc(pgb2->k_bessel,
              pgb2->k_size_bessel * sizeof(double),
              pgb2->error_message);
  printf("Allocated k_bessel array\n");
  // NOTE: this part may need to be more general (i.e. not fixed to scalars)

  double k_min = ppt->k[ppt->index_md_scalars][0];

  double k_max =  ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1];

  for (int i = 0; i < pgb2->k_size_bessel; i++) {

    pgb2->k_bessel[i] = k_min + i*(k_max-k_min)/pgb2->k_size_bessel;
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






  /* Allocate and fill array for the trapezoidal weights for chi integration in the lensing term w_trapz_lens[index_tau] */
  class_alloc(pgb2->w_trapz_lens,
              pgb2->tau_size_bessel * sizeof(double*),
              ppt->error_message);

  /*class_call(array_trapezoidal_weights(pgb2->r,
                                       pgb2->r_size,
                                       w_trapz_r,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);*/
  double * w_trapz_alpha;

  class_alloc(w_trapz_alpha,
              pgb2->alpha_size * sizeof(double),
              ppt->error_message);


  class_call(array_trapezoidal_weights(pgb2->alpha,
                                       pgb2->alpha_size,
                                       w_trapz_alpha,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);

  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    printf("%g      %g\n",pgb2->alpha[index_alpha], w_trapz_alpha);
  }


  class_call(array_trapezoidal_weights(pgb2->tau_sampling_bessel,
                                       pgb2->tau_size_bessel,
                                       pgb2->w_trapz_lens,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);


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
    }

    class_call(array_trapezoidal_mweights(tau0_minus_tau[bin]/*pgb2->tau_sampling_selection[bin]*/,
                                          pgb2->tau_size_selection,
                                          w_trapz[bin],
                                          pgb2->error_message),
                                          ppt2->error_message,
                                          ppt2->error_message);
  }


  //free(tau0_minus_tau);
  /* Computing the window-function */
  /* Declaration of temporary pointer */
  double ** selection;

  /* Allocation of first dimension selection[bin] */
  class_alloc(selection,
              ppt->selection_num * sizeof(double*),
              ppt->error_message);

  printf("selection_num = %d\n",ppt->selection_num);
  class_alloc(pvecback,
              pba->bg_size * sizeof(double),
              ptr->error_message);
  double ** selection_test;

  class_alloc(selection_test,
              ppt->selection_num * sizeof(double*),
              pgb2->error_message);

  /* Allocation of second dimension selection[bin][index_tau] */
  for(int bin = 0; bin < ppt->selection_num; bin++){
    printf("bin = %d\n", bin );
    class_alloc(selection[bin],
                pgb2->tau_size_selection * sizeof(double),
                ppt->error_message);

    class_alloc(selection_test[bin],
                pgb2->tau_size_selection * sizeof(double),
                ppt->error_message);

    /* transfer_selection_compute writes in to selection[bin] */
    class_call(transfer_selection_compute(ppr,
                                          pba,
                                          ppt,
                                          ptr,
                                          selection[bin],
                                          tau0_minus_tau[bin]/*pgb2->tau_sampling_selection[bin]*/,
                                          w_trapz[bin],
                                          pgb2->tau_size_selection,
                                          pvecback,
                                          tau0,
                                          bin),
               pgb2->error_message,
               pgb2->error_message);

  }

  /* Create an array to store the tau_selection indices that contain 66% of the weight of the window function between their
    corresponding time values*/

  int middle_tau_selection = (pgb2->tau_size_selection+1)/2-1;
  double selection_mean_tau;
  int selection_mean_tau_index;

  /* tau_polar_min and tau_polar_max are the maxiumum and minimum conformal time values spanned by any pair of bins using
  the polar parameterisation tau1 = tau0(1-r*sin(alpha)), tau1 = tau0(1-r*sin(alpha)). This parameterisation is used to
  integrate over time to calculate angular power spectra. */

  double ** tau_polar_min;
  double ** tau_polar_max;




  class_alloc(tau_polar_min,
            ppt->selection_num * sizeof(double*),
            pgb2->error_message);

  class_alloc(tau_polar_max,
            ppt->selection_num * sizeof(double*),
            pgb2->error_message);

  for (bin = 0; bin < ppt->selection_num; bin++) {
    class_alloc(tau_polar_min[bin],
              ppt->selection_num * sizeof(double),
              pgb2->error_message);

    class_alloc(tau_polar_max,
              ppt->selection_num * sizeof(double),
              pgb2->error_message);
  }

  /* The bin cut-off corresponds to the smaller and larger index in the pgb2->tau_sampling_selection[bin] grid, which incapsulates
  a chosen proportion of the integrated weight of the Gaussian window function. */
  /*int ** bin_cut_off;

  class_alloc(bin_cut_off,
            ppt->selection_num * sizeof(int*),
            pgb2->error_message);

  for (bin = 0; bin < ppt->selection_num; bin++) {

    class_call(background_tau_of_z(
                            pba,
                            ppt->selection_mean[bin],
                            &selection_mean_tau),
                            ppt->error_message,
                            pgb2->error_message);

    printf("selection_mean_tau = %g\n",selection_mean_tau );
    index_of_tau_sampling_selection(selection_mean_tau,
                    bin,
                    &selection_mean_tau_index,
                    pgb2);
    printf("selection_mean_tau_index = %d\n",selection_mean_tau_index );
    class_alloc(bin_cut_off,
              2 * sizeof(int),    // Upper and lower index
              pgb2->error_message);
    for (int i = 0; i < 2; i++) {
      double left_sum = 0.0;
      double right_sum = 0.0;

      for (int index_tau_left = selection_mean_tau_index; index_tau_left > -1; index_tau_left--) {
        //printf("i = %d\n", index_tau_left );
        left_sum += selection[bin][index_tau_left]*w_trapz[bin][index_tau_left];
        //printf("left_sum = %g\n", left_sum );
        if (left_sum > 0.475){
          bin_cut_off[bin][0] = index_tau_left;
          //printf("*index_tau_left %d*\n", index_tau_left );
          break;
        }
      }

      for (int index_tau_right = selection_mean_tau_index; index_tau_right < pgb2->tau_size_selection; index_tau_right++) {
        //printf("i = %d\n", index_tau_right );
        right_sum += selection[bin][index_tau_right]*w_trapz[bin][index_tau_right];
        //printf("right_sum = %g\n", right_sum );
        if (right_sum > 0.475){
          bin_cut_off[bin][1] = index_tau_right;
          printf("index_tau_right = %d\n", index_tau_right );
          printf("bin_cut_off[%d]][%d] = %d\n", bin, 1, bin_cut_off[bin][1]);
          break;
        }
      }
      //printf("bin_cut_off left: %d, bin_cut_off right: %d\n", bin_cut_off[bin][0], bin_cut_off[bin][1] );
    }
  }*/





  /* Absolute maximum value r can have is when alpha =pi/4 and tau is set to tau_min. If r1 or r2 are above this value, this is because of a sing-
  ularity in the denominator */
  double tau_one;
  double tau_two;
  double tau_one_cls;
  double tau_two_cls;
  int index_tau_one_cls;
  int index_tau_two_cls;
  double  r_abs_max;

  /* define the minimum distance between tau0-tau0 and the tau1_upper_cut_off-tau2_upper_cut_off, for each tau1-tau2 pair;
  that means permuting over bins. */
  double ** r_min3;

  class_alloc(r_min3,
              ppt->selection_num * sizeof(double*),
              pgb2->error_message);

  for (int i = 0; i < ppt->selection_num; i++) {
    class_alloc(r_min3[i],
                ppt->selection_num * sizeof(double),
                pgb2->error_message);
  }

  double **** r_bins;
  double *** alpha2;

  class_alloc(r_bins,
              ppt->selection_num * sizeof(double***),
              pgb2->error_message);

  class_alloc(alpha2,
              ppt->selection_num * sizeof(double**),
              pgb2->error_message);

  for (int i = 0; i < ppt->selection_num; i++) {
    class_alloc(r_bins[i],
                ppt->selection_num * sizeof(double**),
                pgb2->error_message);
    class_alloc(alpha2[i],
                ppt->selection_num * sizeof(double*),
                pgb2->error_message);
    for (int j = 0; j < ppt->selection_num; j++) {
      class_alloc(r_bins[i][j],
                  pgb2->alpha_size * sizeof(double*),
                  pgb2->error_message);
      class_alloc(alpha2[i][j],
                  pgb2->alpha_size * sizeof(double),
                  pgb2->error_message);
      for (int k = 0; k < pgb2->alpha_size; k++) {
        class_alloc(r_bins[i][j][k],
                    pgb2->r_size * sizeof(double),
                    pgb2->error_message);
      }
    }
  }

  double ** alpha_of_r_min;
  double ** r_abs_max3;

  class_alloc(alpha_of_r_min,
              ppt->selection_num * sizeof(double*),
              pgb2->error_message);

  class_alloc(r_abs_max3,
              ppt->selection_num * sizeof(double*),
              pgb2->error_message);

  for (int i = 0; i < ppt->selection_num; i++) {
    class_alloc(alpha_of_r_min[i],
                ppt->selection_num * sizeof(double),
                pgb2->error_message);

    class_alloc(r_abs_max3[i],
                ppt->selection_num * sizeof(double),
                pgb2->error_message);
  }

  int * bin_cut_off_upper;
  int * bin_cut_off_lower;

  class_alloc(bin_cut_off_upper,
            ppt->selection_num * sizeof(int*),
            pgb2->error_message);

  class_alloc(bin_cut_off_lower,
            ppt->selection_num * sizeof(int*),
            pgb2->error_message);

  int upper_cutoff;
  int lower_cutoff;


  /*for (int i = 0; i < ppt->selection_num; i++) {
    class_alloc(bin_cut_off[i],
              2 * sizeof(int),    // Upper and lower index
              pgb2->error_message);
  }*/


  double left_sum = 0.0;
  double right_sum = 0.0;

  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {

    class_call(background_tau_of_z(
                            pba,
                            ppt->selection_mean[bin1],
                            &selection_mean_tau),
                            ppt->error_message,
                            pgb2->error_message);

    index_of_tau_sampling_selection(selection_mean_tau,
                    bin1,
                    &selection_mean_tau_index,
                    pgb2);

    //printf("selection_mean_tau_index = %d\n",selection_mean_tau_index );
    for (int index_tau_left = selection_mean_tau_index; index_tau_left > -1; index_tau_left--) {
      //printf("i = %d\n", index_tau_left );
      left_sum += selection[bin1][index_tau_left]*w_trapz[bin1][index_tau_left];
      //printf("left_sum = %g\n", left_sum );
      if (left_sum > 0.475){
        bin_cut_off_lower[bin1] = index_tau_left;
        lower_cutoff = index_tau_left;
        //printf("*index_tau_left %d*\n", index_tau_left );
        break;
      }
    }
    for (int index_tau_right = selection_mean_tau_index; index_tau_right < pgb2->tau_size_selection; index_tau_right++) {

      right_sum += selection[bin1][index_tau_right]*w_trapz[bin1][index_tau_right];
      //printf("right_sum = %g\n", right_sum );
      if (right_sum > 0.475){
        bin_cut_off_upper[bin1] = index_tau_right;
        upper_cutoff = index_tau_right;

        break;
      }
    }
  }




  /* We are going to loop over the window-bins and create our polar grid. There will be a unique grid between any two pairs of window bins,
    this is to ensure that the time-space is adequately sampled within the window-regions. We first find the full-range that each window-regions function spans */
  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    for (int bin2 = 0; bin2 < bin1+1; bin2++) {
      printf("1. bin_cut_off[%d]][%d] = %d\n", bin1, bin_cut_off_upper[bin1]);
      printf("2. bin_cut_off[%d]][%d] = %d\n", bin2, bin_cut_off_upper[bin2]);
      double tau1_max = pgb2->tau_sampling_selection[bin1][pgb2->tau_size_selection-1];
      double tau2_max = pgb2->tau_sampling_selection[bin2][pgb2->tau_size_selection-1];
      double tau1_min = pgb2->tau_sampling_selection[bin1][0];
      double tau2_min = pgb2->tau_sampling_selection[bin2][0];

      printf("3. bin_cut_off[%d] = %d\n", bin1, bin_cut_off_upper[bin1]);
      double tau1_cutoff_upper = pgb2->tau_sampling_selection[bin1][bin_cut_off_upper[bin1]];
      double tau1_cutoff_lower = pgb2->tau_sampling_selection[bin1][bin_cut_off_lower[bin1]];
      double tau2_cutoff_upper = pgb2->tau_sampling_selection[bin2][bin_cut_off_upper[bin2]];
      double tau2_cutoff_lower = pgb2->tau_sampling_selection[bin2][bin_cut_off_lower[bin2]];


      double alpha_of_r_min = atan((tau0-tau1_cutoff_upper)/(tau0-tau2_cutoff_upper));

      double alpha_of_r_max = atan((tau0-tau1_cutoff_lower)/(tau0-tau2_cutoff_lower));

      double alpha_min = atan((tau0-tau1_cutoff_upper)/(tau0-tau2_cutoff_lower));

      double alpha_max = atan((tau0-tau1_cutoff_lower)/(tau0-tau2_cutoff_upper));



      if(isnan(alpha_of_r_min)) {
        alpha_of_r_min = _PI_/2.;
      }
      if(isnan(alpha_of_r_max)) {
        alpha_of_r_max = _PI_/2.;
      }
      if(isnan(alpha_max)) {
        alpha_max = _PI_/2.;
      }
      if(isnan(alpha_min)) {
        alpha_min = _PI_/2.;
      }
      /*for (int index_alpha = 0; index_alpha < (pgb2->alpha_size+1)/2; index_alpha++) {
        r_max_per_alpha[index_alpha] = (tau0-pgb2->tau_sampling_cls[0])/(tau0*cos(pgb2->alpha[index_alpha]));
      }

      r_max_per_alpha[middle_alpha_index] = r_abs_max;

      for (int index_alpha = middle_alpha_index+1; index_alpha < pgb2->alpha_size; index_alpha++) {
      for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
        alpha2[bin1][bin2][index_alpha] =  alpha_min+index_alpha*(alpha_max-alpha_min)/(pgb2->alpha_size-1);
      }*/

      double r_min = (1.-tau1_cutoff_upper/tau0)/sin(alpha_of_r_min);
      double r_abs_max = (1.-tau1_cutoff_lower/tau0)/sin(alpha_of_r_max);
      double tau_check = tau0*(1.-r_min*cos(alpha_of_r_min));
      printf("tau2_check = %g should be %g\n", tau_check, tau2_cutoff_upper);

      printf("bin1 = %d, bin2 = %d, bin_cut_off = %d\n", bin1, bin2, bin_cut_off_upper[bin1] );
      //alpha_of_r_min[bin1][bin2] = atan((tau0-pgb2->tau_sampling_selection[bin1][bin_cut_off_upper[bin1]])/(tau0-pgb2->tau_sampling_selection[bin2][bin_cut_off_upper[bin2]]));

      //r_abs_max3[bin1][bin2] = (1.-tau_min_cls/tau0)/sin(_PI_/4);
      //r_min3[bin1][bin2] = (1.-pgb2->tau_sampling_selection[bin1][bin_cut_off_upper[bin1]]/tau0)/sin(alpha_of_r_min/*[bin1][bin2]*/);
      printf("r_min3[%d][%d] = %g\n", bin1, bin2, r_min3[bin1][bin2]);

      /*if (tau1_max > tau2_max) {
        tau_polar_max[bin1][bin2] = tau1_max;
        tau_polar_max[bin2][bin1] = tau1_max;
      }

      else{
        tau_polar_max[bin1][bin2] = tau2_max;
        tau_polar_max[bin2][bin1] = tau2_max;
      }

      if (tau1_min < tau2_min) {
        tau_polar_min[bin1][bin2] = tau1_min;
        tau_polar_min[bin2][bin1] = tau1_min;
      }

      else{
        tau_polar_min[bin1][bin2] = tau2_min;
        tau_polar_min[bin2][bin1] = tau2_min;
      }*/
      /* There is a rectangular (square) box that is cut from time-space for cross(auto)-bin correlations using the tau_cutoff
      quantities. This is the region we wish to integrate which contains the weight of the window function. We want to distribute
      the alpha "sun-rays" such that there is a proportionate number of them on the longer/shorter side of the rectangle */

      // NOTE: alpha2 needs to be defined. Must incorporate alpha_of_r_max!
      double ratio =(tau1_cutoff_upper-tau1_cutoff_lower)/(tau2_cutoff_upper-tau2_cutoff_lower);
      if (ratio > 0.0){

        int horizontal_size = ceil((pgb2->alpha_size-1)*(ratio/(ratio+1.0)));
        for (int index_alpha = 0; index_alpha < horizontal_size; index_alpha++) {
          double increment_h = (alpha_of_r_max-alpha_min)/(horizontal_size-1);
          double r_max_per_alpha = (1.0-1.0/tau2_cutoff_lower)/sin(alpha2[bin1][bin2][index_alpha]);
          alpha2[bin1][bin2][index_alpha] = alpha_min + index_alpha*((alpha_of_r_max-increment_h)-alpha_min)/(horizontal_size-1);

          for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            r_bins[bin1][bin2][index_alpha][index_r] = r_min+index_r*(r_max_per_alpha-r_min)/(pgb2->r_size-1);
          }
        }
        /* We want a sun-ray to go directly to (tau1,tau2)=(tau1_cutoff_lower,tau2_cutoff_lower) (along the r_max line for this pair of
      bins. Note that pgb2->alpha_size = horizontal+1+vertical. The 1 referes to this line. r_abs_max refers to the largest r value along this
      given sunray, alpha_of_r_max is the corresponding angle (alpha).*/
        alpha2[bin1][bin2][horizontal_size-1] = alpha_of_r_max;
        for (int index_r = 0; index_r < pgb2->r_size; index_r++) {

          r_bins[bin1][bin2][horizontal_size-1][index_r] = r_min+index_r*(r_abs_max-r_min)/(pgb2->r_size-1);
        }
        int vertical_size = pgb2->alpha_size-1-horizontal_size;
        for (int index_alpha = horizontal_size+1; index_alpha < pgb2->alpha_size; index_alpha++) {
          double increment_v = (alpha_max-alpha_min)/(vertical_size-1);
          alpha2[bin1][bin2][index_alpha] = alpha_of_r_max + increment_v+index_alpha*(alpha_max-(alpha_of_r_max+increment_v))/(vertical_size-1);
          double r_max_per_alpha = (1.0-tau0/tau2_cutoff_lower)/cos(alpha2[bin1][bin2][index_alpha]);

          for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            r_bins[bin1][bin2][index_alpha][index_r] = r_min+index_r*(r_max_per_alpha-r_min)/(pgb2->r_size-1);
          }
        }

      }
      if (ratio == 0.0){
        printf("ERROR ratio of horizontal to vertical side of rectangle cut-out is ZERO. bin = %d with bin = %d\n", bin1, bin2);
        exit(2);
      }

      else{
        int horizontal_size = ceil((pgb2->alpha_size-1)*(1.0-(ratio/(ratio+1.0))));
        for (int index_alpha = 0; index_alpha < horizontal_size; index_alpha++) {
          double r_max_per_alpha = (1.0-1.0/tau2_cutoff_lower)/sin(alpha2[bin1][bin2][index_alpha]);

          for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            r_bins[bin1][bin2][index_alpha][index_r] = r_min+index_r*(r_max_per_alpha-r_min)/(pgb2->r_size-1);
          }
        }

        for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
          r_bins[bin1][bin2][horizontal_size][index_r] = r_min+index_r*(r_abs_max-r_min)/(pgb2->r_size-1);
        }

        for (int index_alpha = horizontal_size; index_alpha < pgb2->alpha_size; index_alpha++) {
          double r_max_per_alpha = (1.0-tau0/tau2_cutoff_lower)/cos(alpha2[bin1][bin2][index_alpha]);

          for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            r_bins[bin1][bin2][index_alpha][index_r] = r_min+index_r*(r_max_per_alpha-r_min)/(pgb2->r_size-1);
          }
        }
      }
    }
  }

  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
      double taux = tau0*(1.0-r_bins[0][0][index_alpha][index_r]*sin(alpha2[0][0][index_alpha]));
      double tauy = tau0*(1.0-r_bins[0][0][index_alpha][index_r]*cos(alpha2[0][0][index_alpha]));
      //printf("alpha2[0][0][%d] = %g\n", index_alpha, alpha2[0][0][index_alpha]);
      //printf("%g    %g\n",taux, tauy);
    }
  }
  // NOTE: the else part doesn't have the increment, the middle alpha index is not correct!
  printf("Completed polar time-parameterisation\n");
  //exit(0);






  printf("tau0 = %g\n",tau0 );
  printf("tau_max_cls = %g\n", tau_max_cls);
  printf("tau_min_cls = %g\n", tau_min_cls);



  r_abs_max = (1.0-(tau_min_cls/tau0))/sin(_PI_/4);

  printf("r_abs_max = %G\n", r_abs_max);


  printf("r_abs_max = %g\n", r_abs_max);
  printf("r_abs_max = %g\n", (1-pgb2->tau_sampling_cls[0]/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1])/cos(_PI_/4));

  double * r_max_per_alpha;

  int middle_alpha_index = (pgb2->alpha_size+1)/2-1;


  class_alloc(r_max_per_alpha,
              pgb2->alpha_size * sizeof(double*),
              pgb2->error_message);


  for (int index_alpha = 0; index_alpha < (pgb2->alpha_size+1)/2; index_alpha++) {
    r_max_per_alpha[index_alpha] = (tau0-pgb2->tau_sampling_cls[0])/(tau0*cos(pgb2->alpha[index_alpha]));
  }

  r_max_per_alpha[middle_alpha_index] = r_abs_max;

  for (int index_alpha = middle_alpha_index+1; index_alpha < pgb2->alpha_size; index_alpha++) {
    int distance_from_middle = index_alpha - middle_alpha_index;
    r_max_per_alpha[index_alpha] = r_max_per_alpha[middle_alpha_index-distance_from_middle];
    //printf("distance_from_middle = %d, set equal to index = %d\n", distance_from_middle, middle_alpha_index-distance_from_middle );
    //printf("r_max_per_alpha[%d] = %g, alpha = %g\n", index_alpha, r_max_per_alpha[index_alpha], pgb2->alpha[index_alpha]);
  }

  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    for (int index_r = 0; index_r < pgb2->r_size; index_r++) {

      pgb2->r2[index_alpha][index_r]= r_min2 + index_r*((r_max_per_alpha[index_alpha]-r_min2)/(pgb2->r_size-1));
    }
  }

  free(r_max_per_alpha);


  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++){
    for (int index_r = 0; index_r < pgb2->r_size; index_r++) {

      tau_one = tau0*(1.0-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
      tau_two = tau0*(1.0-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));

      //printf("%g    %g    %g\n", tau_one, tau_two, pgb2->r2[index_alpha][index_r]);
    }
    /*printf("%g    %g\n",
          pba->conformal_age*(1.0-pgb2->r2[index_alpha][70]*sin(pgb2->alpha[index_alpha])),
          pba->conformal_age*(1.0-pgb2->r2[index_alpha][70]*cos(pgb2->alpha[index_alpha])));*/
  }
  //herehere
  double * w_trapz_r;
  double ** w_trapz_r2;

  /* Allocate and fill array for the trapezoidal weights for Chi integration w_trapz[bin][index_tau] */


  class_alloc(w_trapz_r,
              pgb2->r_size * sizeof(double),
              ppt->error_message);

  class_alloc(w_trapz_r2,
              pgb2->alpha_size * sizeof(double*),
              ppt->error_message);

  for ( int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    class_alloc(w_trapz_r2[index_alpha],
                pgb2->r_size * sizeof(double),
                pgb2->error_message);


    class_call(array_trapezoidal_weights(pgb2->r2[index_alpha],
                                         pgb2->r_size,
                                         w_trapz_r2[index_alpha],
                                         pgb2->error_message),
                                         pgb2->error_message,
                                         pgb2->error_message);
  }






  /* Experimenting with the Simpson's rule */



  double ** simps_weights;

  class_alloc(simps_weights,
              ppt->selection_num * sizeof(double*),
              pgb2->error_message);

  /* Allocation of second dimension selection[bin][index_tau] */
  for(int bin = 0; bin < ppt->selection_num; bin++){

    class_alloc(simps_weights[bin],
                pgb2->tau_size_selection * sizeof(double),
                ppt->error_message);

    compute_simps_weights_selection(bin,
                    simps_weights[bin],
                    pgb2);
                    // ALERT simps_weights used
    /*compute_simps_weights_selection(bin,
                    w_trapz[bin],
                    pgb2);*/


  }






  /*int index_test_l2 = 2;
  double j_start;
  double j_end;
  double j_even;
  double j_odd;
  double j_simps;
  printf("ptr->l[%d] = %d\n", index_test_l2, ptr->l[index_test_l2]);
  //double tau0 = pba->conformal_age;
  class_call(bessel_at_x(pbs, tau0 - pgb2->tau_sampling_selection[0][0] , index_test_l2, &j_start), pbs->error_message, pgb2->error_message);
  class_call(bessel_at_x(pbs, tau0 - pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1] , index_test_l2, &j_end), pbs->error_message, pgb2->error_message);


  double delta_x = (pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1]-pgb2->tau_sampling_selection[0][0])/(pgb2->tau_size_selection-1);
  //double delta_x = (j_end-j_start)/(pgb2->tau_size_selection-1);
  double simpsons_sum = 0.0;
  double trapz_sum = 0.0;
  double simps_weights_sum = 0.0;
  //simpsons_sum += delta_x*(pow(pgb2->tau_sampling_selection[0][0],4)+pow(pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1],4))/3.0;
  simpsons_sum += delta_x*(j_start*j_start+j_end*j_end)/3.0;
  //count even sums
  for (int i = 1; i < (pgb2->tau_size_selection-3)/2+1; i++) {
    class_call(bessel_at_x(pbs, tau0 - pgb2->tau_sampling_selection[0][2*i] , index_test_l2, &j_even), pbs->error_message, pgb2->error_message);
    //simpsons_sum += 2.0*delta_x*(pow(pgb2->tau_sampling_selection[0][2*i],4))/3.0;
    simpsons_sum += 2.0*delta_x*(j_even*j_even)/3.0;
    //printf("even weight = %g\n", 2.0*delta_x/3.0);
  }
  // count odd sums
  for (int i = 1; i < (pgb2->tau_size_selection-3)/2+2; i++) {
    //simpsons_sum += 4.0*delta_x*(pow(pgb2->tau_sampling_selection[0][2*i-1],4))/3.0;
    class_call(bessel_at_x(pbs, tau0 - pgb2->tau_sampling_selection[0][2*i-1] , index_test_l2, &j_odd), pbs->error_message, pgb2->error_message);
    simpsons_sum += 4.0*delta_x*(j_odd*j_odd)/3.0;
    //printf("odd weight = %g\n", 4.0*delta_x/3.0);
  }
  // trapz_sum
  for (int i = 0; i < pgb2->tau_size_selection; i++) {
    //trapz_sum += pow(pgb2->tau_sampling_selection[0][i],4)*w_trapz[0][i];
    class_call(bessel_at_x(pbs, tau0 - pgb2->tau_sampling_selection[0][i] , index_test_l2, &j_simps), pbs->error_message, pgb2->error_message);
    trapz_sum += j_simps*j_simps*w_trapz[0][i];
    simps_weights_sum += j_simps*j_simps*simps_weights[0][i];
  }

  printf("tau_sampling_selection[0][0] = %g, tau_sampling_selection[0][301] = %g\n", pgb2->tau_sampling_selection[0][0], pgb2->tau_sampling_selection[0][300]);
  printf("simpsons_sum = %g\n",simpsons_sum);
  printf("trapz_sum = %g\n",trapz_sum);
  printf("simps_weights_sum = %g\n", simps_weights_sum);*/




  double tau_window;
  int last_index_window;
  last_index_window = 0;
  double * bac_window;
  double normal;
  double z;
  double width;
  double mean;
  double x;
  width = ppt->selection_width[0];
  mean = ppt->selection_mean[0];
  double test_sum;
  test_sum = 0.0;
  class_alloc(bac_window, pba->bg_size * sizeof(double), pba->error_message);





  class_call(background_tau_of_z(
                          pba,
                          ppt->selection_mean[0],
                          &selection_mean),
                          ppt->error_message,
                          pgb2->error_message);


  //free(tau0_minus_tau);
  /* Allocate and fill pgb2->vecback, the array which is filled with background information. This will be used throughout
    the module. It is faster to do this now and call the same array each time, rather than recompute it repeatedly. */
  printf("Computing relevant background information...\n");
  class_alloc(pgb2->vecback, pba->bg_size * sizeof(double), pba->error_message);
  int last_index_bg;

  for (int index_tau_lens = 0; index_tau_lens < pgb2->tau_size_bessel; index_tau_lens++) {

    class_call(background_at_tau(pba,
                                 pgb2->tau_sampling_bessel[index_tau_lens],
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_bg,
                                 pgb2->vecback),
                                 pba->error_message,
                                 pgb2->error_message);
  }
  printf("Background information complete.\n");


  double t1,t2,k1;
  int index_k1, index_t1, index_t2;
  int index_k_bessel;

  printf("ppt->tau_sampling_quadsources has %d points sampled between (%g,%g)\n", ppt->tau_size_quadsources, ppt->tau_sampling_quadsources[0], ppt->tau_sampling_quadsources[ppt->tau_size_quadsources-1] );
  printf("pgb2->tau_sampling_cls has %d points sampled between (%g,%g)\n", pgb2->tau_size_cls, pgb2->tau_sampling_cls[0], pgb2->tau_sampling_cls[pgb2->tau_size_cls-1] );
  printf("pgb2->tau_sampling_selection has %d points sampled between (%g,%g)\n", pgb2->tau_size_selection, pgb2->tau_sampling_selection[0], pgb2->tau_sampling_selection[pgb2->tau_size_selection-1]);

  int dump = 0;
  double f,g;
  int i2;
  int index;
  int index_source;
  double tau;
  int last_index;
  int last_index_k;

  //NOTE: Fix this
  index_type = 0;
  index_source = 0;
  double k5 = 5.0;

/* NOTE: The following if statements should be dependent on the user input */
// if want delta_cdm

/*------------------------------------------------
================ Type Declaration ================
-------------------------------------------------*/

  /* Here we initialise and define the different 'types' which are the various contributions to galaxy number counts at first and second order.
    We initialise each type to a negative integer. It is easier this way to skip cointributions that are not required: loops over type start from
    index_type = 0 and will therefore avoid any of the types which are not requested. TODO: The turning on and off of different contributions should
    be the requested by the user in the input file. At present it is manual in the following lines.*/

  pgb2->index_type_density = -1;
  pgb2->index_type_rsd = -1;
  pgb2->index_type_lens = -1;
  pgb2->index_type_d1 = -1;
  pgb2->index_type_d2 = -1;
  pgb2->index_type_g1 = -1;
  pgb2->index_type_g2 = -1;
  pgb2->index_type_g3 = -1;
  pgb2->index_type_g4 = -1;
  pgb2->index_type_g5 = -1;

  pgb2->index_source_v = -1;
  pgb2->index_source_theta = -1;
  pgb2->index_source_delta_cdm = -1;
  pgb2->index_source_phi_plus_psi = -1;
  pgb2->index_source_phi_plus_psi_prime = -1;
  pgb2->index_source_phi = -1;
  pgb2->index_source_psi = -1;
  pgb2->index_source_phi_prime = -1;


  /* Next we want to turn on any source types (needed for interpolation) and make up the contribution types
    needed by the integral function. For example pgb2->index_source_theta is needed for interpolation but not for
    the integral (over k) function, therefore it is set to something other than -1 but is not counted in the pgb2->type_size.
    This will speed up the code since it will not need to be checked in the integral function. */

  // turn on for delta_cdm
  /*if (k5 == 5.0){
    pgb2->index_source_delta_cdm = index_source;
    index_source++;
    pgb2->index_type_density = index_type;
    index_type++;
    //pgb2->index_source_phi = index_source;
    //index_source++;
    /*pgb2->index_source_v = index_type;
    index_source++;
    pgb2->index_source_phi = index_type;
    index_source++;*/
  //}
  // turn on for rsd

  if (k5 == 5.0){
    //pgb2->index_source_v = index_source;
    //index_source++;
    pgb2->index_source_theta = index_source;
    index_source++;

    pgb2->index_type_rsd = index_type;
    index_type++;
    //pgb2->index_type_d1 = index_type;
    //index_type++;
    //pgb2->index_type_d2 = index_type;
    //index_type++;
  }

  /*if (k5 == 5.0){
    pgb2->index_source_psi = index_source;
    index_source++;
    pgb2->index_source_phi = index_source;
    index_source++;
    pgb2->index_source_phi_prime = index_source;
    index_source++;
    pgb2->index_source_phi_plus_psi = index_source;
    index_source++;
    pgb2->index_source_phi_plus_psi_prime = index_source;
    index_source++;*/


    /*pgb2->index_type_g1 = index_type;
    index_type++;
    pgb2->index_type_g2 = index_type;
    index_type++;
    pgb2->index_type_g3 = index_type;
    index_type++;
    pgb2->index_type_g4 = index_type;
    index_type++;
    pgb2->index_type_g5 = index_type;
    index_type++;*/
  //}



  /*if (k5 == 5.0){
    pgb2->index_type_lens = index_type;
    index_type++;
  }*/

  pgb2->source_size = index_source;
  pgb2->type_size = index_type;



  printf("type size = %d\n", pgb2->type_size );
  printf("source size = %d\n", pgb2->source_size );
  printf("tau_size_selection = %d\n", pgb2->tau_size_selection);
  printf("k_size_bessel = %d\n",pgb2->k_size_bessel);

  /* Define an array of values of first order transfer functions:
              pgb2->first_order_sources[index_type][index_tau][index_k_bessel] */
  class_alloc(pgb2->first_order_sources, pgb2->source_size * sizeof(double **), ppt->error_message);
    for (int index_type = 0; index_type < pgb2->source_size; index_type++) {
      /* Allocate memory for pgb2->first_order_sources[index_type][index_tau] */
      class_alloc(pgb2->first_order_sources[index_type],
                  pgb2->tau_size_cls * sizeof(double *),
                  ppt->error_message);
      /* Allocate memory for pgb2->first_order_sources[index_type] */
      for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {
          /* Loop over type and tau. For each of them, allocate memory
           for pgb2->first_order_sources[index_type][index_tau][index_k_bessel]  */
        class_alloc(pgb2->first_order_sources[index_type][index_tau],
                    pgb2->k_size_bessel * sizeof(double),
                    ppt->error_message);
    }
  }

  printf("First order sources allocated\n" );

  class_alloc(pgb2->first_order_sources_integrand, pgb2->source_size * sizeof(double **), ppt->error_message);
    for (int index_type = 0; index_type < pgb2->source_size; index_type++) {
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

  printf("First order sources_integrand allocated\n" );

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
                    pgb2->tau_size_cls * sizeof(double *),
                    ppt->error_message);
      /* Allocate memory for pgb2->first_order_sources_integ[index_type] */
        for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {

            /* Loop over type, l and tau. For each of them, allocate memory
             for pgb2->first_order_sources_integ[index_type][index_l][index_tau][index_k_bessel]  */
           class_alloc(pgb2->first_order_sources_integ[index_type][index_l][index_tau],
                        pgb2->k_size_bessel * sizeof(double),
                        ppt->error_message);
        }
      }
    }

    printf("First order sources_integ allocated\n" );



  //printf("%g\n",ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][100 * ppt->k_size[ppt->index_md_scalars] + 20]);


  /***************************************************************************
  ==============   Interpolate on first_order_sources array =================
  ****************************************************************************/

/* Now fill and interpolate this new pointer-array pgb2->first_order_sources[index_type][index_tau][index_k]
  with information from the pre computed ppt->quadsources[index_md][index_ic*ppt->tp_size[index_md]+index_type]
  [index_tau * ppt->k_size[index_md] + index_k].The galbispectra2 module uses finer k and tau sampling at late times
  as compared with the pt2 module. We therefore have to interpolate any of the transfer functions we wish to use.*/



    printf("Source array allocation completed.\n" );

  /* NOTE: The following first_order_sources interpolation is only necessary for each of the perturbations delta, v, phi etc. ,
    terms in the galaxy number over-density which share these terms, can just reuse the associated */
    double cdm, b, pho;
    double intermediate;
    double intermediate_plus;
    double ** matter;

    /*class_alloc(matter, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
      for (int index = 0; index < ppt->tau_size_quadsources; index++) {
        /* Allocate memory for pgb2->first_order_sources[index_type][index_tau] */
        /*class_alloc(matter[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);
    }

    for (int index = 0; index < ppt->tau_size_quadsources; index++) {
      for (int index_k = 0; index_k < ppt->k_size[ppt->index_md_scalars]; index_k++) {
        matter[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k]
          + ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_b][index * ppt->k_size[ppt->index_md_scalars] + index_k];
      }
    }*/
    //pgb2->first_order_sources[0][0][0] = 0.0;
    if (pgb2->index_source_delta_cdm != -1) {
      printf("Preparing density source term..\n");
      for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
        index_k = 0;

        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
          index = 0;

          double tau = pgb2->tau_sampling_cls[index_tau];

          class_call(index_of_tau_sampling_quadsources(tau, &index, ppt),ppt->error_message,pgb2->error_message);

          double k = pgb2->k_bessel[index_k_bessel];

          class_call(index_of_k(k, &index_k, ppt), ppt->error_message, pgb2->error_message);

          f = (pgb2->tau_sampling_cls[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

          intermediate  = (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k]+
              (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k]);

          intermediate_plus =  (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k+1]+
              (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k+1]);

          g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

          pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_bessel] = g*intermediate_plus +(1-g)*intermediate;


        }
      }
      double j1;
      for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
        for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
          for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
            double x1 = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);
            class_call(bessel_at_x(pbs, x1 , index_l, &j1), pbs->error_message, pgb2->error_message);

            pgb2->first_order_sources_integ[pgb2->index_source_delta_cdm][index_l][index_tau][index_k_bessel] =
              pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_bessel] * j1;
          }
        }
      }
    }


    // set up sources for rsd


    intermediate = 0;
    intermediate_plus = 0;
    int last_index_rsd;
    int last_index_k_rsd;
    int index_k_rsd;
    int index_test = 4;
    int index_test_qs=0;

    class_call(index_of_tau_sampling_quadsources(pgb2->tau_sampling_cls[index_test], &index_test_qs, ppt), pgb2->error_message, pgb2->error_message);

    printf("ppt->tau_sampling_quadsources[%d] = %g\n",index_test_qs, ppt->tau_sampling_quadsources[index_test_qs]);

    if (pgb2->index_source_v != -1) {
      for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
        // NOTE: ppt->index_qs_theta_cdm is equal to -ppt->index_qs_v_cdm/k*k velocity. v = -theta/(k*k)
        index_k_rsd = 0;

        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
          //Coarse time index
          index = 0;
          //Coarse tau value
          double tau = pgb2->tau_sampling_cls[index_tau];

          class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);

          double k = pgb2->k_bessel[index_k_bessel];

          class_call(index_of_k(k, &index_k_rsd, ppt), pgb2->error_message, pgb2->error_message);

          f = (pgb2->tau_sampling_cls[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

          intermediate  = (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]+
              (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]);

          intermediate_plus =  (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k_rsd+1]+
              (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k_rsd+1]);

          g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k_rsd])/(ppt->k[ppt->index_md_scalars][index_k_rsd+1]-ppt->k[ppt->index_md_scalars][index_k_rsd]);
          //printf("%d    %d    %d  %g\n", pgb2->index_source_v, index_tau, index_k_bessel, pgb2->first_order_sources[pgb2->index_source_v][index_tau][index_k_bessel]);
          pgb2->first_order_sources[pgb2->index_source_v][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);

          /* Debug: Check the interpolation is working correctly. */
          //printf("%g  %g  %g  %g\n",
            //pgb2->k_bessel[index_k_bessel],
            //pgb2->first_order_sources[pgb2->index_type_rsd][index_test][index_k_bessel],
            //ppt->k[ppt->index_md_scalars][index_k_rsd],
            //ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][index_test_qs * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]);
          //printf(" \n");
          //printf("ppt->quadsources[ppt->index_md_scalars][ppt->index_qs_theta_cdm][%d * ppt->k_size[ppt->index_md_scalars] + k[%d]=%g) = %g\n",index, index_k_rsd, ppt->k[ppt->index_md_scalars][index_k_rsd], ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_theta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]);
          //printf("ppt->quadsources[ppt->index_md_scalars][ppt->index_qs_v_cdm][%d * ppt->k_size[ppt->index_md_scalars] + k[%d]=%g) = %g\n",index, index_k_rsd, ppt->k[ppt->index_md_scalars][index_k_rsd],ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]);
          //printf("ppt->quadsources[ppt->index_md_scalars][ppt->index_qs_delta_cdm][%d * ppt->k_size[ppt->index_md_scalars] + k[%d]=%g) = %g\n",index, index_k_rsd, ppt->k[ppt->index_md_scalars][index_k_rsd], ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]);

        }
      }
    }

    if (pgb2->index_source_theta != -1) {
      printf("Preparing RSD source term..\n");
      for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
        // NOTE: ppt->index_qs_theta_cdm is equal to -ppt->index_qs_v_cdm/k*k velocity. v = -theta/(k*k)
        index_k_rsd = 0;

        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
          index = 0;

          double tau = pgb2->tau_sampling_cls[index_tau];

          class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);

          double k = pgb2->k_bessel[index_k_bessel];

          class_call(index_of_k(k, &index_k_rsd, ppt), pgb2->error_message, pgb2->error_message);

          f = (pgb2->tau_sampling_cls[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

          intermediate  = (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_theta_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]+
              (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_theta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]);

          intermediate_plus =  (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_theta_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k_rsd+1]+
              (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_theta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k_rsd+1]);

          g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k_rsd])/(ppt->k[ppt->index_md_scalars][index_k_rsd+1]-ppt->k[ppt->index_md_scalars][index_k_rsd]);

          pgb2->first_order_sources[pgb2->index_source_theta][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);


          /* Debug: Check the interpolation is working correctly. */
          //printf("%g  %g  %g  %g\n",
            //pgb2->k_bessel[index_k_bessel],
            //pgb2->first_order_sources[pgb2->index_type_rsd][index_test][index_k_bessel],
            //ppt->k[ppt->index_md_scalars][index_k_rsd],
            //ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][index_test_qs * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]);
          //printf(" \n");
          //printf("ppt->quadsources[ppt->index_md_scalars][ppt->index_qs_theta_cdm][%d * ppt->k_size[ppt->index_md_scalars] + k[%d]=%g) = %g\n",index, index_k_rsd, ppt->k[ppt->index_md_scalars][index_k_rsd], ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_theta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]);
          //printf("ppt->quadsources[ppt->index_md_scalars][ppt->index_qs_v_cdm][%d * ppt->k_size[ppt->index_md_scalars] + k[%d]=%g) = %g\n",index, index_k_rsd, ppt->k[ppt->index_md_scalars][index_k_rsd],ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_v_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]);
          //printf("ppt->quadsources[ppt->index_md_scalars][ppt->index_qs_delta_cdm][%d * ppt->k_size[ppt->index_md_scalars] + k[%d]=%g) = %g\n",index, index_k_rsd, ppt->k[ppt->index_md_scalars][index_k_rsd], ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k_rsd]);

        }
      }
    }


    // ALERT MAJOR CHANGES BELOW //

      /********************************************/
      /*========= Source Term Preparation ========*/
      /*__________________________________________*/

      double * pvecback11;
      double * pvecback22;
      int last_index_rsd2 = 0;
      double lensing_result = 0.;
      double j, j_first_deriv, j_second_deriv;
      double f_evo1, f_evo2;
      double * pvecback_rsd;
      class_alloc(pvecback_rsd, pba->bg_size * sizeof(double), pba->error_message);


      printf("Entering source term preparation.\n");
      /* /* Here we prepare the source term, ready for the first-integration (over k). We pair input type with corresponding source type. Each type is made up of,
        the transfer function of the perturbation, a prefactor, some k-factor, and a Bessel function (or a derivative of); it is written into the
        pgb2->first_order_sources_integ[pgb2->index_type_density][index_l][index_tau][index_k_bessel]
        array. */

      int index_tau_lens;
      double prefactor1, prefactor2;
      double lensing_result1, lensing_result2;
      double velocity1, velocity2;

        /* Density */
        if(pgb2->index_type_density != -1){
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);

                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);

                pgb2->first_order_sources_integ[pgb2->index_type_density][index_l][index_tau][index_k_bessel] = pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_bessel] * j;
              }
            }
          }
        }

        /* RSD */
        if (pgb2->index_type_rsd != -1) {
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              //printf("indexing index_l x index_tau = %dx%d\n", index_l, index_tau);
              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_rsd,
                                           pvecback_rsd),
                                           pba->error_message,
                                           pgb2->error_message);

              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);

                class_call(bessel_at_x_second_deriv(pgb2, pbs, x, index_l, &j_second_deriv), pbs->error_message, pgb2->error_message);

                double prefactor_rsd = 1.0
                             /pvecback_rsd[pba->index_bg_H]
                             /pvecback_rsd[pba->index_bg_a];

              /* Write in to the array */
                pgb2->first_order_sources_integ[pgb2->index_type_rsd][index_l][index_tau][index_k_bessel] = prefactor_rsd
                                *pgb2->first_order_sources[pgb2->index_source_theta][index_tau][index_k_bessel]
                                *j_second_deriv;
              }
            }
          }
        }
        free(pvecback_rsd);

        int last_index_d1;
        double * pvecback_d1;
        class_alloc(pvecback_d1, pba->bg_size * sizeof(double), pba->error_message);

        /* First Type: Doppler1 */
        if(pgb2->index_type_d1 != -1 ){
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_d1,
                                           pvecback_d1),
                                           pba->error_message,
                                           pgb2->error_message);

              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);
                class_call(bessel_at_x_first_deriv(pgb2, pbs, x, index_l, &j_first_deriv), pbs->error_message, pgb2->error_message);

                f_evo1 = 2.
                         /pvecback_d1[pba->index_bg_H]
                         /pvecback_d1[pba->index_bg_a]
                         /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau])
                         +pvecback_d1[pba->index_bg_H_prime]
                         /pvecback_d1[pba->index_bg_H]
                         /pvecback_d1[pba->index_bg_H]
                         /pvecback_d1[pba->index_bg_a];

              double  prefactor_d1 = (1.
                              +pvecback_d1[pba->index_bg_H_prime]
                              /pvecback_d1[pba->index_bg_a]
                              /pvecback_d1[pba->index_bg_H]
                              /pvecback_d1[pba->index_bg_H]
                              +(2.-5.*ptr->s_bias)
                              /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau])
                              /pvecback_d1[pba->index_bg_a]
                              /pvecback_d1[pba->index_bg_H]
                              +5.*ptr->s_bias
                              -f_evo1);

                pgb2->first_order_sources_integ[pgb2->index_type_d1][index_l][index_tau][index_k_bessel] = prefactor_d1
                        *pgb2->first_order_sources[pgb2->index_source_theta][index_tau][index_k_bessel]
                        *j_first_deriv
                        /pgb2->k_bessel[index_k_bessel];
              }
            }
          }
        }
        double * pvecback_d2;
        class_alloc(pvecback_d2, pba->bg_size * sizeof(double), pba->error_message);
        /* First Type: Doppler2 */
        int last_index_d2;
        if(pgb2->index_type_d2 != -1){
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){

              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_d2,
                                           pvecback_d2),
                                           pba->error_message,
                                           pgb2->error_message);

              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);
                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);

                f_evo1 = 2.
                         /pvecback_d2[pba->index_bg_H]
                         /pvecback_d2[pba->index_bg_a]
                         /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau])
                         +pvecback_d2[pba->index_bg_H_prime]
                         /pvecback_d2[pba->index_bg_H]
                         /pvecback_d2[pba->index_bg_H]
                         /pvecback_d2[pba->index_bg_a];
                         //alert f_evo1 skipped
                //prefactor1 = -3.0*pvecback_d2[pba->index_bg_a]*pvecback_d2[pba->index_bg_H];

                pgb2->first_order_sources_integ[pgb2->index_type_d2][index_l][index_tau][index_k_bessel] = (f_evo1-3.0)
                        *pvecback_d2[pba->index_bg_a]
                        *pvecback_d2[pba->index_bg_H]
                        *pgb2->first_order_sources[pgb2->index_source_theta][index_tau][index_k_bessel]
                        *j
                        /pgb2->k_bessel[index_k_bessel]
                        /pgb2->k_bessel[index_k_bessel];
              }
            }
          }
        }
        free(pvecback_d2);
        /* First Type: g1 (first of the GR terms) */
        if(pgb2->index_type_g1 != -1){
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);

                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);


                double prefactor_g1 = -2.0+5.0*ptr->s_bias;

                pgb2->first_order_sources_integ[pgb2->index_type_g1][index_l][index_tau][index_k_bessel] = prefactor_g1
                        *pgb2->first_order_sources[pgb2->index_source_phi][index_tau][index_k_bessel]
                        *j;
              }
            }
          }
        }

        /* First Type: g2 */
        double * pvecback_g2;
        class_alloc(pvecback_g2, pba->bg_size * sizeof(double), pba->error_message);
        if(pgb2->index_type_g2 != -1 ){
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);
                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);

                f_evo1 = 2.
                         /pvecback_g2[pba->index_bg_H]
                         /pvecback_g2[pba->index_bg_a]
                         /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau])
                         +pvecback_g2[pba->index_bg_H_prime]
                         /pvecback_g2[pba->index_bg_H]
                         /pvecback_g2[pba->index_bg_H]
                         /pvecback_g2[pba->index_bg_a];

                double prefactor_g2 = (2.0
                             +pvecback_g2[pba->index_bg_H_prime]
                             /pvecback_g2[pba->index_bg_H]
                             /pvecback_g2[pba->index_bg_H]
                             /pvecback_g2[pba->index_bg_a]
                             +(2.0-5.0*ptr->s_bias)
                             /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau])
                             /pvecback_g2[pba->index_bg_H]
                             /pvecback_g2[pba->index_bg_a]
                             +5*ptr->s_bias
                             -f_evo1);

                pgb2->first_order_sources_integ[pgb2->index_type_g2][index_l][index_tau][index_k_bessel] = prefactor_g2
                        *pgb2->first_order_sources[pgb2->index_source_psi][index_tau][index_k_bessel]
                        *j;
              }
            }
          }
        }
        free(pvecback_g2);
        /* First Type: g3 */
        double * pvecback_g3;
        class_alloc(pvecback_g2, pba->bg_size * sizeof(double), pba->error_message);
        if(pgb2->index_type_g3 != -1 ){
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);
                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);


                prefactor1 = 1.0
                             /pvecback_g3[pba->index_bg_H]
                             /pvecback_g3[pba->index_bg_a];

                pgb2->first_order_sources_integ[pgb2->index_type_g3][index_l][index_tau][index_k_bessel] = prefactor1
                        *pgb2->first_order_sources[pgb2->index_source_phi_prime][index_tau][index_k_bessel]
                        *j;
              }
            }
          }
        }
        free(pvecback_g3);
    //herehere



    /* Lensing source term preparation */

    intermediate = 0;
    intermediate_plus = 0;
    int last_index_lens;
    int last_index_k_lens;
    int index_k_lens;
    double ** phi_plus_psi;
    double ** phi;
    double ** psi;
    /* Allocate a temporary array phi_plus_psi[index_tau][index_k] that stores the quadsources transfer function for phi+psi
        purely for brevity*/

    class_alloc(phi_plus_psi, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
      for (int index = 0; index < ppt->tau_size_quadsources; index++) {

        class_alloc(phi_plus_psi[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);
      }

    class_alloc(phi, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
    class_alloc(psi, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
      for (int index = 0; index < ppt->tau_size_quadsources; index++) {
        /* Allocate memory for pgb2->first_order_sources[index_type][index_tau] */
        class_alloc(phi[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);

        class_alloc(psi[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);
      }

    for (int index = 0; index < ppt->tau_size_quadsources; index++) {
      for (int index_k = 0; index_k < ppt->k_size[ppt->index_md_scalars]; index_k++) {
        phi_plus_psi[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_phi][index * ppt->k_size[ppt->index_md_scalars] + index_k]
          + ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_psi][index * ppt->k_size[ppt->index_md_scalars] + index_k];

        phi[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_phi][index * ppt->k_size[ppt->index_md_scalars] + index_k];

        psi[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_psi][index * ppt->k_size[ppt->index_md_scalars] + index_k];
      }
    }


    printf("starting pgb2->index_type_lens\n" );

    if (pgb2->index_source_phi_plus_psi != -1) {
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

            f = (pgb2->tau_sampling_bessel[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

            intermediate  = f*phi_plus_psi[index+1][index_k]+(1-f)*phi_plus_psi[index][index_k];

            intermediate_plus =  f*phi_plus_psi[index+1][index_k+1]+(1-f)*phi_plus_psi[index][index_k+1];

            g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

            pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);
          }
      }
    }

    if (pgb2->index_source_phi != -1) {
      for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){

          last_index_k = 0;
          index_k_lens = 0;
          int last_index;
          int index_tau_lens;
          double lensing_result;
          f = 0.;
          g = 0.;

          for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
            //printf("index_tau = %d, index_k_bessel = %d\n",index_tau, index_k_bessel);
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

            f = (pgb2->tau_sampling_bessel[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

            intermediate  = f*phi[index+1][index_k]+(1-f)*phi[index][index_k];

            intermediate_plus =  f*phi[index+1][index_k+1]+(1-f)*phi[index][index_k+1];

            g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

            pgb2->first_order_sources[pgb2->index_source_phi][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);
          }
      }
    }

    /* Interpolate Psi */
    if (pgb2->index_source_psi != -1) {
      for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){

          last_index_k = 0;
          index_k_lens = 0;
          int last_index;
          int index_tau_lens;
          double lensing_result;
          f = 0.;
          g = 0.;

          for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
            index = 0;

            double tau = pgb2->tau_sampling_bessel[index_tau];

            class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);

            double k = pgb2->k_bessel[index_k_bessel];

            class_call(index_of_k_old(k,
                           &index_k,
                           &last_index_k,
                           ppt),
                           pgb2->error_message,
                           pgb2->error_message);

            f = (pgb2->tau_sampling_bessel[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

            intermediate  = f*psi[index+1][index_k]+(1-f)*psi[index][index_k];

            intermediate_plus =  f*psi[index+1][index_k+1]+(1-f)*psi[index][index_k+1];

            g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

            pgb2->first_order_sources[pgb2->index_source_psi][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);
          }
      }
    }

    /* Now using the interpolated phi transfer function, we will take the time derivatives to fill the pgb2->index_source_phi_prime
      transfer array, we take the derivatives of the end points and internal points separately. */
    if (pgb2->index_source_phi_prime != -1 || pgb2->index_source_phi_plus_psi_prime != -1) {
      for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
        pgb2->first_order_sources[pgb2->index_source_phi_prime][0][index_k_bessel] = (pgb2->first_order_sources[pgb2->index_source_phi][1][index_k_bessel]-pgb2->first_order_sources[pgb2->index_source_phi][0][index_k_bessel])/
          (pgb2->tau_sampling_cls[1]-pgb2->tau_sampling_cls[0]);

        pgb2->first_order_sources[pgb2->index_source_phi_prime][pgb2->tau_size_cls-1][index_k_bessel] = (pgb2->first_order_sources[pgb2->index_source_phi][pgb2->tau_size_cls-1][index_k_bessel]-pgb2->first_order_sources[pgb2->index_source_phi][pgb2->tau_size_cls-2][index_k_bessel])/
          (pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]-pgb2->tau_sampling_cls[pgb2->tau_size_cls-2]);


        pgb2->first_order_sources[pgb2->index_source_phi_plus_psi_prime][0][index_k_bessel] = (pgb2->first_order_sources[pgb2->index_source_phi_plus_psi][1][index_k_bessel]-pgb2->first_order_sources[pgb2->index_source_phi_plus_psi][0][index_k_bessel])/
          (pgb2->tau_sampling_cls[1]-pgb2->tau_sampling_cls[0]);

        pgb2->first_order_sources[pgb2->index_source_phi_plus_psi_prime][pgb2->tau_size_cls-1][index_k_bessel] = (pgb2->first_order_sources[pgb2->index_source_phi_plus_psi][pgb2->tau_size_cls-1][index_k_bessel]-pgb2->first_order_sources[pgb2->index_source_phi_plus_psi][pgb2->tau_size_cls-2][index_k_bessel])/
          (pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]-pgb2->tau_sampling_cls[pgb2->tau_size_cls-2]);
      }


      for (int index_tau = 1; index_tau < pgb2->tau_size_cls-1; index_tau++){
        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
          pgb2->first_order_sources[pgb2->index_source_phi_prime][index_tau][index_k_bessel] = (pgb2->first_order_sources[pgb2->index_source_phi][index_tau+1][index_k_bessel]-pgb2->first_order_sources[pgb2->index_source_phi][index_tau-1][index_k_bessel])/
            (pgb2->tau_sampling_cls[index_tau+1]-pgb2->tau_sampling_cls[index_tau-1]);

          pgb2->first_order_sources[pgb2->index_source_phi_plus_psi_prime][index_tau][index_k_bessel] = (pgb2->first_order_sources[pgb2->index_source_phi_plus_psi][index_tau+1][index_k_bessel]-pgb2->first_order_sources[pgb2->index_source_phi_plus_psi][index_tau-1][index_k_bessel])/
            (pgb2->tau_sampling_cls[index_tau+1]-pgb2->tau_sampling_cls[index_tau-1]);
        }
      }
    }


    /* Lensing integral, at present this is the slowest part of the code. We check that this computation is required by the first if
      statement. */

    if (pgb2->index_type_lens != -1) {
      for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
        for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){

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

            for (index_tau_lens = (index_tau * bessel_boost); index_tau_lens < pgb2->tau_size_bessel-1; index_tau_lens++) {
              // NOTE that we skip the last index "pgb2->tau_size_bessel-1" to avoid the singularity. This would analytically drop out anyways
              // for l not 1.

              tau_lens = pgb2->tau_sampling_bessel[index_tau_lens];

              if(x>pbs->x_max ){
              printf("ALERT! x= %g \n",x);continue;}

              class_call(bessel_at_x(pbs,pgb2->k_bessel[index_k_bessel]*(pba->conformal_age -tau_lens), index_l, &j), pbs->error_message, pgb2->error_message);

              index_of_tau_sampling_cls(tau_lens, &index_tau_cls, pgb2);

              first_index_tau_in_lens = index_tau * bessel_boost;

              /* If we are taking the first or last trapezoidal weight then take the first weight, else take the other weight */
              if (index_tau_lens == first_index_tau_in_lens || index_tau_lens == (pgb2->tau_size_bessel-1) ){
                weight = (pgb2->tau_sampling_bessel[first_index_tau_in_lens]-pgb2->tau_sampling_bessel[pgb2->tau_size_bessel-1])/(2.0*((pgb2->tau_size_bessel-1)-first_index_tau_in_lens));
              }

              else{
                weight = (pgb2->tau_sampling_bessel[first_index_tau_in_lens]-pgb2->tau_sampling_bessel[pgb2->tau_size_bessel-1])/((pgb2->tau_size_bessel-1)-first_index_tau_in_lens);
              }

              lensing_result += ((pba->conformal_age - pgb2->tau_sampling_bessel[index_tau_lens])-(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]))
                * j * weight * pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][index_tau_lens][index_k_bessel]/((pba->conformal_age - pgb2->tau_sampling_bessel[index_tau_lens])*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]));


            }

            pgb2->first_order_sources_integ[pgb2->index_type_lens][index_l][index_tau][index_k_bessel] = lensing_result;

          }
        }
      }
    }
    double * pvecbackg;
    int last_index_g = 0;
    class_alloc(pvecbackg, pba->bg_size*sizeof(double), pba->error_message);

    /* Prepare the two integrated GR terms (g4 and g5) for the integral function */
    if ((pgb2->index_type_g5 != -1) || (pgb2->index_type_g4 != -1)) {
      printf("Integrating g4 and g5\n");

      for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
        for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){

        last_index_k_lens = 0;
        index_k_lens = 0;

        int first_index_tau_in_lens;
        double f_minus, f_plus, w_minus, w_plus, f_interpolate, tau_lens;
        double weight, end_weights;
        double x;
        double g5;
        double g4;
        double g5_prefactor;
        double f_evo;
        double j;

        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
            g5 = 0.;
            g4 = 0.;
            first_index_tau_in_lens = index_tau * bessel_boost;

            for (int index_tau_lens = (index_tau * bessel_boost); index_tau_lens < pgb2->tau_size_bessel-1; index_tau_lens++) {
              // NOTE that we skip the last index "pgb2->tau_size_bessel-1" to avoid the singularity. This would analytically drop out anyways
              // for l not 1.

              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_bessel[index_tau_lens],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_g,
                                           pvecbackg),
                                           pba->error_message,
                                           pgb2->error_message);

              f_evo = 2.
                      /pvecbackg[pba->index_bg_H]
                      /pvecbackg[pba->index_bg_a]
                      /(pba->conformal_age - pgb2->tau_sampling_bessel[index_tau_lens])
                      +pvecbackg[pba->index_bg_H_prime]
                      /pvecbackg[pba->index_bg_H]
                      /pvecbackg[pba->index_bg_H]
                      /pvecbackg[pba->index_bg_a];

              tau_lens = pgb2->tau_sampling_bessel[index_tau_lens];

              if(x>pbs->x_max ){
              printf("ALERT! x= %g \n",x);continue;}

              class_call(bessel_at_x(pbs,pgb2->k_bessel[index_k_bessel]*(pba->conformal_age -tau_lens), index_l, &j), pbs->error_message, pgb2->error_message);

              //index_of_tau_sampling_cls(tau_lens, &index_tau_cls,  pgb2);

              first_index_tau_in_lens = index_tau * bessel_boost;

              /* If we are taking the first or last trapezoidal weight then take the first weight, else take the other weight */
              if (index_tau_lens == first_index_tau_in_lens || index_tau_lens == (pgb2->tau_size_bessel-1) ){
                weight = (pgb2->tau_sampling_bessel[first_index_tau_in_lens]-pgb2->tau_sampling_bessel[pgb2->tau_size_bessel-1])/(2.0*((pgb2->tau_size_bessel-1)-first_index_tau_in_lens));
              }

              else{
                weight = (pgb2->tau_sampling_bessel[first_index_tau_in_lens]-pgb2->tau_sampling_bessel[pgb2->tau_size_bessel-1])/((pgb2->tau_size_bessel-1)-first_index_tau_in_lens);
              }

              g5_prefactor = (1.0
                            +pvecbackg[pba->index_bg_H_prime]
                            /pvecbackg[pba->index_bg_a]
                            /pvecbackg[pba->index_bg_H]
                            /pvecbackg[pba->index_bg_H]
                            +(2.0-5.0*ptr->s_bias)
                            /(pba->conformal_age - pgb2->tau_sampling_bessel[index_tau_lens])
                            /pvecbackg[pba->index_bg_a]
                            /pvecbackg[pba->index_bg_H]
                            +5.0*ptr->s_bias
                            -f_evo);

              g4 += (2.0-5.0*ptr->s_bias)
                    *pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][index_tau_lens][index_k_bessel]
                    *j
                    *weight
                    /(pba->conformal_age - pgb2->tau_sampling_bessel[index_tau_lens]);

              g5 += g5_prefactor
                    *pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi_prime][index_tau_lens][index_k_bessel]
                    *weight
                    *j;
            }
            pgb2->first_order_sources_integ[pgb2->index_type_g4][index_l][index_tau][index_k_bessel] = g4;
            pgb2->first_order_sources_integ[pgb2->index_type_g5][index_l][index_tau][index_k_bessel] = g5;
          }
        }
      }
    }
  printf("First order source array filled.\n");


  /* Turn off types which are not needed for the integral function */
  /***************************************************************************
  ==========================   k-Integration =================================
  ****************************************************************************/

/*  class_alloc(pgb2->Dl, pgb2->type_size * sizeof(double ****), pgb2->error_message);
  for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    class_alloc(pgb2->Dl[index_type_first], pgb2->type_size * sizeof(double ***), pgb2->error_message);
    for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
      class_alloc(pgb2->Dl[index_type_first][index_type_second], ptr->l_size[ppt->index_md_scalars] * sizeof(double **), pgb2->error_message);
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
        class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l], pgb2->tau_size_cls * sizeof(double *), pgb2->error_message);
        for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_cls; index_tau_first++) {
          class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first], pgb2->tau_size_cls * sizeof(double), pgb2->error_message);
        }
      }
    }
  }*/
  printf("Starting k-integration\n");

  /* Allocate array for Dl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] */
/*  class_alloc(pgb2->Dl, pgb2->type_size * sizeof(double ****), pgb2->error_message);
  for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    class_alloc(pgb2->Dl[index_type_first], pgb2->type_size * sizeof(double ***), pgb2->error_message);
    for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
      class_alloc(pgb2->Dl[index_type_first][index_type_second], ptr->l_size[ppt->index_md_scalars] * sizeof(double **), pgb2->error_message);
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
        class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l], pgb2->tau_size_cls * sizeof(double *), pgb2->error_message);
        for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_cls; index_tau_first++) {
          class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first], pgb2->tau_size_cls * sizeof(double), pgb2->error_message);
        }
      }
    }
  }*/
/* Allocate array for Dl[index_type_first][index_type_second][index_l][index_alpha][index_r] */
  class_alloc(pgb2->Dl2, pgb2->type_size * sizeof(double ****), pgb2->error_message);
  for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    class_alloc(pgb2->Dl2[index_type_first], pgb2->type_size * sizeof(double ***), pgb2->error_message);
    for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
      class_alloc(pgb2->Dl2[index_type_first][index_type_second], ptr->l_size[ppt->index_md_scalars] * sizeof(double **), pgb2->error_message);
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
        class_alloc(pgb2->Dl2[index_type_first][index_type_second][index_l], pgb2->alpha_size * sizeof(double *), pgb2->error_message);
        for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
          class_alloc(pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha], pgb2->r_size * sizeof(double), pgb2->error_message);
        }
      }
    }
  }
  double ***** Dl2_integrand;

  class_alloc(Dl2_integrand, pgb2->type_size * sizeof(double ****), pgb2->error_message);
  for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    class_alloc(Dl2_integrand[index_type_first], pgb2->type_size * sizeof(double ***), pgb2->error_message);
    for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
      class_alloc(Dl2_integrand[index_type_first][index_type_second], ptr->l_size[ppt->index_md_scalars] * sizeof(double **), pgb2->error_message);
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
        class_alloc(Dl2_integrand[index_type_first][index_type_second][index_l], pgb2->alpha_size * sizeof(double *), pgb2->error_message);
        for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
          class_alloc(Dl2_integrand[index_type_first][index_type_second][index_l][index_alpha], pgb2->r_size * sizeof(double), pgb2->error_message);
        }
      }
    }
  }
  printf("Allocating size %ix%ix%ix%ix%i bytes \n", pgb2->type_size, pgb2->type_size, ptr->l_size[ppt->index_md_scalars], pgb2->tau_size_cls, pgb2->tau_size_cls);

  /* Allocate array for Cl[index_type_first][index_type_second][index_l][bin_first][bin_second] */
  class_alloc(pgb2->Cl, pgb2->type_size * sizeof(double ****), pgb2->error_message);
  for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    class_alloc(pgb2->Cl[index_type_first], pgb2->type_size * sizeof(double ***), pgb2->error_message);
    for (int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
      class_alloc(pgb2->Cl[index_type_first][index_type_second], ptr->l_size[ppt->index_md_scalars] * sizeof(double **), ppt->error_message);
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
        class_alloc(pgb2->Cl[index_type_first][index_type_second][index_l], ppt->selection_num * sizeof(double *), ppt->error_message);
        for (int bin_first = 0; bin_first < ppt->selection_num; bin_first++) {
          class_alloc(pgb2->Cl[index_type_first][index_type_second][index_l][bin_first], ppt->selection_num * sizeof(double), ppt->error_message);
        }
      }
    }
  }




  int index1=0;
  int index2=0;

  /*for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
     //for(int index_r = 0; index_r < pgb2->r_size; index_r++){
     int index_r = pgb2->r_size-1;
     double tau1 = pba->conformal_age*(1.0-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
     double tau2 = pba->conformal_age*(1.0-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));

     index_of_tau_sampling_cls(tau1, &index1, pgb2);
     index_of_tau_sampling_cls(tau2, &index2, pgb2);

     printf("%g   %g    %g    %g\n",pgb2->tau_sampling_cls[index1],pgb2->tau_sampling_cls[index2], tau1, tau2);
     //}
   }*/


  double integ, integ_dens_rsd;
  double result2;


  printf("integrating k between %g and %g\n",ppt->k[ppt->index_md_scalars][0],ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1] );

  double p1,p2,f1,f2,j1,j2;
  int k_test_index;
  double pvecback1;
  double pvecback2;

  //TODO: make the background pointers be input parameters in to the integral function
  // ALERT: the second type is counted up to the first type so (1,2) is counted but (2,1) is not.
    // cannot make this short cut when z1 and z2 are not equivalent.


  double type1, type2;
  double Pk;

  FILE *file = fopen("Dl_Data.txt", "w");
if (file == NULL)
{
    printf("Error opening file!\n");
    exit(1);
}

/* print some text */
//const char *text = "Write this to the file";
//fprintf(file, "Some text: %s\n", text);



/* printing single chatacters */
//char c = 'A';
//fprintf(file, "A character: %c\n", c);
  int test_index_l = 5;
  fprintf(file, "#k_size_bessel = %d\n", pgb2->k_size_bessel);
  fprintf(file, "#tau_size_selection = %d\n", pgb2->tau_size_selection);
  fprintf(file, "#Gaussian z1=%g, z2= %g (%g)\n", ppt->selection_mean[0], ppt->selection_mean[1], ppt->selection_width[0] );
  fprintf(file, "#rsd-rsd\n");
  //fprintf(file,"#%s-%s", getName(pgb2->index_type_rsd),getName(pgb2->index_type_rsd));
  fprintf(file, "#tau1      tau2         l(l+1)Dl2(alpha,r)\n");
  fprintf(file, "#l = %d\n", ptr->l[test_index_l] );
  printf("#l = %d\n", ptr->l[test_index_l] );
  double norm_test;
  int index_tau_first;
  int index_tau_second;

  double type1_interp, type2_interp;



  /* Fill the array pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] by calling the integral() function */
  for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    for (int k = 0; k < index_type_first+1; k++){
      printf("X");
    }
    printf("\n");
// alert last tau width is missing
  /* Allocate backgroundmake song vectors */
    double * backvec_z1;
    double * backvec_z2;

    class_alloc(backvec_z1, pba->bg_size * sizeof(double), pba->error_message);
    class_alloc(backvec_z2, pba->bg_size * sizeof(double), pba->error_message);

    int index_tau_first_test = 1000;
    for(int index_type_second = 0; index_type_second < index_type_first+1; index_type_second++){
      //alert we don't run over last index_l
      // ALERT LOOPING ONL UP TO FIRST TIME INDEX
      // ALERT l running
      //ALERT l start
      //for(int index_l = test_index_l; index_l < test_index_l+1/*ptr->l_size[ppt->index_md_scalars]-1*/; index_l++){
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
        printf("index_l = %d\n", index_l );
        //printf("\r          ");
        //printf("%d/%d\n",index_l, ptr->l_size[ppt->index_md_scalars]-1);
        //for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_cls; index_tau_first++){
        for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
          printf("index_alpha = %d\n",index_alpha);
          for(int index_r = 1; index_r < pgb2->r_size; index_r++){
            tau_one = tau0*(1.0-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
            tau_two = tau0*(1.0-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));

            index_of_tau_sampling_cls(tau_one, &index_tau_first, pgb2);
            index_of_tau_sampling_cls(tau_two, &index_tau_second, pgb2);




        //for (int index_tau_first = index_tau_first_test; index_tau_first < index_tau_first_test+1; index_tau_first++){
            class_call(background_at_tau(pba,
                                         pgb2->tau_sampling_cls[index_tau_first],
                                         pba->long_info,
                                         pba->inter_normal,
                                         &last_index,
                                         backvec_z1),
                                         pba->error_message,
                                         pgb2->error_message);

            /* infer redhsift */
            double z1 = pba->a_today/backvec_z1[pba->index_bg_a]-1.;
            //for(int index_tau_second = 0; index_tau_second < index_tau_first+1/*pgb2->tau_size_cls*/; index_tau_second++){
            //for(int index_tau_second = 0; index_tau_second < pgb2->tau_size_cls; index_tau_second++){
            /* Zero the integration sum */
            class_call(background_at_tau(pba,
                                         pgb2->tau_sampling_cls[index_tau_second],
                                         pba->long_info,
                                         pba->inter_normal,
                                         &last_index,
                                         backvec_z2),
                                         pba->error_message,
                                         ptr->error_message);

            /* infer redhsift */
            double z2 = pba->a_today/backvec_z2[pba->index_bg_a]-1.;
            integ = 0.0;

            for(int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++){

              /* Call primordial power spectrum (Fourier space)*/
              class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, pgb2->k_bessel[index_k_bessel], &Pk), ppm->error_message, pgb2->error_message);
              //herehere
              //printf("tau_one = %g\n",tau_one );
              /*printf("%dx%dx%dx%dx%dx%d\n",
              index_type_first,
              index_type_second,
              index_l,
              index_alpha,
              index_r,
              index_k_bessel);
              printf("index_tau_first = %d, index_tau_second = %d \n", index_tau_first, index_tau_second);*/

              type1 = pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_first][index_k_bessel];
              //printf("type1 = %g\n", type1);
              //printf("type2 = %g\n", type2);
              type2 = pgb2->first_order_sources_integ[index_type_second][index_l][index_tau_second][index_k_bessel];

              /* Interpolate between the two indices found on the tau_sampling_cls grid to find the exact tau_one and tau_two */

              type1_interp = pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_first-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_tau_first]-tau_one)
                            + pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_first][index_k_bessel]*(tau_one-pgb2->tau_sampling_cls[index_tau_first-1]);
              type1_interp /= (pgb2->tau_sampling_cls[index_tau_first] - pgb2->tau_sampling_cls[index_tau_first-1]);

              type2_interp = pgb2->first_order_sources_integ[index_type_second][index_l][index_tau_second-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_tau_second]-tau_two)
                            + pgb2->first_order_sources_integ[index_type_second][index_l][index_tau_second][index_k_bessel]*(tau_two-pgb2->tau_sampling_cls[index_tau_second-1]);

              type2_interp /= (pgb2->tau_sampling_cls[index_tau_second] - pgb2->tau_sampling_cls[index_tau_second-1]);

              /*printf("minus(%g) = %g, interp(%g) = %g, plus(%g) = %g\n",
              pgb2->tau_sampling_cls[index_tau_first-1],
              pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_first-1][index_k_bessel],
              tau_one,
              type1_interp,
              pgb2->tau_sampling_cls[index_tau_first],
              pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_first][index_k_bessel]);*/

              /* Using area of trapezoid, we can sum a number of areas of trapezoids to approximate the integral */
              integ += pow(pgb2->k_bessel[index_k_bessel],-1.0) * 4. * _PI_ *  Pk * type1_interp * type2_interp  * pgb2->w_trapz_k[index_k_bessel];
            }
            pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r] = integ;
            //printf("integ = %g\n", integ );
          }

            /*class_call(integral(pba,
                                ppm,
                                pbs,
                                pgb2,
                                ppt,
                                ptr,
                                index_type_first,
                                index_type_second,
                                index_tau_first,
                                index_tau_second,
                                index_l,
                                &pvecback1,
                                &pvecback2,
                                &integ),
                       pgb2->error_message,
                       pgb2->error_message);*/



            //pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] = integ;

            //printf("%d    %g\n",ptr->l[index_l], ptr->l[index_l]*(ptr->l[index_l]+1)*pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second]/2/_PI_);

          }
        }
      }
    }
    /*for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
    //  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
        for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
          printf("%d    %g\n", ptr->l[index_l],  ptr->l[index_l]*( ptr->l[index_l]+1)*pgb2->Dl2[0][0][index_l][33][index_r]/2./_PI_ );
        }
    //  }
    }
    exit(0);*/

    /* Printing the integrand of the double time-integration, fixing r and varying alpha */
    //herehere
    int index_l_test1 = 0;
    int index_l_test2 = 4;
    int index_l_test3 = 9;
    int index_r_test1 = 0;
    int index_r_test2 = ceil(pgb2->r_size/2);
    int index_r_test3 = pgb2->r_size-1;

    int index_tau_first_polar1;
    int index_tau_first_polar2;
    int index_tau_first_polar3;

    int index_tau_second_polar1;
    int index_tau_second_polar2;
    int index_tau_second_polar3;

    double window_interp11, window_interp12;
    double window_interp21, window_interp22;
    double window_interp31, window_interp32;



    printf("# alpha    l%dr%g   l%dr%g    l%dr%g   l%dr%g    l%dr%g   l%dr%g    l%dr%g   l%dr%g    l%dr%g\n",
    ptr->l[index_l_test1],
    pgb2->r2[index_r_test1],
    ptr->l[index_l_test1],
    pgb2->r2[index_r_test2],
    ptr->l[index_l_test1],
    pgb2->r2[index_r_test3],
    ptr->l[index_l_test2],
    pgb2->r2[index_r_test1],
    ptr->l[index_l_test2],
    pgb2->r2[index_r_test2],
    ptr->l[index_l_test2],
    pgb2->r2[index_r_test3],
    ptr->l[index_l_test3],
    pgb2->r2[index_r_test1],
    ptr->l[index_l_test3],
    pgb2->r2[index_r_test2],
    ptr->l[index_l_test3],
    pgb2->r2[index_r_test3]);

    for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
      for(int index_r = 1; index_r < pgb2->r_size; index_r++){
        double tau_one_polar1 = tau0*(1.0-pgb2->r2[index_alpha][index_r_test1]*sin(pgb2->alpha[index_alpha]));
        double tau_two_polar1 = tau0*(1.0-pgb2->r2[index_alpha][index_r_test1]*cos(pgb2->alpha[index_alpha]));

        double tau_one_polar2 = tau0*(1.0-pgb2->r2[index_alpha][index_r_test2]*sin(pgb2->alpha[index_alpha]));
        double tau_two_polar2 = tau0*(1.0-pgb2->r2[index_alpha][index_r_test2]*cos(pgb2->alpha[index_alpha]));

        double tau_one_polar3 = tau0*(1.0-pgb2->r2[index_alpha][index_r_test3]*sin(pgb2->alpha[index_alpha]));
        double tau_two_polar3 = tau0*(1.0-pgb2->r2[index_alpha][index_r_test3]*cos(pgb2->alpha[index_alpha]));

        index_of_tau_sampling_selection(tau_one_polar1, 0, &index_tau_first_polar1, pgb2);
        index_of_tau_sampling_selection(tau_two_polar1, 0, &index_tau_second_polar1, pgb2);

        index_of_tau_sampling_selection(tau_one_polar2, 0, &index_tau_first_polar2, pgb2);
        index_of_tau_sampling_selection(tau_two_polar2, 0, &index_tau_second_polar2, pgb2);

        index_of_tau_sampling_selection(tau_one_polar3, 0, &index_tau_first_polar3, pgb2);
        index_of_tau_sampling_selection(tau_two_polar3, 0, &index_tau_second_polar3, pgb2);



        /* Interpolate between the two indices found on the tau_sampling_cls grid to find the exact tau_one and tau_two */

        window_interp11 = selection[0][index_tau_first_polar1-1]*(pgb2->tau_sampling_selection[0][index_tau_first_polar1]-tau_one_polar1)
                      + selection[0][index_tau_first_polar1]*(tau_one_polar1-pgb2->tau_sampling_selection[0][index_tau_first_polar1-1]);
        window_interp11 /= (pgb2->tau_sampling_selection[0][index_tau_first_polar1] - pgb2->tau_sampling_selection[0][index_tau_first_polar1-1]);

        window_interp12 = selection[0][index_tau_second_polar1-1]*(pgb2->tau_sampling_selection[0][index_tau_second_polar1]-tau_two_polar1)
                      + selection[0][index_tau_second_polar1]*(tau_two_polar1-pgb2->tau_sampling_selection[0][index_tau_second_polar1-1]);
        window_interp12 /= (pgb2->tau_sampling_selection[0][index_tau_second_polar1] - pgb2->tau_sampling_selection[0][index_tau_second_polar1-1]);

        window_interp21 = selection[0][index_tau_first_polar2-1]*(pgb2->tau_sampling_selection[0][index_tau_first_polar2]-tau_one_polar2)
                      + selection[0][index_tau_first_polar2]*(tau_one_polar2-pgb2->tau_sampling_selection[0][index_tau_first_polar2-1]);
        window_interp21 /= (pgb2->tau_sampling_selection[0][index_tau_first_polar2] - pgb2->tau_sampling_selection[0][index_tau_first_polar2-1]);

        window_interp22 = selection[0][index_tau_second_polar2-1]*(pgb2->tau_sampling_selection[0][index_tau_second_polar2]-tau_two_polar2)
                      + selection[0][index_tau_second_polar2]*(tau_two_polar2-pgb2->tau_sampling_selection[0][index_tau_second_polar2-1]);
        window_interp22 /= (pgb2->tau_sampling_selection[0][index_tau_second_polar2] - pgb2->tau_sampling_selection[0][index_tau_second_polar2-1]);

        window_interp31 = selection[0][index_tau_first_polar3-1]*(pgb2->tau_sampling_selection[0][index_tau_first_polar3]-tau_one_polar3)
                      + selection[0][index_tau_first_polar3]*(tau_one_polar3-pgb2->tau_sampling_selection[0][index_tau_first_polar3-1]);
        window_interp31 /= (pgb2->tau_sampling_selection[0][index_tau_first_polar3] - pgb2->tau_sampling_selection[0][index_tau_first_polar3-1]);

        window_interp32 = selection[0][index_tau_second_polar3-1]*(pgb2->tau_sampling_selection[0][index_tau_second_polar3]-tau_two_polar3)
                      + selection[0][index_tau_second_polar3]*(tau_two_polar3-pgb2->tau_sampling_selection[0][index_tau_second_polar3-1]);
        window_interp32 /= (pgb2->tau_sampling_selection[0][index_tau_second_polar3] - pgb2->tau_sampling_selection[0][index_tau_second_polar3-1]);


      }
    }



    // Interpolate to read exact Dl(z1,z2) given
    double z_test = 1.0;
    double tau1_test;
    double tau2_test;
    int index_alpha;
    int index_r;

    class_call(background_tau_of_z(
                            pba,
                            z_test,
                            &tau1_test),
                            ppt->error_message,
                            pgb2->error_message);

    class_call(background_tau_of_z(
                            pba,
                            z_test,
                            &tau2_test),
                            ppt->error_message,
                            pgb2->error_message);



    printf("tau1_test = %g\n", tau1_test);
    printf("tau2_test = %g\n", tau2_test);


    double alpha = atan((1.-(tau1_test/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))/(1.-(tau2_test/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1])));

    double r = sqrt((1.-(tau1_test/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))*(1.-(tau1_test/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))+(1.-(tau2_test/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))*(1.-(tau2_test/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1])));

    index_of_alpha(alpha, &index_alpha, pgb2);
    index_of_r(r, index_alpha, &index_r, pgb2);
    printf("alpha* = %g\n",alpha );
    printf("alpha_found_plus = %g\n",pgb2->alpha[index_alpha]);
    printf("alpha_found_minus = %g\n",pgb2->alpha[index_alpha-1]);
    printf("r* = %g\n",r );
    printf("r_found_neutral = %g\n",pgb2->r2[index_alpha][index_r]);
    printf("r_found_minus = %g\n",pgb2->r2[index_alpha][index_r-1]);
    printf("r_found_plus = %g\n",pgb2->r2[index_alpha][index_r+1]);
    double tau1_found = pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]*(1.-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
    double tau2_found = pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]*(1.-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));
    printf("tau1_found = %g\n", tau1_found);
    printf("tau2_found = %g\n", tau2_found);
    printf("tau1_found minus = %g\n", pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]*(1.-pgb2->r2[index_alpha][index_r-1]*sin(pgb2->alpha[index_alpha])));
    printf("tau2_found minus = %g\n", pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]*(1.-pgb2->r2[index_alpha][index_r-1]*cos(pgb2->alpha[index_alpha])));
    for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
      double norm2 = ptr->l[index_l]*(ptr->l[index_l]+1)/(2.*_PI_);
      double temp_minus = pgb2->Dl2[0][0][index_l][index_alpha-1][index_r-1]*(pgb2->alpha[index_alpha]-alpha)+
                            pgb2->Dl2[0][0][index_l][index_alpha][index_r-1]*(alpha-pgb2->alpha[index_alpha-1]);

      temp_minus /= (pgb2->alpha[index_alpha]-pgb2->alpha[index_alpha-1]);

      double temp_plus = pgb2->Dl2[0][0][index_l][index_alpha-1][index_r]*(pgb2->alpha[index_alpha]-alpha)+
                            pgb2->Dl2[0][0][index_l][index_alpha][index_r]*(alpha-pgb2->alpha[index_alpha-1]);

      temp_plus /= (pgb2->alpha[index_alpha]-pgb2->alpha[index_alpha-1]);


      double result = temp_minus*(pgb2->r2[index_alpha][index_r]-r)+temp_plus*(r-pgb2->r2[index_alpha][index_r-1]);

      result /= (pgb2->r2[index_alpha][index_r]-pgb2->r2[index_alpha][index_r-1]);


      printf("%d    %g\n",ptr->l[index_l], norm2*result);
    }


    /*double inner_wind_check;
    double outer_wind_check;
    double result_second_int;
    double result_first_int;

    for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
      for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
        for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
          for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
            for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
              printf("index_l = %d\n",index_l);
              /*if (index_l > 0) {
                printf("%d    %g\n",ptr->l[index_l-1], ptr->l[index_l-1]*(ptr->l[index_l-1]+1)*pgb2->Cl[0][0][index_l-1][0][0]/(2*_PI_));
              }*/

              /*result_second_int = 0.;
              outer_wind_check = 0.;
              for(int index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++){
                double tau2 = pgb2->tau_sampling_selection[bin2][index_tau_second];
                double result_first_int = 0.0;
                inner_wind_check = 0.0;

                for(int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++){
                  double tau1 = pgb2->tau_sampling_selection[bin1][index_tau_first];

                  //double alpha = atan((1.-(tau1/pba->conformal_age))/(1.-(tau2/pba->conformal_age)));

                  //NOTE should it be tau_sampling_cls?
                  double argument = (1.-(tau1/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))/(1.-(tau2/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]));
                  double alpha = atan(argument);

                  double r = sqrt((1.-(tau1/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))*(1.-(tau1/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))+(1.-(tau2/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))*(1.-(tau2/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1])));

                  if(isnan(argument)) {
                    printf("alpha = %g\n",alpha);
                    printf("r = %g\n",r);
                    printf("index_r = %d\n", index_r);
                    printf("index_alpha = %d\n", index_alpha);
                    printf("argument = %g\n", argument);
                    printf("is this one? = %g\n", (tau2/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]));
                    alpha = _PI_/2.;
                    printf("alpha = %g\n",alpha);
                  }

                  index_of_alpha(alpha, &index_alpha, pgb2);
                  index_of_r(r, index_alpha, &index_r, pgb2);
                  //printf("index_r =%d\n",index_r);
                  //printf("index_alpha = %d\n", index_alpha);

                  //printf("r_given = %g, r_found = %g, error = %g\n",r,pgb2->r2[index_alpha][index_r],100.*(r-pgb2->r2[index_alpha][index_r])/r);
                  //printf("first_alpha = %g, alpha_found = %g, error = %g\n",alpha, pgb2->alpha[index_alpha], 100.*(alpha-pgb2->alpha[index_alpha])/alpha);

                  double temp_minus = pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha-1][index_r-1]*(pgb2->alpha[index_alpha]-alpha)+
                                        pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r-1]*(alpha-pgb2->alpha[index_alpha-1]);

                  temp_minus /= (pgb2->alpha[index_alpha]-pgb2->alpha[index_alpha-1]);

                  double temp_plus = pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha-1][index_r]*(pgb2->alpha[index_alpha]-alpha)+
                                        pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r]*(alpha-pgb2->alpha[index_alpha-1]);

                  temp_plus /= (pgb2->alpha[index_alpha]-pgb2->alpha[index_alpha-1]);


                  double result = temp_minus*(pgb2->r2[index_alpha][index_r]-r)+temp_plus*(r-pgb2->r2[index_alpha][index_r-1]);

                  result /= (pgb2->r2[index_alpha][index_r]-pgb2->r2[index_alpha][index_r-1]);

                  /*printf("Dl2[%d][%d] = %g, Dl2[%d][%d] = %g, interp = %g\n",
                          index_alpha,
                          index_r,
                          pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r],
                          index_alpha-1,
                          index_r-1,
                          pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha-1][index_r-1],
                          result);
                  /*printf("alpha = %g\n",alpha);
                  printf("result = %g\n", result);
                  printf("r = %g\n",r );*/
                  /*printf("r-denominator =%g\n",(pgb2->r2[index_alpha][index_r]-pgb2->r2[index_alpha][index_r-1]));
                  printf("alpha-denominator = %g\n", (pgb2->alpha[index_alpha]-pgb2->alpha[index_alpha-1]));
                  */

              /*    result_first_int += result * w_trapz[bin1][index_tau_first]*selection[bin1][index_tau_first];

                  inner_wind_check += result * simps_weights[bin1][index_tau_first]*selection[bin1][index_tau_first];

                }
                //CHECK THE BINs
                outer_wind_check += inner_wind_check*simps_weights[bin2][index_tau_second]*selection[bin2][index_tau_second];

                result_second_int += result_first_int * w_trapz[bin2][index_tau_second]*selection[bin2][index_tau_second];
              }
              pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2] = result_second_int;
              printf("%d      %g\n", ptr->l[index_l], ptr->l[index_l]*(ptr->l[index_l]+1)*outer_wind_check/(2*_PI_));
            }
          }
        }
      }
    }
    printf("#alpha size = %d\n", pgb2->alpha_size);
    printf("#r size = %d\n", pgb2->r_size );
    printf("#tau_size_cls = %d\n", pgb2->tau_size_cls);
    printf("#tau_size_selection = %d\n", pgb2->tau_size_selection);
    printf("#z1=z2=%g (%g) Gaussian\n", ppt->selection_mean[0],ppt->selection_width[0]);
    printf("#l    l(l+1)Cl/2pi\n");
    for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
      printf("%d    %g\n",ptr->l[index_l], ptr->l[index_l]*(ptr->l[index_l]+1)*pgb2->Cl[0][0][index_l][0][0]/(2*_PI_));
    }

    printf("Reached here\n");*/



  int index_test_alpha = 6;
  int index_test_r = 1;
  int index_test_l =0;
  //double result;
  /*printf("#alpha[%d] = %g\n", index_test_alpha, pgb2->alpha[index_test_alpha]);
  printf("#l[%d] = %d\n", index_test_l, ptr->l[index_test_l] );
  printf("#dens-dens\n" );
  printf("#r     Dl");
  for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
    //tau_one = pba->conformal_age*(1.0-pgb2->r[index_r]*sin(pgb2->alpha[index_test_alpha]));
    //tau_two = pba->conformal_age*(1.0-pgb2->r[index_r]*cos(pgb2->alpha[index_test_alpha]));

    //index_of_tau_sampling_cls(tau_one, &index_tau_one_cls, pgb2);
    //index_of_tau_sampling_cls(tau_two, &index_tau_two_cls, pgb2);

    //tau_one_cls = pgb2->tau_sampling_cls[index_tau_one_cls];
    //tau_two_cls = pgb2->tau_sampling_cls[index_tau_two_cls];

    double result = pgb2->Dl2[pgb2->index_type_density][pgb2->index_type_density][index_test_l][index_test_alpha][index_r];
    printf("%g      %g\n",pgb2->r[index_r], result);
  }*/

  printf("#dens-rsd\n" );
  printf("#alpha     Dl");

  // Warning: This only works for auto-term correlations index_type_first = index_type_second.
  /*for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    for(int index_type_second = 0; index_type_second < index_type_first+1; index_type_second++){
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
        for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_cls; index_tau_first++){
          for(int index_tau_second = index_tau_first+1; index_tau_second < pgb2->tau_size_cls; index_tau_second++){
            pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] = pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_second][index_tau_first];
          }
        }
      }
    }
  }*/

  /*for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_cls; index_tau_first++) {
    for (int index_tau_second = 0; index_tau_second < pgb2->tau_size_cls; index_tau_second++) {
      norm_test = ptr->l[test_index_l]*(ptr->l[test_index_l]+1.)/(2*_PI_);
      fprintf(file,
              "%g     %g      %g\n",
              pgb2->tau_sampling_cls[index_tau_first],
              pgb2->tau_sampling_cls[index_tau_second]);
    }
  }*/
  //fprintf(file, "K_MAX = %g \n",pgb2->k_bessel[pgb2->k_size_bessel-1] );
  //fclose(file);
  //printf("Written Dl data into file.\n");
  printf("Done k-integration.\n");



  /***************************************************************************
  ==========================   Time-Integration ==============================
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
  double tau1, tau2;

  /* If the window function is a Dirac-delta function then there is no need to integrate over time,
  intead we just interpolate the would-be integrand pgb2->Dl and write its value at the tau values
  set by the redshift bins into pgb2->Cl. */

  if (ppt->selection == dirac) {
    for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {

      class_call(background_tau_of_z(
                              pba,
                              ppt->selection_mean[bin1],
                              &tau1),
                              ppt->error_message,
                              pgb2->error_message);
      printf("tau1(z=1.0) = %g\n", tau1 );
      // ALERT CROSS BINS
      for (int bin2 = 1; bin2 < ppt->selection_num; bin2++) {

        class_call(background_tau_of_z(
                                pba,
                                ppt->selection_mean[bin2],
                                &tau2),
                                ppt->error_message,
                                pgb2->error_message);
        printf("tau2(z=1.0) = %g\n", tau2);

        for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
          for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
            for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
              pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2] = 0.;

              index_of_tau_sampling_cls(tau1, &index_tau1,  pgb2);
              index_of_tau_sampling_cls(tau2, &index_tau2,  pgb2);

              temp_minus = pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1-1][index_tau2-1]*(pgb2->tau_sampling_cls[index_tau1]-tau1)
                            + pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1][index_tau2-1]*(tau1-pgb2->tau_sampling_cls[index_tau1-1]);
              temp_minus /= (pgb2->tau_sampling_cls[index_tau1] - pgb2->tau_sampling_cls[index_tau1-1]);

              temp_plus = pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1-1][index_tau2]*(pgb2->tau_sampling_cls[index_tau1]-tau1)
                            + pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1][index_tau2]*(tau1-pgb2->tau_sampling_cls[index_tau1-1]);
              temp_plus /= (pgb2->tau_sampling_cls[index_tau1] - pgb2->tau_sampling_cls[index_tau1-1]);

              temp = temp_minus * (pgb2->tau_sampling_cls[index_tau2]-tau2)
                            + temp_plus * (tau2-pgb2->tau_sampling_cls[index_tau2-1]);

              temp /= (pgb2->tau_sampling_cls[index_tau2] - pgb2->tau_sampling_cls[index_tau2-1]);

              pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2] = temp;
            }
          }
        }
      }
    }
  }


  else{
    double inner_polar = 0.0;
    double outer_polar = 0.0;
    double window1_interp;
    double window2_interp;
    double window1;
    double window2;

    double tau_one_polar, tau_two_polar;
    int index_tau_first_polar, index_tau_second_polar;
    printf("Entered non-Dirac window function integration.\n");
    printf("Starting polar time-integration\n");
    for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++) {
      for (int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++) {
        for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
          for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
            for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
              double outer_integration = 0.0;
              double window_int_check2 = 0.0;
              printf("index_l = %d\n", index_l);
              for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
                printf("index_alpha = %d\n", index_alpha );
                double inner_integration = 0.0;
                double temp2 = 0.0;
                for (int index_r = 0; index_r < pgb2->r_size; index_r++) {

                  tau_one_polar = tau0*(1.0-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
                  tau_two_polar = tau0*(1.0-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));

                  index_of_tau_sampling_selection(tau_one_polar, 0, &index_tau_first_polar, pgb2);
                  index_of_tau_sampling_selection(tau_two_polar, 0, &index_tau_second_polar, pgb2);

                  window1 = selection[0][index_tau_first_polar];
                  window2 = selection[0][index_tau_second_polar];


                  /* Interpolate between the two indices found on the tau_sampling_cls grid to find the exact tau_one and tau_two */

                  window1_interp = selection[0][index_tau_first_polar-1]*(pgb2->tau_sampling_selection[0][index_tau_first_polar]-tau_one_polar)
                                + selection[0][index_tau_first_polar]*(tau_one_polar-pgb2->tau_sampling_selection[0][index_tau_first_polar-1]);
                  window1_interp /= (pgb2->tau_sampling_selection[0][index_tau_first_polar] - pgb2->tau_sampling_selection[0][index_tau_first_polar-1]);

                  window2_interp = selection[0][index_tau_second_polar-1]*(pgb2->tau_sampling_selection[0][index_tau_second_polar]-tau_two_polar)
                                + selection[0][index_tau_second_polar]*(tau_two_polar-pgb2->tau_sampling_selection[0][index_tau_second_polar-1]);
                  window2_interp /= (pgb2->tau_sampling_selection[0][index_tau_second_polar] - pgb2->tau_sampling_selection[0][index_tau_second_polar-1]);


                  inner_integration += tau0*tau0*pgb2->r2[index_alpha][index_r]*window1_interp*window2_interp*w_trapz_r2[index_alpha][index_r]
                                          *pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r];

                  Dl2_integrand[index_type_first][index_type_second][index_l][index_alpha][index_r] = tau0*tau0
                                        *pgb2->r2[index_alpha][index_r]*window1_interp*window2_interp
                                          *pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r];
                  temp2 += tau0*tau0*pgb2->r2[index_alpha][index_r]*window1_interp*window2_interp*w_trapz_r2[index_alpha][index_r];
                  /*type2_interp = pgb2->first_order_sources_integ[index_type_second][index_l][index_tau_second-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_tau_second]-tau_two)
                                + pgb2->first_order_sources_integ[index_type_second][index_l][index_tau_second][index_k_bessel]*(tau_two-pgb2->tau_sampling_cls[index_tau_second-1]);

                  type2_interp /= (pgb2->tau_sampling_cls[index_tau_second] - pgb2->tau_sampling_cls[index_tau_second-1]);*/
                  /*printf("window_start1(%g) = %g, window_exact1(%g) = %g, window_end1(%g) = %g\n",
                      pgb2->tau_sampling_selection[0][index_tau_first_polar-1],
                      selection[0][index_tau_first_polar-1],
                      tau_one_polar,
                      window1_interp,
                      pgb2->tau_sampling_selection[0][index_tau_first_polar],
                      selection[0][index_tau_first_polar]);

                  printf("window_start2(%g) = %g, window_exact2(%g) = %g, window_end2(%g) = %g\n",
                      pgb2->tau_sampling_selection[0][index_tau_second_polar-1],
                      selection[0][index_tau_second_polar-1],
                      tau_two_polar,
                      window2_interp,
                      pgb2->tau_sampling_selection[0][index_tau_second_polar],
                      selection[0][index_tau_second_polar]);


                  printf("inner_integration = %g, Dl = %g, w_trapz_r2 = %g \n",
                  inner_integration,
                  pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r],
                  w_trapz_r2[index_alpha][index_r]);*/

                  //printf("%g      %g      %g      %g\n",tau_one_polar, window1_interp, tau_two_polar, window2_interp);
                }
                outer_integration += inner_integration*w_trapz_alpha[index_alpha];
                window_int_check2 += temp2*w_trapz_alpha[index_alpha];
              }
              pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2] = outer_integration;
              //printf("*** l = %d, Cl = %g\n",ptr->l[index_l], pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2]);


              printf("window_int_check2 = %g\n", window_int_check2);
            }
          }
        }
      }
    }
    int test_index_l1 = 0;
    int test_index_l2 = 4;
    int test_index_l3 = 8;

    int test_index_r1 = 2;
    int test_index_r2 = 7;
    int test_index_r3 = 13;
    /*printf("#alpha  (l=%d,r=%g)  (l=%d,r=%g)  (l=%d,r=%g)  (l=%d,r=%g)  (l=%d,r=%g)  (l=%d,r=%g)  (l=%d,r=%g)  (l=%d,r=%g)  (l=%d,r=%g)\n",
            ptr->l[test_index_l1], (test_index_r1+1.)/pgb2->r_size, ptr->l[test_index_l1], (test_index_r2+1.)/pgb2->r_size, ptr->l[test_index_l1], (test_index_r3+1.)/pgb2->r_size,
            ptr->l[test_index_l2], (test_index_r1+1.)/pgb2->r_size, ptr->l[test_index_l2], (test_index_r2+1.)/pgb2->r_size, ptr->l[test_index_l2], (test_index_r3+1.)/pgb2->r_size,
            ptr->l[test_index_l3], (test_index_r1+1.)/pgb2->r_size, ptr->l[test_index_l3], (test_index_r2+1.)/pgb2->r_size, ptr->l[test_index_l3], (test_index_r3+1.)/pgb2->r_size);

    for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
      /*printf("%g    %g    %g    %g    %g    %g    %g    %g    %g    %g\n",
            pgb2->alpha[index_alpha],
            Dl2_integrand[0][0][test_index_l1][index_alpha][test_index_r1],
            Dl2_integrand[0][0][test_index_l1][index_alpha][test_index_r2],
            Dl2_integrand[0][0][test_index_l1][index_alpha][test_index_r3],
            Dl2_integrand[0][0][test_index_l2][index_alpha][test_index_r1],
            Dl2_integrand[0][0][test_index_l2][index_alpha][test_index_r2],
            Dl2_integrand[0][0][test_index_l2][index_alpha][test_index_r3],
            Dl2_integrand[0][0][test_index_l3][index_alpha][test_index_r1],
            Dl2_integrand[0][0][test_index_l3][index_alpha][test_index_r2],
            Dl2_integrand[0][0][test_index_l3][index_alpha][test_index_r3]);*/
  //  }

    FILE *file2 = fopen("Polar_time_integration_data.dat", "w");
    if (file2 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }


    int test_index_l = 5;
    fprintf(file2, "#First two columns are conformal time, the other three columns are tau0*tau0*r*Dl*Window1*Window2\n");

    fprintf(file2, "#tau1     tau2      D%d_integrand      D%d_integrand      D%d_integrand\n",ptr->l[test_index_l1], ptr->l[test_index_l2], ptr->l[test_index_l3] );
    for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
      for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
        double tau_pol1 = tau0*(1.0-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
        double tau_pol2 = tau0*(1.0-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));
        fprintf(file2,
                "%g      %g      %g      %g      %g\n",
                tau_pol1,
                tau_pol2,
                Dl2_integrand[0][0][test_index_l1][index_alpha][index_r],
                Dl2_integrand[0][0][test_index_l2][index_alpha][index_r],
                Dl2_integrand[0][0][test_index_l3][index_alpha][index_r]);
      }
    }
    fclose(file2);

    printf("#alpha size = %d\n", pgb2->alpha_size);
    printf("#log spacing base = %g\n",base );
    printf("#r size = %d\n", pgb2->r_size );
    printf("#tau_size_cls = %d\n", pgb2->tau_size_cls);
    printf("#tau_size_selection = %d\n", pgb2->tau_size_selection);
    printf("#z1=z2=%g (%g) Gaussian\n", ppt->selection_mean[0],ppt->selection_width[0]);
    printf("#l    l(l+1)Cl/2pi\n");
    for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {


      printf("%d    %g\n",ptr->l[index_l], ptr->l[index_l]*(ptr->l[index_l]+1)*pgb2->Cl[0][0][index_l][0][0]/(2*_PI_));
    }
    exit(0);





    double dummy_inner;
    double dummy_outer;

    for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
      for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
        for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
          for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
            for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){

              pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2] = 0.;
              //dummy_outer = 0.0;
              double temp456 = 0.0;
      /* Final result Cl[index_type_first][index_type_second][index_l][bin1][bin2] is defined on the pgb2->tau_sampling_cls, so an interpolation is required
            from Cl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] which is defined on pgb2->tau_sampling_cls_selection[index_tau]*/
              for(index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++){

                double temp123 = 0.;
                dummy_inner = 0.0;
                double tau2 = pgb2->tau_sampling_selection[bin2][index_tau_second];
                //index_of_tau_sampling_cls(tau2, &index_tau2, pgb2);


                for(index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++){
                  //printf("indexing %dx%dx%dx%dx%dx%dx%d\n",index_type_first, index_type_second, bin1, bin2, index_l,index_tau_second,index_tau_first);

                  double tau1 = pgb2->tau_sampling_selection[bin1][index_tau_first];

                  /*index_of_tau_sampling_cls(tau1, &index_tau1, pgb2);

                  temp_minus = pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1-1][index_tau2-1]*(pgb2->tau_sampling_cls[index_tau1]-tau1)
                                + pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1][index_tau2-1]*(tau1-pgb2->tau_sampling_cls[index_tau1-1]);
                  temp_minus /= (pgb2->tau_sampling_cls[index_tau1] - pgb2->tau_sampling_cls[index_tau1-1]);

                  temp_plus = pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1-1][index_tau2]*(pgb2->tau_sampling_cls[index_tau1]-tau1)
                                + pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1][index_tau2]*(tau1-pgb2->tau_sampling_cls[index_tau1-1]);
                  temp_plus /= (pgb2->tau_sampling_cls[index_tau1] - pgb2->tau_sampling_cls[index_tau1-1]);

                  temp = temp_minus * (pgb2->tau_sampling_cls[index_tau2]-tau2)
                                + temp_plus * (tau2-pgb2->tau_sampling_cls[index_tau2-1]);

                  temp /= (pgb2->tau_sampling_cls[index_tau2] - pgb2->tau_sampling_cls[index_tau2-1]);

                  temp123 += temp * w_trapz[bin1][index_tau_first]*selection[bin1][index_tau_first];

                  dummy_inner += w_trapz[bin1][index_tau_first]*selection[bin1][index_tau_first];*/

                  /* r & alpha interpolation to find specified tau1 and tau2 defined by pgb2->tau_sampling_selection */


                  double alpha2 = atan((1.-(tau1/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))/(1.-(tau2/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1])));

                  double r2 = sqrt((1.-(tau1/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))*(1.-(tau1/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))+(1.-(tau2/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))*(1.-(tau2/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1])));

                  index_of_alpha(alpha2, &index_alpha, pgb2);
                  index_of_r(r2, index_alpha, &index_r, pgb2);


                  //double tau1_found = pba->conformal_age*(1.-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
                  //double tau2_found = pba->conformal_age*(1.-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));


                  double temp_minus2 = pgb2->Dl2[bin1][bin2][index_l][index_alpha-1][index_r-1]*(pgb2->alpha[index_alpha]-alpha)+
                                        pgb2->Dl2[bin1][bin2][index_l][index_alpha][index_r-1]*(alpha-pgb2->alpha[index_alpha-1]);

                  temp_minus2 /= (pgb2->alpha[index_alpha]-pgb2->alpha[index_alpha-1]);

                  double temp_plus2 = pgb2->Dl2[bin1][bin2][index_l][index_alpha-1][index_r]*(pgb2->alpha[index_alpha]-alpha)+
                                        pgb2->Dl2[bin1][bin2][index_l][index_alpha][index_r]*(alpha-pgb2->alpha[index_alpha-1]);

                  temp_plus2 /= (pgb2->alpha[index_alpha]-pgb2->alpha[index_alpha-1]);


                  double result2 = temp_minus*(pgb2->r2[index_alpha][index_r]-r)+temp_plus*(r-pgb2->r2[index_alpha][index_r-1]);

                  result2 /= (pgb2->r2[index_alpha][index_r]-pgb2->r2[index_alpha][index_r-1]);



                  temp123 += result2 * w_trapz[bin1][index_tau_first]*selection[bin1][index_tau_first];

                }

                dummy_outer += dummy_inner * w_trapz[bin2][index_tau_second] * selection[bin2][index_tau_second];

                temp456 += temp123 * w_trapz[bin2][index_tau_second] * selection[bin2][index_tau_second];
                pgb2->integral_over_single_window[index_type_first][index_type_second][index_l][bin1][bin2][index_tau_second] = temp123;
              }
              pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2] = temp456;
            }
          }
        }
      }
    }
  }
//New exit!!


  /* alpha and r time parameterisation integration */
  /*else{
    printf("Entered non-Dirac window function integration.\n");
    printf("Computing Cls via r and alpha integration...\n");

    double test_integ;
    double temp456;
    double temp123;
    printf("#l    r   int tau0*tau0*r*Dl da dr\n");
    for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
      for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
        //for (int bin1 = 0; bin1 < 1/*ppt->selection_num*///; bin1++) {
          //for (int bin2 = 1; bin2 < 2/*ppt->selection_num*/; bin2++) {
            //for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){

              //pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2] = 0.;

              //temp456 = 0.0;

      /* Final result Cl[index_type_first][index_type_second][index_l][bin1][bin2] is defined on the pgb2->tau_sampling_cls, so an interpolation is required
            from Cl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] which is defined on pgb2->tau_sampling_cls[index_tau]*/
            /*  for (int index_r_up = 1; index_r_up < pgb2->r_size; index_r_up++) {
                double result = 0.0;


                for(int index_r = 0; index_r < index_r_up+1; index_r++){
                  temp123 = 0.0;
                  test_integ = 0.0;

                  for(int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++){

                    double tau1 = pba->conformal_age*(1.-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
                    double tau2 = pba->conformal_age*(1.-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));

                    index_of_tau_sampling_selection(tau1, bin1, &index_tau_first, pgb2);
                    index_of_tau_sampling_selection(tau2, bin2, &index_tau_second, pgb2);

                    test_integ +=pba->conformal_age*pba->conformal_age*pgb2->r2[index_alpha][index_r]
                                    *pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r]
                                      *w_trapz_alpha[index_alpha];



                    /* Interpolation */

                    /*temp_minus = pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1-1][index_tau2-1]*(pgb2->tau_sampling_cls[index_tau1]-tau1)
                                  + pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1][index_tau2-1]*(tau1-pgb2->tau_sampling_cls[index_tau1-1]);
                    temp_minus /= (pgb2->tau_sampling_cls[index_tau1] - pgb2->tau_sampling_cls[index_tau1-1]);

                    temp_plus = pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1-1][index_tau2]*(pgb2->tau_sampling_cls[index_tau1]-tau1)
                                  + pgb2->Dl[index_type_first][index_type_second][index_l][index_tau1][index_tau2]*(tau1-pgb2->tau_sampling_cls[index_tau1-1]);
                    temp_plus /= (pgb2->tau_sampling_cls[index_tau1] - pgb2->tau_sampling_cls[index_tau1-1]);

                    temp = temp_minus * (pgb2->tau_sampling_cls[index_tau2]-tau2)
                                  + temp_plus * (tau2-pgb2->tau_sampling_cls[index_tau2-1]);

                    temp /= (pgb2->tau_sampling_cls[index_tau2] - pgb2->tau_sampling_cls[index_tau2-1])

                    //printf("pgb2->Dl2 = %g\n", pgb2->Dl2[pgb2->index_type_density][pgb2->index_type_density][index_l][index_alpha][index_r]);
                    temp123 += pba->conformal_age*pba->conformal_age*pgb2->r[index_r]* pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r]
                                * w_trapz_r[index_r] * selection[bin1][index_tau_first]*selection[bin2][index_tau_second];



                  }

                  //temp456 += temp123 * w_trapz_alpha[index_alpha]*w_trapz_r[index_r];
                  result += test_integ*w_trapz_r[index_r];

                }
                printf("%d  %d  %g\n",ptr->l[index_l], index_r_up, result);
              }

              //pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2] = temp456;
            }
          //}
        //}
      }
    }
  }

*/
printf("integrating between k =%g and %g\n",pgb2->k_bessel[index_k_bessel],pgb2->k_bessel[pgb2->k_size_bessel-1]);
printf("#l       l(l+1)C_l/2pi\n");

  double norm;
  for  (index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
    norm = ptr->l[index_l]*(ptr->l[index_l]+1.)/(2*_PI_);

    /*  printf("%d    %g    %g    %g    %g    %g    %g\n",
              ptr->l[index_l],
              norm*pgb2->Cl[pgb2->index_type_density][pgb2->index_type_density][index_l][0][0],
              norm*pgb2->Cl[pgb2->index_type_rsd][pgb2->index_type_rsd][index_l][0][0],
              norm*pgb2->Cl[pgb2->index_type_d1][pgb2->index_type_d1][index_l][0][0],
              norm*pgb2->Cl[pgb2->index_type_d2][pgb2->index_type_d2][index_l][0][0],
              norm*pgb2->Cl[pgb2->index_type_lens][pgb2->index_type_lens][index_l][0][0],
              norm*(pgb2->Cl[pgb2->index_type_g1][pgb2->index_type_g1][index_l][0][0],
              +pgb2->Cl[pgb2->index_type_g2][pgb2->index_type_g2][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g3][pgb2->index_type_g3][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g4][pgb2->index_type_g4][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g5][pgb2->index_type_g5][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g5][pgb2->index_type_g1][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g4][pgb2->index_type_g1][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g3][pgb2->index_type_g1][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g2][pgb2->index_type_g1][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g5][pgb2->index_type_g2][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g4][pgb2->index_type_g2][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g3][pgb2->index_type_g2][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g5][pgb2->index_type_g3][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g4][pgb2->index_type_g3][index_l][0][0]
              +pgb2->Cl[pgb2->index_type_g5][pgb2->index_type_g4][index_l][0][0]));*/

    printf("%d    %g\n",
           ptr->l[index_l],
           norm*pgb2->Cl[pgb2->index_type_rsd][pgb2->index_type_rsd][index_l][0][0]);
    }



  printf("Starting bispectrum computation.\n");


/***************************************************************************
=====================   Bispectrum Computation =============================
****************************************************************************/

  int index_l1, index_l2, index_l3;

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
