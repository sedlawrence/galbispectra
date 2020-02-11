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

                 double output = (tau - pgb2->tau_sampling_cls[0])/(pgb2->tau_sampling_cls[pgb2->tau_size_selection-1] - pgb2->tau_sampling_cls[0]) * (pgb2->tau_size_selection-1);
                 //printf("output = %g\n",output );
                 * index_tau = (int) ceil(output);

                 if (* index_tau < 1) {
                   * index_tau = 1;
                 }

                 else if(* index_tau > pgb2->tau_size_selection -1){
                   * index_tau = pgb2->tau_size_selection -1;

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
                * index_tau = (int) ceil (output);

                if (* index_tau < 1) {
                  * index_tau = 1;
                }

                else if(* index_tau > pgb2->tau_size_selection -1){
                  * index_tau = pgb2->tau_size_selection -1;

                }
              }

int index_of_tau_rp(double tau1,
                int * index_r_optimal,
                int * index_alpha_optimal,
                struct galbispectra2 * pgb2){

                int index_r = -1;
                int index_alpha = -1;
                double  min_diff;
                min_diff = 15000.0;
                double diff;

                for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
                  for (int index_alpha = pgb2->alpha_size; index_alpha < 2*pgb2->alpha_size-1; index_alpha++) {
                    diff = abs(tau1-pgb2->tau_rp[index_r][index_alpha]);
                    printf("tau_rp = %g\n", pgb2->tau_rp[index_r][index_alpha]);
                    printf("diff = %g\n",diff);
                    if (diff < min_diff) {
                      printf("new diff = %g\n", min_diff);
                      min_diff = diff;
                      *index_r_optimal = index_r;
                      *index_alpha_optimal = index_alpha;
                    }
                    else{
                      continue;
                    }
                  }
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
  /* Please ensure the pgb2->tau_size_selection grid is of a higher resolution than pgb2->r_size * pgb2->alpha_size = 13*/
  pgb2->tau_size_selection = 1000;

  pgb2->r_size = 75; //previously on 75

  /* alpha_size must be an odd positive-integer in order to correctly fill the pgb2->r and pgb2->alpha grids using symmetries */
  pgb2->alpha_size = 13; //previously on 12

  if (pgb2->alpha_size % 2 == 0) {
    printf("ERROR! alpha_size must be an ODD positive-integer\n" );
    exit(0);
  }

  /* We wish to double the number of slots such that points in the two grids tau_sampling_cls and tau_sampling_bessel align. */
  pgb2->tau_size_bessel = bessel_boost * (pgb2->tau_size_selection-1) + 1;
  pgb2->k_size_bessel = 25000;  /*tau_size_bessel*/   /*k_size * boost*/
  printf("Starting galaxy bispectra module...\n");
  if (ppt->selection==gaussian) {
    printf("We have a Gaussian Window Function on our hands.\n");
  }
  if (ppt->selection==dirac) {
    printf("We have a Dirac Window Function on our hands.\n");
  }

  double * result;

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
              pgb2->tau_size_selection * sizeof(double),
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
    // ALERT SHOULD THE TIME GRIDS RUN TO conformal_age or tau_max
    for (index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
      pgb2->tau_sampling_selection[bin][index_tau] = tau_min + index_tau*(tau_max-tau_min)/(pgb2->tau_size_selection-1);
      //printf("tau[%d][%d] = %g\n", bin, index_tau, pgb2->tau_sampling_selection[bin][index_tau]);
    }
  }



  for (index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
    pgb2->tau_sampling_cls[index_tau] = overall_tau_min + index_tau*(pba->conformal_age-overall_tau_min)/(pgb2->tau_size_selection-1);

  }
  /* Find the ranges of alpha to span such that we sample within tau_sampling_cls[0] to tau_sampling_cls[pgb2->tau_size_selection-1]*/

  double r_max;


  /* Under new parameterisation tau1 = tau0(1-rsin(alpha)). tau2 = tau0(1-rcos(alpha)). Max r value can be found by
  rearranging tau_min = tau1(r_max,pi/2)=tau2(r_max,0). This is assuming alpha runs over 0 to pi/2. */

  r_max = (1.0-pgb2->tau_sampling_cls[0]/pba->conformal_age);


  printf("pgb2->tau_sampling_cls[0]=%g\n",pgb2->tau_sampling_cls[0]);
  printf("conformal_age = %g\n",pba->conformal_age);

  // Time-grid values for tau1 = tau0(1-r*sin(alpha)), tau2 = tau0(1-r*cos(alpha)). The motivation for this parameterisation is to
  // effectively sample highly oscillatory regions of upcoming integrands. We sample one dimension with fine precision and the other
  // with coarse precision. This way when we integrate over alpha and r we have fewer iterations since alpha_size * r_size < tau_size_selection^2.
  printf("alpha size = %d\n", pgb2->alpha_size);
  printf("r size = %d\n", pgb2->r_size );



  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {

    pgb2->alpha[index_alpha] =  index_alpha*(_PI_/2.0)/(pgb2->alpha_size-1);
    //printf("pgb2->alpha[%d] = %g,\n", index_alpha, pgb2->alpha[index_alpha]);
  }

  for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
    pgb2->r[index_r]= index_r*(r_max/(pgb2->r_size-1));
    //printf("pgb2->r[%d] = %g\n", index_r, pgb2->r[index_r]);
  }
  /* Absolute maximum value r can have is when alpha =pi/4 and tau is set to tau_min. If r1 or r2 are above this value, this is because of a sing-
  ularity in the denominator */
  double tau_one;
  double tau_two;
  double tau_one_cls;
  double tau_two_cls;
  int index_tau_one_cls;
  int index_tau_two_cls;

  double  r_abs_max;

  r_abs_max = (1.0 - pgb2->tau_sampling_cls[0]/pba->conformal_age)/sin(_PI_/4.);
  printf("r_abs_max = %g\n", r_abs_max);
  printf("r_abs_max = %g\n", (pba->conformal_age-pgb2->tau_sampling_cls[0])/pba->conformal_age*cos(_PI_/4));

  double *r_max_per_alpha;
  int middle_alpha_index = (pgb2->alpha_size+1)/2-1;


  class_alloc(r_max_per_alpha,
              pgb2->alpha_size * sizeof(double*),
              pgb2->error_message);

  for (int index_alpha = 0; index_alpha < (pgb2->alpha_size+1)/2; index_alpha++) {
    r_max_per_alpha[index_alpha] = (pba->conformal_age-pgb2->tau_sampling_cls[0])/(pba->conformal_age*cos(pgb2->alpha[index_alpha]));
  }

  r_max_per_alpha[middle_alpha_index] = r_abs_max;

  for (int index_alpha = middle_alpha_index+1; index_alpha < pgb2->alpha_size; index_alpha++) {
    int distance_from_middle = index_alpha - middle_alpha_index;
    r_max_per_alpha[index_alpha] = r_max_per_alpha[middle_alpha_index-distance_from_middle];
    //printf("distance_from_middle = %d, set equal to index = %d\n", distance_from_middle, middle_alpha_index-distance_from_middle );
    //printf("r_max_per_alpha[%d] = %g, alpha = %g\n", index_alpha, r_max_per_alpha[index_alpha], pgb2->alpha[index_alpha]);
  }

  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    //printf("r_max_per_alpha[%d] = %g, alpha = %g\n", index_alpha, r_max_per_alpha[index_alpha], pgb2->alpha[index_alpha]);
  }





  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
      double r1 = (1.0 -  pgb2->tau_sampling_cls[0]/pba->conformal_age)/sin(pgb2->alpha[index_alpha]);
      double r2 = (1.0 - pgb2->tau_sampling_cls[0]/pba->conformal_age)/cos(pgb2->alpha[index_alpha]);

      //printf("r1 = %g, r2 = %g\n", r1, r2 );
      //printf("alpha = %g, *r_max = %g*\n", pgb2->alpha[index_alpha], r_max );
      pgb2->r2[index_alpha][index_r]= index_r*(r_max_per_alpha[index_alpha]/(pgb2->r_size-1));
      //printf("pgb2->r[%d] = %g\n", index_r, pgb2->r[index_r]);
      //printf("%g    %g\n",);
    }
  }
  free(r_max_per_alpha);
  //exit(0);




  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++){
    //printf("index_r = %d\n",index_r );

    for (int index_r = 0; index_r < pgb2->r_size; index_r++) {

      tau_one = pba->conformal_age*(1.0-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
      tau_two = pba->conformal_age*(1.0-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));

      index_of_tau_sampling_cls(tau_one, &index_tau_one_cls,  pgb2);
      index_of_tau_sampling_cls(tau_two, &index_tau_two_cls, pgb2);

      tau_one_cls = pgb2->tau_sampling_cls[index_tau_one_cls];
      tau_two_cls = pgb2->tau_sampling_cls[index_tau_two_cls];
      //printf("(tau_one,tau_cls,error)=(%g,%g,%g)\n",tau_one,tau_one_cls,(tau_one-tau_one_cls)*100./tau_one);
      //printf("(tau_two,tau_cls,error)=(%g,%g,%g)\n",tau_two,tau_two_cls,(tau_two-tau_two_cls)*100./tau_two);
      //printf("%g    %g\n", tau_one, tau_two);
    }
  }


  printf("tau_max = %g\n", pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1] );

  int test_tau_index1;
  int test_tau_index2;
  double test_time1;
  double test_time2;
  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
      test_time1 = pba->conformal_age*(1.0-pgb2->r[index_r]*sin(pgb2->alpha[index_alpha]));
      test_time2 = pba->conformal_age*(1.0-pgb2->r[index_r]*cos(pgb2->alpha[index_alpha]));

      index_of_tau_sampling_selection(test_time1, 0, &test_tau_index1, pgb2);
      index_of_tau_sampling_selection(test_time2, 0, &test_tau_index2, pgb2);

      /*printf("test_time1, tau_selection1, error1 = %g, %g, %g\n",
                test_time1,
                pgb2->tau_sampling_selection[0][test_tau_index1],
                100.0*(test_time1-pgb2->tau_sampling_selection[0][test_tau_index1])/test_time1);

      printf("test_time2, tau_selection2, error2 = %g, %g, %g\n",
                test_time2,
                pgb2->tau_sampling_selection[0][test_tau_index2],
                100.0*(test_time2-pgb2->tau_sampling_selection[0][test_tau_index2])/test_time2);*/
    }
  }



  /*class_alloc(pgb2->tau_rp,
              pgb2->r_size * sizeof(double*),
              ppt->error_message);

  for (int index_r = 0; index_r < pgb2->r_size; index_r++) {

    class_alloc(pgb2->tau_rp[index_r],
              pgb2->alpha_size * sizeof(double),
              ppt->error_message);

    for (int index_alpha = 0; index_alpha < 2*pgb2->alpha_size-1; index_alpha++) {
      pgb2->tau_rp[index_r][index_alpha] = pba->conformal_age*(1.0-pgb2->r[index_r]*sin(pgb2->alpha[index_alpha]));
      //printf("%g, ",tau_rp[index_r][index_alpha]);
    }
  }
  for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
    for (int index_alpha = pgb2->alpha_size; index_alpha < 2*pgb2->alpha_size-1; index_alpha++) {
      pgb2->tau_rp[index_r][index_alpha] = pba->conformal_age*(1.0-pgb2->r[index_r]*sin(pgb2->alpha[index_alpha]));
      //printf("%g, %d, %d\n ",pgb2->tau_rp[index_r][index_alpha],index_r,index_alpha);
    }
  }*/




  for (int index_tau_bessel = 0; index_tau_bessel < pgb2->tau_size_bessel; index_tau_bessel++) {
    pgb2->tau_sampling_bessel[index_tau_bessel] = overall_tau_min + index_tau_bessel*(pba->conformal_age-overall_tau_min)/(pgb2->tau_size_bessel-1);
  }
  //exit(0);
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

  double * w_trapz_r;
  double ** w_trapz_r2;
  double * w_trapz_alpha;
  /* Allocate and fill array for the trapezoidal weights for Chi integration w_trapz[bin][index_tau] */
  class_alloc(w_trapz,
              ppt->selection_num * sizeof(double*),
              ppt->error_message);

  class_alloc(w_trapz_r,
              pgb2->r_size * sizeof(double),
              ppt->error_message);

  class_alloc(w_trapz_alpha,
              pgb2->alpha_size * sizeof(double),
              ppt->error_message);

  for ( bin = 0; bin < ppt->selection_num; bin++) {
    class_alloc(w_trapz[bin],
                pgb2->tau_size_selection * sizeof(double),
                ppt->error_message);
  }

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



  /* Allocate and fill array for the trapezoidal weights for chi integration in the lensing term w_trapz_lens[index_tau] */
  class_alloc(pgb2->w_trapz_lens,
              pgb2->tau_size_bessel * sizeof(double*),
              ppt->error_message);

  class_call(array_trapezoidal_weights(pgb2->r,
                                       pgb2->r_size,
                                       w_trapz_r,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);

  class_call(array_trapezoidal_weights(pgb2->alpha,
                                       pgb2->alpha_size,
                                       w_trapz_alpha,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);

  class_call(array_trapezoidal_weights(pgb2->tau_sampling_bessel,
                                       pgb2->tau_size_bessel,
                                       pgb2->w_trapz_lens,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);


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
                                          tau0_minus_tau[bin],
                                          w_trapz[bin],
                                          pgb2->tau_size_selection,
                                          pvecback,
                                          tau0,
                                          bin),
               pgb2->error_message,
               pgb2->error_message);
  }

  printf("selection[0][max] = %g\n", selection[0][pgb2->tau_size_selection-1]);


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

  for (int bin = 0; bin < ppt->selection_num; bin++) {
    for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
      tau_window = pgb2->tau_sampling_selection[bin][index_tau];
      /* get background quantitites at this time */
      class_call(background_at_tau(pba,
                                   pgb2->tau_sampling_selection[bin][index_tau],
                                   pba->long_info,
                                   pba->inter_normal,
                                   &last_index_window,
                                   bac_window),
                 pba->error_message,
                 pgb2->error_message);

      z = pba->a_today/bac_window[pba->index_bg_a]-1.;

      selection_test[bin][index_tau] = exp(-0.5*pow(fabs(z-ppt->selection_mean[bin])/ppt->selection_width[bin],2))/(2.0*_PI_)/ppt->selection_width[bin];
      //printf("selection_test[%d][%d] = %g, should be \n", bin, index_tau, selection_test[bin][index_tau]/*, exp(-0.5*pow(fabs(z-ppt->selection_mean[bin])/ppt->selection_width[bin],2))/(2.0*_PI_)/ppt->selection_width[bin]*/);
    }
  }
  double test_sum_selection = 0.0;
  for (int bin = 0; bin < ppt->selection_num; bin++) {
    normal = 0.0;
    class_call(array_trapezoidal_integral(selection_test[bin],
                                          pgb2->tau_size_selection,
                                          w_trapz[bin],
                                          &normal,
                                          ptr->error_message),
               pgb2->error_message,
               pgb2->error_message);

    for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {



      selection_test[bin][index_tau] = selection_test[bin][index_tau]/normal;

      //printf("tau = %g, selection[%d][%d] = %g\n", pgb2->tau_sampling_selection[bin][index_tau], bin, index_tau, selection[bin][index_tau]);
      //printf("selection_test[%d] at tau = %g is = **%g**\n", index_tau, pgb2->tau_sampling_selection[bin][index_tau], selection_test[bin][index_tau]);
      //printf("(selection_test[%d][%d], selection[%d][%d]) = (%g,%g)\n", bin, index_tau, bin, index_tau,selection_test[bin][index_tau], selection[bin][index_tau]);
      test_sum += selection_test[bin][index_tau] * w_trapz[bin][index_tau];
      test_sum_selection += selection[bin][index_tau] *w_trapz[bin][index_tau];
    }
  }

  //exit(0);
  class_call(background_tau_of_z(
                          pba,
                          ppt->selection_mean[0],
                          &selection_mean),
                          ppt->error_message,
                          pgb2->error_message);

  //printf("selection mean tau = %g\n", selection_mean);
  printf("test_sum = %g\n",test_sum);
  printf("test_sum_selection = %g\n",test_sum_selection);
  //exit(0);
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
  printf("pgb2->tau_sampling_cls has %d points sampled between (%g,%g)\n", pgb2->tau_size_selection, pgb2->tau_sampling_cls[0], pgb2->tau_sampling_cls[pgb2->tau_size_selection-1] );
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
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
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
        for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
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
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
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
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
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
            for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
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
            for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
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
            for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
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
            for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){

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
            for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
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
            for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
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
            for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
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
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){

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
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){

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

        pgb2->first_order_sources[pgb2->index_source_phi_prime][pgb2->tau_size_selection-1][index_k_bessel] = (pgb2->first_order_sources[pgb2->index_source_phi][pgb2->tau_size_selection-1][index_k_bessel]-pgb2->first_order_sources[pgb2->index_source_phi][pgb2->tau_size_selection-2][index_k_bessel])/
          (pgb2->tau_sampling_cls[pgb2->tau_size_selection-1]-pgb2->tau_sampling_cls[pgb2->tau_size_selection-2]);


        pgb2->first_order_sources[pgb2->index_source_phi_plus_psi_prime][0][index_k_bessel] = (pgb2->first_order_sources[pgb2->index_source_phi_plus_psi][1][index_k_bessel]-pgb2->first_order_sources[pgb2->index_source_phi_plus_psi][0][index_k_bessel])/
          (pgb2->tau_sampling_cls[1]-pgb2->tau_sampling_cls[0]);

        pgb2->first_order_sources[pgb2->index_source_phi_plus_psi_prime][pgb2->tau_size_selection-1][index_k_bessel] = (pgb2->first_order_sources[pgb2->index_source_phi_plus_psi][pgb2->tau_size_selection-1][index_k_bessel]-pgb2->first_order_sources[pgb2->index_source_phi_plus_psi][pgb2->tau_size_selection-2][index_k_bessel])/
          (pgb2->tau_sampling_cls[pgb2->tau_size_selection-1]-pgb2->tau_sampling_cls[pgb2->tau_size_selection-2]);
      }


      for (int index_tau = 1; index_tau < pgb2->tau_size_selection-1; index_tau++){
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
        for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){

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
        class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l], pgb2->tau_size_selection * sizeof(double *), pgb2->error_message);
        for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++) {
          class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first], pgb2->tau_size_selection * sizeof(double), pgb2->error_message);
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
        class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l], pgb2->tau_size_selection * sizeof(double *), pgb2->error_message);
        for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++) {
          class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first], pgb2->tau_size_selection * sizeof(double), pgb2->error_message);
        }
      }
    }
  }*/
/* Allocate array for Dl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] */
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
  printf("Allocating size %ix%ix%ix%ix%i bytes \n", pgb2->type_size, pgb2->type_size, ptr->l_size[ppt->index_md_scalars], pgb2->tau_size_selection, pgb2->tau_size_selection);

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





  /* Fill the array pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] by calling the integral() function */
  for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    for (int k = 0; k < index_type_first+1; k++){
      printf("X");
    }
    printf("\n");
// alert last tau width is missing
  /* Allocate background vectors */
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
        //printf("\r          ");
        //printf("%d/%d\n",index_l, ptr->l_size[ppt->index_md_scalars]-1);
        //for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++){
        for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
          for(int index_r = 0; index_r < pgb2->r_size; index_r++){
            tau_one = pba->conformal_age*(1.0-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
            tau_two = pba->conformal_age*(1.0-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));

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
            //for(int index_tau_second = 0; index_tau_second < index_tau_first+1/*pgb2->tau_size_selection*/; index_tau_second++){
            //for(int index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++){
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
              type1 = pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_first][index_k_bessel];
              //printf("type1 = %g\n", type1);
              //printf("type2 = %g\n", type2);
              type2 = pgb2->first_order_sources_integ[index_type_second][index_l][index_tau_second][index_k_bessel];


              /* Using area of trapezoid, we can sum a number of areas of trapezoids to approximate the integral */
              integ += pow(pgb2->k_bessel[index_k_bessel],-1.0) * 4. * _PI_ *  Pk * type1 * type2  * pgb2->w_trapz_k[index_k_bessel];
            }
            pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r] = integ;
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
    //exit(0);

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
  printf("#r[%d] = %g\n", index_test_r, pgb2->r[index_test_r]);
  printf("#l[%d] = %d\n", index_test_l, ptr->l[index_test_l]);
  printf("#dens-rsd\n" );
  printf("#alpha     Dl");
  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    //tau_one = pba->conformal_age*(1.0-pgb2->r[index_r]*sin(pgb2->alpha[index_test_alpha]));
    //tau_two = pba->conformal_age*(1.0-pgb2->r[index_r]*cos(pgb2->alpha[index_test_alpha]));

    //index_of_tau_sampling_cls(tau_one, &index_tau_one_cls, pgb2);
    //index_of_tau_sampling_cls(tau_two, &index_tau_two_cls, pgb2);

    //tau_one_cls = pgb2->tau_sampling_cls[index_tau_one_cls];
    //tau_two_cls = pgb2->tau_sampling_cls[index_tau_two_cls];

    /*double result = 0.5*(pgb2->Dl2[pgb2->index_type_density][pgb2->index_type_rsd][index_test_l][index_alpha][index_test_r]
                          +pgb2->Dl2[pgb2->index_type_rsd][pgb2->index_type_density][index_test_l][index_alpha][index_test_r]);*/

    //double result =  pgb2->Dl2[pgb2->index_type_density][pgb2->index_type_density][index_test_l][index_alpha][index_test_r];
    //printf("%g      %g\n",pgb2->alpha[index_alpha], result);
  }
  //exit(0);
  // Warning: This only works for auto-term correlations index_type_first = index_type_second.
  /*for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    for(int index_type_second = 0; index_type_second < index_type_first+1; index_type_second++){
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
        for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++){
          for(int index_tau_second = index_tau_first+1; index_tau_second < pgb2->tau_size_selection; index_tau_second++){
            pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] = pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_second][index_tau_first];
          }
        }
      }
    }
  }*/

  for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++) {
    for (int index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++) {
      norm_test = ptr->l[test_index_l]*(ptr->l[test_index_l]+1.)/(2*_PI_);
      fprintf(file,
              "%g     %g      %g\n",
              pgb2->tau_sampling_cls[index_tau_first],
              pgb2->tau_sampling_cls[index_tau_second]);
    }
  }
  fprintf(file, "K_MAX = %g \n",pgb2->k_bessel[pgb2->k_size_bessel-1] );
  fclose(file);
  printf("Written Dl data into file.\n");
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


  /*else{
    printf("Entered non-Dirac window function integration.\n");
    double dummy_inner;
    double dummy_outer;

    for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
      for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
        for (int bin1 = 0; bin1 < 1/*ppt->selection_num*///; bin1++) {
        /*  for (int bin2 = 1; bin2 < 2/*ppt->selection_num*///; bin2++) {
          //  for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){

            //  pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2] = 0.;
              //dummy_outer = 0.0;
              //double temp456 = 0.0;
      /* Final result Cl[index_type_first][index_type_second][index_l][bin1][bin2] is defined on the pgb2->tau_sampling_cls, so an interpolation is required
            from Cl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] which is defined on pgb2->tau_size_selection[bin][index_tau]*/
              /*for(index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++){

                double temp123 = 0.;
                dummy_inner = 0.0;
                double tau2 = pgb2->tau_sampling_selection[bin2][index_tau_second];
                index_of_tau_sampling_cls(tau2, &index_tau2, pgb2);


                for(index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++){
                  //printf("indexing %dx%dx%dx%dx%dx%dx%d\n",index_type_first, index_type_second, bin1, bin2, index_l,index_tau_second,index_tau_first);

                  double tau1 = pgb2->tau_sampling_selection[bin1][index_tau_first];

                  index_of_tau_sampling_cls(tau1, &index_tau1, pgb2);

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

                  dummy_inner += w_trapz[bin1][index_tau_first]*selection[bin1][index_tau_first];

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
  }*/


  /* alpha and r time parameterisation integration */
  else{
    printf("Entered non-Dirac window function integration.\n");
    printf("Computing Cls via r and alpha integration...\n");

    double test_integ;
    double temp456;
    double temp123;
    printf("#l    r   int tau0*tau0*r*Dl da dr\n");
    for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
      for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
        //for (int bin1 = 0; bin1 < 1/*ppt->selection_num*/; bin1++) {
          //for (int bin2 = 1; bin2 < 2/*ppt->selection_num*/; bin2++) {
            for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){

              //pgb2->Cl[index_type_first][index_type_second][index_l][bin1][bin2] = 0.;

              temp456 = 0.0;

      /* Final result Cl[index_type_first][index_type_second][index_l][bin1][bin2] is defined on the pgb2->tau_sampling_cls, so an interpolation is required
            from Cl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] which is defined on pgb2->tau_size_selection[bin][index_tau]*/
              for (int index_r_up = 1; index_r_up < pgb2->r_size; index_r_up++) {
                double result = 0.0;


                for(int index_r = 0; index_r < index_r_up+1; index_r++){
                  temp123 = 0.0;
                  test_integ = 0.0;

                  for(int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++){

                    double tau1 = pba->conformal_age*(1.-pgb2->r2[index_alpha][index_r]*sin(pgb2->alpha[index_alpha]));
                    double tau2 = pba->conformal_age*(1.-pgb2->r2[index_alpha][index_r]*cos(pgb2->alpha[index_alpha]));

                    //index_of_tau_sampling_selection(tau1, bin1, &index_tau_first, pgb2);
                    //index_of_tau_sampling_selection(tau2, bin2, &index_tau_second, pgb2);

                    test_integ +=pba->conformal_age*pba->conformal_age*pgb2->r2[index_alpha][index_r]
                                    *pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r]
                                      *w_trapz_alpha[index_alpha];

                    //printf("pgb2->Dl2 = %g\n", pgb2->Dl2[pgb2->index_type_density][pgb2->index_type_density][index_l][index_alpha][index_r]);
                    /*temp123 += pba->conformal_age*pba->conformal_age*pgb2->r[index_r]* pgb2->Dl2[index_type_first][index_type_second][index_l][index_alpha][index_r]
                                * w_trapz_r[index_r] * selection[bin1][index_tau_first]*selection[bin2][index_tau_second];*/


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
  printf("Reached here\n");
  exit(0);

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
           norm*pgb2->Cl[pgb2->index_type_rsd][pgb2->index_type_rsd][index_l][0][1]);
    }

  exit(0);


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
