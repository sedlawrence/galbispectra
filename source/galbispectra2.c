/* Module to compute galaxy bispectra */

#include <math.h>
#include "galbispectra2.h"
#include "perturbations2.h"
#define getName(var)  #var

/*==================================================================
====================================================================
=============================   TOOLS ==============================
====================================================================
===================================================================*/

double heaviside(double x,
                 double * result){

                  if (x < 0.0) {
                    *result = 0.0;
                  }
                  if (x == 0.0) {
                    *result = 0.5;
                  }
                  if(x > 0.0){
                    *result = 1.0;
                  }

}

Al1l2l3(int l2,
        int l3,
        double * result,
        struct galbispectra2 * pgb2){

  double minI, minII, minIII;
  double maxI, maxII, maxIII;
  double AI, AII, AIII;


  class_call(drc3jj (l2,
                     l3,
                     1,
                     -1,
                     &minI,
                     &maxI,
                     &AI,
                     1000,
                     pgb2->error_message),
                     pgb2->error_message,
                     pgb2->error_message);

  class_call(drc3jj (l2,
                     l3,
                     -1,
                     1,
                     &minII,
                     &maxII,
                     &AII,
                     1000,
                     pgb2->error_message),
                     pgb2->error_message,
                     pgb2->error_message);

  class_call(drc3jj (l2,
                     l3,
                     0,
                     0,
                     &minIII,
                     &maxIII,
                     &AIII,
                     1000,
                     pgb2->error_message),
                     pgb2->error_message,
                     pgb2->error_message);


  double A = 0.5*(AI+AII)/AIII;

  *result = A;

}



/*********************** SPHERICAL BESSEL FUNCTIONS *****************/

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
      *result = 0.0;
      out = 0.0;
    }

    // NOTE THIS NEEDS TO BE FIXED
    /*else if(z< 1e-4){

      class_call(bessel_at_x(pbs, z, index_l, &j1), pbs->error_message, pgb2->error_message);

      class_call(bessel_at_x_first_deriv(pgb2,
                             pbs,
                             z,
                             index_l,
                             &j2),
                             pbs->error_message,
                             pgb2->error_message);
      out = (-2.*j2)/z-j1;
    }*/
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

      A = (-2.*j2)/z;
      B = -j1;
      C = (l*(l+1)*j1)/(z*z);
      l = pbs->l[index_l];
      /* The following out works */
      out = (-2*z*j2-(z*z-l*(l+1))*j1)/(z*z);
      /* Modified */
    //  out = (-2.*j2)/z-j1+(l*(l+1)*j1)/(z*z);
      //out = A+B+C;


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


double bessel_at_x_third_deriv(struct galbispectra2 * pgb2,
                       struct bessels * pbs,
                       double x,
                       int index_l,
                       double * result){

    double out;
    double jl;
    double jlplus1;

    int l = pbs->l[index_l];

    if (x == 0.0) {
      *result = 0.0;
    }

    else{
     /*using recurrence relation https://dlmf.nist.gov/10.51 eq. 10.51.2 */
      class_call(bessel_j(pbs, pbs->l[index_l], x, &jl), pbs->error_message, pgb2->error_message);
      class_call(bessel_j(pbs, pbs->l[index_l]+1, x, &jlplus1), pbs->error_message, pgb2->error_message);

      out = (l-2)
            *(l*l-l-x*x)
            *jl
            /x
            /x
            /x
            -(l*l+l-x*x+6)
            *jlplus1
            /x
            /x;

      *result = out;
    }

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
    printf("search out of bounds!\n");
    printf("%s\n",'Search out of bounds in index_of_k()' );
  }


  return _SUCCESS_;
}

int index_of_l(int l,
               int * index_l,
               int * last_index_l,
               struct perturbs * ppt,
               struct transfers * ptr){
  double l_start;
  double l_end;

  /* Initialise index_k to a value that does not exist within the grid, this will circumvent any chance of the
    wrong index being assigned */
  *index_l = -1;
  //printf("inside iol: last_index_l = %d\n", *last_index_l);
  if (*last_index_l > 0){

    l_start = ptr->l[*last_index_l - 1];
    l_end = ptr->l[*last_index_l];

    if(l_end == l){
      *index_l=*last_index_l;
      //printf("1. setting *index_l = *last_index_l\n" );
    }
  }

  else{
    l_start = ptr->l[* last_index_l];
  }

  double l_ptr;


  /* If l< l_start, search grid from index = 0 */
  if (l < l_start){
    for (int index = 0; index < *last_index_l; index++) {
      l_ptr = ptr->l[index];

      if (l < l_ptr && *index_l == -1) {

        *index_l=index;
        *last_index_l=index;
        break;
      }

    }
  }

  else {

    for (int index = *last_index_l; index < ptr->l_size[ppt->index_md_scalars]; index++) {
     l_ptr = ptr->l[index];

     if (l <= l_ptr && *index_l == -1) {

       *index_l=index;
       *last_index_l=index;
       break;
     }
    }
  }

  if (l == ptr->l[*last_index_l]) {

    *index_l=*last_index_l;

  }

  if (*index_l == -1){
    printf("search out of bounds!\n");
    printf("%s\n",'Search out of bounds in index_of_l()' );
  }


  return _SUCCESS_;
}

int interpolate_Dl_for_l(int l,
               int * index_type_first,
               int * index_type_second,
               int * bin1,
               int * bin2,
               int * last_index_l,
               int * index_tau_first,
               int * index_tau_second,
               double * spline,
               struct galbispectra2 * pgb2,
               struct transfers * ptr,
               struct perturbs * ppt){

  int exact_found = -1;
  int index_l_found;
  //printf("last_index_l = %d\n", last_index_l);
  double l_start;
  double l_end;
  int index_l;

  /* Search for closest index_l that gives a corresponding l_found >= l_find. Initialise index_l to a value that does not exist within the grid, this will circumvent any chance of the
    wrong index being assigned */
  index_l = -1;
  //printf("inside iol: last_index_l = %d\n", *last_index_l);
  if (*last_index_l > 0){

    l_start = ptr->l[*last_index_l - 1];
    l_end = ptr->l[*last_index_l];

    if(l_end == l){
      index_l=*last_index_l;
      //printf("1. setting index_l = *last_index_l\n" );
    }
  }

  else{
    l_start = ptr->l[* last_index_l];
  }

  double l_ptr;


  /* If l< l_start, search grid from index = 0 */
  if (l < l_start){
    for (int index = 0; index < *last_index_l; index++) {
      l_ptr = ptr->l[index];

      if (l < l_ptr && index_l == -1) {
      //  printf("1. setting index_l = *index\n" );
        index_l=index;
        *last_index_l=index;
        break;
      }

    }
  }

  else {

    for (int index = *last_index_l; index < ptr->l_size[ppt->index_md_scalars]; index++) {
     l_ptr = ptr->l[index];

     if (l <= l_ptr && index_l == -1) {
       //printf("2. setting index_l = *index\n" );
       index_l=index;
       *last_index_l=index;
       break;
     }
    }
  }

  if (l == ptr->l[*last_index_l]) {
    //printf("exact found: l = %d, ptr->l[*last_index_l] = %d\n", l, ptr->l[*last_index_l]);
  //  printf("setting index_l = *last_index_l\n" );
    exact_found = 1;
    index_l=*last_index_l;
    //printf(" index_l = %d\n", index_l );

  }

  if (index_l == -1){
    //printf("search out of bounds!\n");
    //printf("%s\n",'Search out of bounds in index_of_l()' );
  }

  //printf("end of search: index_l = %d\n", index_l);
  //printf("l = %d closest = %d\n", l, ptr->l[index_l_found]);

  if (exact_found == 1) {
    //printf("entering exact found\n");
    *spline = pgb2->Dl[*index_type_first][*index_type_second][index_l][*bin1][*bin2][*index_tau_first][*index_tau_second];

  }

  else{
    *spline = pgb2->Dl[*index_type_first][*index_type_second][index_l-1][*bin1][*bin2][*index_tau_first][*index_tau_second]
      *(ptr->l[index_l]-l)+pgb2->Dl[*index_type_first][*index_type_second][index_l][*bin1][*bin2][*index_tau_first][*index_tau_second]
        *(l-ptr->l[index_l-1]);
    *spline /= 1.0*(ptr->l[index_l] - ptr->l[index_l-1]);
  }
  /*double spline = restrict_array[index_l_found-1]*(ptr->l[index_l_found]-l)
                + restrict_array[index_l_found]*(l-ptr->l[index_l_found-1]);*/


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


  if(k<k_start){printf("inside index_of_k: THING THAT SHOULD NOT HAPPEN HAPPENED k =%g < k_start = %g, index_k = %d\n", k, k_start, index_k);exit(2);}
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
               int * index_k,
               struct galbispectra2 * pgb2){
  double k_start;
  double k_pgb2;

  k_start = pgb2->k_bessel[*index_k];

  //printf("k_required = %g\n",k );
  //printf("k_start = %g\n",k_start );
  if(k<k_start){printf("inside index_of_k_bessel: THING THAT SHOULD NOT HAPPEN HAPPENED\n");exit(2);}
  for (int index = *index_k+1; index < pgb2->k_size_bessel; index++) {
   k_pgb2 = pgb2->k_bessel[index];
   //printf("k_pgb2 scan: %g\n", k_pgb2 );

   if (k_pgb2>k) {
     *index_k=index-1;
     break;
   }

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

                   if(tau<tau_start){printf("index_of_tau_sampling_quadsources: THING THAT SHOULD NOT HAPPEN HAPPENED\n");exit(2);}
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

  /*==================================================================
  ====================================================================
  =============================   SAMPLING  ==========================
  ====================================================================
  ===================================================================*/
  printf("Starting sampling..\n" );




    /* Define local (stack) variables */
  int index_l;
  int index_k;
  int index_tau;
  int index_tau_ppt2;
  int index_type;
  int index_md;
  int bin_first;
  int bin_second;
  double * pvecback;
  double * w_trapz_tau;
  double tau0;
  double ** w_trapz;
  double ** tau0_minus_tau;
  double ** tau0_minus_tau_lens_bi;
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

  /* Define the SIZES/resolution of grids */

  /*the bessel grid is boosted by an integer factor to account for rapid oscillations.*/
  int bessel_boost = 1;
  double e = 2.71828;
  double g_bias = 1.0;
  tau0 = pba->conformal_age;
  double heaviside_test;
  double ** A;
  class_alloc(A,
              10 * sizeof(double*),
              ppt->error_message);
  for (int i = 0; i < 10; i++) {
    class_alloc(A[i],
                10 * sizeof(double),
                ppt->error_message);

    for (int j = 0; j < 10; j++) {
      A[i][j] = i+j;

    }
  }


  for (int i = 0; i < 10; i++) {
    free(A[i]);
  }
  free(A);

  printf("A[1][5] = %g\n", A[1][5]);


  printf("tau0 = %g\n",tau0);




  /* Please ensure the pgb2->tau_size_selection grid is of a higher resolution than pgb2->r_size * pgb2->alpha_size = 13*/
  pgb2->tau_size_cls = 500; //prev on 500
  pgb2->tau_size_selection = 500; //prev on 601



  /* the tau-bessel grid is a hi-res time grid used for the lening convergence term. We wish to double the number of slots such that points in the two grids tau_sampling_cls
   and tau_sampling_bessel align. */

  pgb2->tau_size_bessel = bessel_boost * (pgb2->tau_size_cls-1) + 1;
  pgb2->k_size_bessel = 2500; // prev 2500 4999 /*tau_size_bessel*/   /*k_size * boost*/ // formerly on 2500
  printf("Starting galaxy bispectra module...\n");
  if (ppt->selection==gaussian) {
    printf("We have a Gaussian Window Function on our hands.\n");
  }
  if (ppt->selection==dirac) {
    printf("We have a Dirac Window Function on our hands.\n");
  }

  double * result;



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

  class_alloc(tau0_minus_tau_lens_bi,
              ppt->selection_num * sizeof(double*),
              ppt->error_message);

  for ( bin = 0; bin < ppt->selection_num; bin++) {

    class_alloc(tau0_minus_tau_lens_bi[bin],
              pgb2->tau_size_selection * sizeof(double),
              ppt->error_message);
  }

  double *** w_trapz_lens_bi;

  class_alloc(pgb2->tau_sampling_selection,
              ppt->selection_num * sizeof(double*),
              pgb2->error_message);

  class_alloc(pgb2->tau_sampling_lens_bi,
              ppt->selection_num * sizeof(double**),
              pgb2->error_message);

  class_alloc(w_trapz_lens_bi,
              ppt->selection_num * sizeof(double**),
              pgb2->error_message);


  for (bin = 0; bin < ppt->selection_num; bin++) {

    class_alloc(pgb2->tau_sampling_selection[bin],
                pgb2->tau_size_selection * sizeof(double),
                pgb2->error_message);

    class_alloc(pgb2->tau_sampling_lens_bi[bin],
                pgb2->tau_size_selection * sizeof(double*),
                pgb2->error_message);

    class_alloc(w_trapz_lens_bi[bin],
                pgb2->tau_size_selection * sizeof(double*),
                pgb2->error_message);

    for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
      class_alloc(pgb2->tau_sampling_lens_bi[bin][index_tau],
                  pgb2->tau_size_selection * bessel_boost * sizeof(double),
                  pgb2->error_message);

      class_alloc(w_trapz_lens_bi[bin][index_tau],
                  bessel_boost*(index_tau+1) * sizeof(double),
                  ppt->error_message);
    }
  }

  class_alloc(pgb2->tau_sampling_cls,
              pgb2->tau_size_cls * sizeof(double),
              pgb2->error_message);






  double overall_tau_min = 160000.;
  double overall_tau_max = -1.;
  int * bin_mean_index_selection;
  class_alloc(bin_mean_index_selection,
              ppt->selection_num * sizeof(int),
              pgb2->error_message);

  int * bin_mean_index_cls;
  class_alloc(bin_mean_index_cls,
              ppt->selection_num * sizeof(int),
              pgb2->error_message);

  /* Define new tau_sampling_selection. This sampling is bin-dependent */
  double selection_mean;
  for (bin = 0; bin < ppt->selection_num; bin++) {
    //finer sampling of bins
    double tau_min;
    double tau_max;
    double z_max, z_min;
    int index_selection_mean;



     //printf("ppt->selection_mean[%d]=%g, width = %g\n",bin,ppt->selection_mean[bin], ppt->selection_width[bin] );
     // ALERT HARDCODED
    z_max =  ppt->selection_mean[bin]+ 5. * ppt->selection_width[bin];
    z_min =  /*ppt->selection_mean[bin]- 5. * ppt->selection_width[bin];//*/ 0.001;
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
    printf("selection_mean = %g\n", selection_mean );

    if (tau_max > ppt->tau_sampling_quadsources[ppt->tau_size_quadsources-1]) {
      tau_max = ppt->tau_sampling_quadsources[ppt->tau_size_quadsources-1];
    }



    if (tau_min < 400.) {
      printf("Rubbish, Window function extends to before recombination\n");
      tau_min = 400.;
    }

    overall_tau_min = MIN(tau_min,overall_tau_min);
    overall_tau_max = MAX(tau_max,overall_tau_max);
    printf("selection_mean = %g\n", selection_mean );
    for (index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
      /* Warning Dirac-case */
      //pgb2->tau_sampling_selection[bin][index_tau] = selection_mean;
      pgb2->tau_sampling_selection[bin][index_tau] = tau_min + index_tau*(tau_max-tau_min)/(pgb2->tau_size_selection-1);



    }
      index_of_tau_sampling_selection(selection_mean,
                    bin,
                    &index_selection_mean,
                    pgb2);
      pgb2->tau_sampling_selection[bin][index_selection_mean] = selection_mean;

      bin_mean_index_selection[bin] = index_selection_mean;
  }






  // ALERT HARDCODED!!
  for (index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {
    /* Warning Dirac-case */

    pgb2->tau_sampling_cls[index_tau] = overall_tau_min + index_tau*(overall_tau_max-overall_tau_min)/(pgb2->tau_size_cls-1);
    //pgb2->tau_sampling_cls[index_tau] = selection_mean;
    //printf("pgb2->tau_sampling_cls[%d] = %g\n" , index_tau, pgb2->tau_sampling_cls[index_tau]);
    //printf("pgb2->tau_sampling_selection[%d] = %g\n" , index_tau, pgb2->tau_sampling_selection[0][index_tau]);

  }


  for (bin = 0; bin < ppt->selection_num; bin++) {
    int index_selection_mean;
    class_call(background_tau_of_z(
                            pba,
                            ppt->selection_mean[bin],
                            &selection_mean),
                            ppt->error_message,
                            pgb2->error_message);

    printf("selection_mean[%d] = %g\n", selection_mean);

    index_of_tau_sampling_cls(selection_mean,
                              &index_selection_mean,
                              pgb2);

    bin_mean_index_cls[bin] = index_selection_mean;

    printf("bin_mean_index_cls[%d] = %d\n", bin, bin_mean_index_cls[bin]);

    pgb2->tau_sampling_cls[index_selection_mean] = selection_mean;
  }






  /* Find the ranges of alpha to span such that we sample within tau_sampling_cls[0] to tau_sampling_cls[pgb2->tau_size_selection-1]*/

  printf("tau_sampling_cls has %d points that span (%g, %g)\n", pgb2->tau_size_cls, pgb2->tau_sampling_cls[0], pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]);
  //printf("tau_sampling_bessel has %d points that span (%g, %g)\n", pgb2->tau_size_bessel, pgb2->tau_sampling_bessel[index_tau_bessel], pgb2->tau_sampling_cls[pgb2->tau_size_bessel-1]);

  printf("conformal_age = %g\n",pba->conformal_age);







  printf("tau_max_cls = %g\n", pgb2->tau_sampling_cls[pgb2->tau_size_selection-1] );
  /* Allocate and define the pgb2->tau_sampling_bessel2[index_tau_cls][index_tau_selection] array which is used to integrate the innermost (time/chi) integral of the lensing term.
    There is bessel_boost times more resolution in this grid compared with tau_sampling cl. The bessel function is called using the bessel2 grid,
    so it is not necessary to have as high resolution in tau_cls compared with most other first order terms (as long as the bessel_boost makes the
    resolution in the tau_sampling_bessel2 grid sufficient.). pgb2->tau_sampling_cls[index_tau] is
    equivalent to pgb2->tau_sampling_bessel2[index_tau][bessel_boost*index_tau]*/


  //herehere
  double ** w_trapz_lens_bessel2;
  class_alloc(pgb2->tau_sampling_bessel2,
              pgb2->tau_size_cls * sizeof(double*),
              ppt->error_message);

  class_alloc(w_trapz_lens_bessel2,
              pgb2->tau_size_cls * sizeof(double*),
              ppt->error_message);



  class_alloc(pgb2->tau_sampling_lens_bi,
              pgb2->tau_size_selection * sizeof(double*),
              ppt->error_message);



  for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {

    class_alloc(pgb2->tau_sampling_bessel2[index_tau],
                bessel_boost*(index_tau+1) * sizeof(double),
                ppt->error_message);
    class_alloc(w_trapz_lens_bessel2[index_tau],
                bessel_boost*(index_tau+1) * sizeof(double),
                ppt->error_message);


  }















  for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(0+1); index_tau_bessel++) {

    pgb2->tau_sampling_bessel2[0][index_tau_bessel] = pgb2->tau_sampling_cls[0];

  }


  for (int index_tau = 1; index_tau < pgb2->tau_size_cls-1; index_tau++) {

    for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(index_tau+1); index_tau_bessel++) {
      double no_of_wedges = bessel_boost*(index_tau+1.)-1.;

      pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel] = pgb2->tau_sampling_cls[index_tau] + index_tau_bessel*((overall_tau_max/*tau0-0.001*/)-pgb2->tau_sampling_cls[index_tau])/(no_of_wedges);
      //printf("pgb2->tau_sampling_bessel2[%d][%d] = %g\n", index_tau, index_tau_bessel, pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel]);


    }
    //printf("pgb2->tau_sampling_selection[0][%d] = %g\n", index_tau, pgb2->tau_sampling_selection[0][index_tau]);
    //printf("pgb2->tau_sampling_bessel2[%d][%d] = %g\n", index_tau, 0, pgb2->tau_sampling_bessel2[index_tau][0]);
  }



  for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(index_tau+1); index_tau_bessel++) {
    pgb2->tau_sampling_bessel2[pgb2->tau_size_cls-1][index_tau_bessel] = overall_tau_max;
  }


  for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {
    class_call(array_trapezoidal_weights(pgb2->tau_sampling_bessel2[index_tau],
                                         bessel_boost*(index_tau+1),
                                         w_trapz_lens_bessel2[index_tau],
                                         pgb2->error_message),
                                         pgb2->error_message,
                                         pgb2->error_message);

  }




  /* Bessel k-sampling to capture features of the Bessel oscillations */

  class_alloc(pgb2->k_bessel,
              pgb2->k_size_bessel * sizeof(double),
              pgb2->error_message);
  printf("Allocated k_bessel array\n");


  /* Choose k-sampling: linear is zero, log is 1. (Currently set by hand). */
  int k_sampling = 0;

  double k_min = ppt->k[ppt->index_md_scalars][0];

  double k_max = ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1]; //0.0494957; //0.63613;// ;
  printf("k_max = %g\n", k_max );






  if (k_max < ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1]) {
    printf("ERROR! SONG doesn't comput k_max as high as galbispectra2.c demands. Please increase k_max_for_pk or reduce k_max in galbispectra2.c\n" );
    exit(0);
  }

  /* Linear case */
  if (k_sampling == 0) {
    for (int i = 0; i < pgb2->k_size_bessel; i++) {
      pgb2->k_bessel[i] = k_min + i*(k_max-k_min)/(pgb2->k_size_bessel-1);

    }
  }
  printf("pgb2->k_bessel[%d] = %g\n", pgb2->k_size_bessel-1, pgb2->k_bessel[pgb2->k_size_bessel-1]);


  /* Log case */
  if (k_sampling == 1 ) {
    double k_logbase = 1.01;
    for (int i = 0; i < pgb2->k_size_bessel; i++) {
      pgb2->k_bessel[i] =  k_min+pow(k_logbase,i)*(k_max-k_min)/(pow(k_logbase,pgb2->k_size_bessel-1));
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



  /* Allocate and fill array for the trapezoidal weights for chi integration in the lensing term w_trapz_lens[index_tau] */
  /*class_alloc(pgb2->w_trapz_lens,
              pgb2->tau_size_bessel * sizeof(double*),
              ppt->error_message);


  class_call(array_trapezoidal_weights(pgb2->tau_sampling_bessel,
                                       pgb2->tau_size_bessel,
                                       pgb2->w_trapz_lens,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);*/


   class_alloc(w_trapz,
               ppt->selection_num * sizeof(double*),
               ppt->error_message);

   for ( bin = 0; bin < ppt->selection_num; bin++) {
     class_alloc(w_trapz[bin],
                 pgb2->tau_size_selection * sizeof(double),
                 ppt->error_message);
   }

   class_alloc(w_trapz_lens_bi,
               ppt->selection_num * sizeof(double*),
               ppt->error_message);

   for ( bin = 0; bin < ppt->selection_num; bin++) {
     class_alloc(w_trapz_lens_bi[bin],
                 pgb2->tau_size_selection * sizeof(double),
                 ppt->error_message);
   }
  /* Fill the array ptw->tau0_minus_tau[index_tau] */
  int last_index_test;
  class_alloc(pvecback,
              pba->bg_size * sizeof(double),
              ptr->error_message);

  for (bin = 0; bin < ppt->selection_num ; bin++) {

    for(int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++){
      class_call(background_at_tau(pba,
                                   pgb2->tau_sampling_selection[bin][index_tau],
                                   pba->long_info,
                                   pba->inter_normal,
                                   &last_index_test,
                                   pvecback),
                                   pba->error_message,
                                   pgb2->error_message);

        /*infer redhsift */
      double z_test = pba->a_today/pvecback[pba->index_bg_a]-1.;
      tau0_minus_tau[bin][index_tau] = tau0 - pgb2->tau_sampling_selection[bin][index_tau];
      tau0_minus_tau_lens_bi[bin][index_tau] = (tau0 - pgb2->tau_sampling_selection[bin][0]) + index_tau*(tau0-(tau0 - pgb2->tau_sampling_selection[bin][0]))/(pgb2->tau_size_selection-1);

      //printf("%g      %g       %g\n", z_test, pgb2->tau_sampling_selection[bin][index_tau], tau0 - pgb2->tau_sampling_selection[bin][index_tau] );
      //printf("tau0 minus_tau_lens_bi[%d][%d] = %g\n", bin, index_tau, tau0_minus_tau_lens_bi[bin][index_tau]);
      //printf("tau0_minus_tau_lens_bi[%d][%d] = %g\n", bin, index_tau, tau0_minus_tau_lens_bi[bin][index_tau]);
    }


    class_call(array_trapezoidal_mweights(tau0_minus_tau[bin]/*pgb2->tau_sampling_selection[bin]*/,
                                          pgb2->tau_size_selection,
                                          w_trapz[bin],
                                          pgb2->error_message),
                                          ppt2->error_message,
                                          ppt2->error_message);

    class_call(array_trapezoidal_mweights(tau0_minus_tau_lens_bi[bin]/*pgb2->tau_sampling_selection[bin]*/,
                                          pgb2->tau_size_selection,
                                          w_trapz_lens_bi[bin],
                                          pgb2->error_message),
                                          ppt2->error_message,
                                          ppt2->error_message);
  }



  /* Computing the window-function */
  /* Declaration of temporary pointer */
  double ** selection;

  /* Allocation of first dimension selection[bin] */
  class_alloc(selection,
              ppt->selection_num * sizeof(double*),
              ppt->error_message);

  printf("selection_num = %d\n",ppt->selection_num);



  /* Allocation of second dimension selection[bin][index_tau] */
  for(int bin = 0; bin < ppt->selection_num; bin++){
    printf("bin = %d\n", bin );
    class_alloc(selection[bin],
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




  /* Absolute maximum value r can have is when alpha =pi/4 and tau is set to tau_min. If r1 or r2 are above this value, this is because of a sing-
  ularity in the denominator */
  //double tau_one;
  double tau_two;
  double tau_one_cls;
  double tau_two_cls;
  int index_tau_one_cls;
  int index_tau_two_cls;


  /* define the minimum distance between tau0-tau0 and the tau1_upper_cut_off-tau2_upper_cut_off, for each tau1-tau2 pair;
  that means permuting over bins. */


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
  pgb2->index_type_delta = -1;
  pgb2->index_type_lens = -1;
  pgb2->index_type_d1 = -1;
  pgb2->index_type_d2 = -1;
  pgb2->index_type_g1 = -1;
  pgb2->index_type_g2 = -1;
  pgb2->index_type_g3 = -1;
  pgb2->index_type_g4 = -1;
  pgb2->index_type_g5 = -1;
  pgb2->index_type_quad_density_p = -1;
  pgb2->index_type_quad_v = -1;
  pgb2->index_type_quad_v_p = -1; /*p denotes prime. The quad terms are linear but only appear at second order with another quad source */
  pgb2->index_type_quad_v_pp = -1;

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
  if (k5 == 5.0){
    pgb2->index_source_delta_cdm = index_source;
    index_source++;
    pgb2->index_type_density = index_type;
    index_type++;
    //pgb2->index_type_delta = index_type; /* Eq. 3.28 in [1812.09297] */
    //index_type++;

    //pgb2->index_type_quad_density_p = index_type;
    //index_type++;
    //pgb2->index_source_phi = index_source;
    //index_source++;
    //pgb2->index_source_v = index_type;
    //index_source++;
    //pgb2->index_source_phi = index_type;
    //index_source++;*/
  }
  // turn on for rsd

 //if (k5 == 5.0){
    //pgb2->index_source_v = index_source;
    //index_source++;
    //pgb2->index_type_vp = index_source;
    //index_type++;
    //pgb2->index_source_theta = index_source;
    //index_source++;

    //pgb2->index_type_rsd = index_type;
    //index_type++;

    //pgb2->index_type_quad_v = index_type;
    //index_type++;

    //pgb2->index_type_quad_v_pp = index_type;
    //index_type++;
    //pgb2->index_type_d1 = index_type;
    //index_type++;
    //pgb2->index_type_d2 = index_type;
    //index_type++;
  //}

  if (k5 == 5.0){
    pgb2->index_source_psi = index_source;
    index_source++;
    pgb2->index_source_phi = index_source;
    index_source++;
    //pgb2->index_source_phi_prime = index_source;
    //index_source++;
    pgb2->index_source_phi_plus_psi = index_source;
    index_source++;
    //pgb2->index_source_phi_plus_psi_prime = index_source;
    //index_source++;



    //pgb2->index_type_g1 = index_type;
    //index_type++;
    //pgb2->index_type_g2 = index_type;
    //index_type++;
    //pgb2->index_type_g3 = index_type;
    //index_type++;

    //pgb2->index_type_g5 = index_type;
    //index_type++;

    pgb2->index_type_lens = index_type;
    index_type++;

    pgb2->index_type_g4 = index_type;
    index_type++;
  }



  /*if (k5 == 5.0){
    pgb2->index_source_psi = index_source;
    index_source++;
    pgb2->index_source_phi = index_source;
    index_source++;
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

  /* Define an array of values of first order transfer functions:
              pgb2->second_order_sources_eq[index_type][index_tau][index_k_bessel] */
  /*class_alloc(pgb2->second_order_sources_eq, pgb2->source_size * sizeof(double **), ppt->error_message);
    for (int index_type = 0; index_type < pgb2->source_size; index_type++) {
      /* Allocate memory for pgb2->second_order_sources_eq[index_type][index_tau] */
    /*  class_alloc(pgb2->second_order_sources_eq[index_type],
                  pgb2->tau_size_cls * sizeof(double *),
                  ppt->error_message);
      /* Allocate memory for pgb2->second_order_sources_eq[index_type] */
    /*  for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {
          /* Loop over type and tau. For each of them, allocate memory
           for pgb2->second_order_sources_eq[index_type][index_tau][index_k_bessel]  */
    /*    class_alloc(pgb2->second_order_sources_eq[index_type][index_tau],
                    pgb2->k_size_bessel * sizeof(double),
                    ppt->error_message);
    }
  }*/

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



    double *** w_trapz_selection_hires;

    class_alloc(pgb2->tau_sampling_selection_hires,
                ppt->selection_num * sizeof(double**),
                ppt->error_message);

    class_alloc(w_trapz_selection_hires,
                ppt->selection_num * sizeof(double**),
                ppt->error_message);

    for ( bin = 0; bin < ppt->selection_num; bin++) {

      class_alloc(pgb2->tau_sampling_selection_hires[bin],
                  pgb2->tau_size_selection * sizeof(double*),
                  ppt->error_message);

      class_alloc(w_trapz_selection_hires[bin],
                  pgb2->tau_size_selection * sizeof(double*),
                  ppt->error_message);
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {

        class_alloc(pgb2->tau_sampling_selection_hires[bin][index_tau],
                    bessel_boost*((pgb2->tau_size_selection-1)-index_tau+1) * sizeof(double),
                    ppt->error_message);

        class_alloc(w_trapz_selection_hires[bin][index_tau],
                    bessel_boost*((pgb2->tau_size_selection-1)-index_tau+1) * sizeof(double),
                    ppt->error_message);


      }
    }
      /* fine here */




    /* We next deine pgb2->tau_sampling_selection_hires, which is a higher resolution version of tau_sampling_selection
      that runs between, tau and tau_max for each bin and each tau value in the tau_sampling_selection grid. This is useful for any
      integrals over this range. */

    for ( bin = 0; bin < ppt->selection_num; bin++) {
      for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*((pgb2->tau_size_selection-1)-0+1) /*bessel_boost*(0+1)*/; index_tau_bessel++) {

        pgb2->tau_sampling_selection_hires[bin][0][index_tau_bessel] = pgb2->tau_sampling_selection[bin][0]+index_tau_bessel*((pgb2->tau_sampling_selection[bin][pgb2->tau_size_selection-1])-pgb2->tau_sampling_selection[bin][0])/(bessel_boost*((pgb2->tau_size_selection-1)-0.+1.)-1.);
        //printf("pgb2->tau_sampling_selection_hires[%d][%d][%d] = %g\n", bin, 0, index_tau_bessel, pgb2->tau_sampling_selection_hires[bin][0][index_tau_bessel]);
      }




      /*selection is fine at this point*/
      for (int index_tau = 1; index_tau < pgb2->tau_size_selection; index_tau++) {
        //printf("index_tau = %d\n", index_tau);
        //printf("##2. k_min = %g\n", ppt->k[ppt->index_md_scalars][0]);

        for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*((pgb2->tau_size_selection-1)-index_tau+1); index_tau_bessel++) {
          double no_of_wedges = bessel_boost*((pgb2->tau_size_selection-1)-index_tau+1.)-1.;

          pgb2->tau_sampling_selection_hires[bin][index_tau][index_tau_bessel] = pgb2->tau_sampling_selection[bin][index_tau] + index_tau_bessel*((pgb2->tau_sampling_selection[bin][pgb2->tau_size_selection-1])-pgb2->tau_sampling_selection[bin][index_tau])/(no_of_wedges);
          printf("pgb2->tau_sampling_selection_hires[%d][%d][%d] = %g\n", bin, index_tau, index_tau_bessel, pgb2->tau_sampling_selection_hires[bin][index_tau][index_tau_bessel]);


        }
        //printf("pgb2->tau_sampling_selection[0][%d] = %g\n", index_tau, pgb2->tau_sampling_selection[bin][index_tau]);
        //printf("pgb2->tau_sampling_selection_hires[0][%d][%d] = %g\n", index_tau, 0, pgb2->tau_sampling_selection_hires[bin][index_tau][0]);
      }


      /*selection is ruined at this point*/


    printf("First order sources_integ allocated\n" );

    /*  Initialise pgb2->Dl[index_type_first][index_type_second][index_l][bin_first][bin_second][index_tau_first][index_tau_second];
      this is the array that holds the Dirac-angular power spectra */

    class_alloc(pgb2->Dl, pgb2->type_size * sizeof(double ******), pgb2->error_message);
    for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
      class_alloc(pgb2->Dl[index_type_first], pgb2->type_size * sizeof(double *****), pgb2->error_message);
      for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
        class_alloc(pgb2->Dl[index_type_first][index_type_second], ptr->l_size[ppt->index_md_scalars] * sizeof(double ****), pgb2->error_message);
        for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
          class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l], ppt->selection_num * sizeof(double ***), ppt->error_message);
          for (int bin_first = 0; bin_first < ppt->selection_num; bin_first++) {
            class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l][bin_first], ppt->selection_num * sizeof(double**), ppt->error_message);
            for  (int bin_second = 0; bin_second < ppt->selection_num; bin_second++) {
              class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l][bin_first][bin_second], pgb2->tau_size_selection * sizeof(double *), pgb2->error_message);
              for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++) {
                class_alloc(pgb2->Dl[index_type_first][index_type_second][index_l][bin_first][bin_second][index_tau_first], pgb2->tau_size_selection * sizeof(double), pgb2->error_message);
              }
            }
          }
        }
      }
    }





      /* k_min is wrong here */





      for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost; index_tau_bessel++) {
        pgb2->tau_sampling_selection_hires[bin][pgb2->tau_size_selection-1][index_tau_bessel] = pgb2->tau_sampling_selection[bin][pgb2->tau_size_selection-1];
        printf("pgb2->tau_sampling_selection_hires[%d][%d][%d] = %g\n", bin, pgb2->tau_size_selection-1, index_tau_bessel, pgb2->tau_sampling_selection_hires[bin][pgb2->tau_size_selection-1][index_tau_bessel]);
      }
    }



    for ( bin = 0; bin < ppt->selection_num; bin++) {
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
        class_call(array_trapezoidal_weights(pgb2->tau_sampling_selection_hires[bin][index_tau],
                                             bessel_boost*((pgb2->tau_size_selection-1)-index_tau+1),
                                             w_trapz_selection_hires[bin][index_tau],
                                             pgb2->error_message),
                                             pgb2->error_message,
                                             pgb2->error_message);

      }
    }

    /*for ( bin = 0; bin < ppt->selection_num; bin++) {
      for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
        for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*((pgb2->tau_size_selection-1)-index_tau+1); index_tau_bessel++) {
          printf("w_trapz_selection_hires[%d][%d][%d] = %g\n", bin , index_tau, index_tau_bessel, w_trapz_selection_hires[bin][index_tau][index_tau_bessel] );
        }
      }
    }*/







      printf("Allocating size %ix%ix%ix%ix%i bytes \n", pgb2->type_size, pgb2->type_size, ptr->l_size[ppt->index_md_scalars], pgb2->tau_size_cls, pgb2->tau_size_cls);


      int index1=0;
      int index2=0;
      double integ, integ_dens_rsd;
      double result2;


      printf("integrating k between %g and %g\n",ppt->k[ppt->index_md_scalars][0],ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1] );

      double p1,p2,f1,f2,j1,j2;
      int k_test_index;
      double pvecback1;
      double pvecback2;


      double type1, type2;
      double Pk;

      /* Allocate array for Cl3[index_type_first][index_type_second][index_l][bin_first][bin_second] */
      class_alloc(pgb2->Cl3, pgb2->type_size * sizeof(double ****), pgb2->error_message);
      for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
        class_alloc(pgb2->Cl3[index_type_first], pgb2->type_size * sizeof(double ***), pgb2->error_message);
        for (int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
          class_alloc(pgb2->Cl3[index_type_first][index_type_second], ptr->l_size[ppt->index_md_scalars] * sizeof(double **), ppt->error_message);
          for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
            class_alloc(pgb2->Cl3[index_type_first][index_type_second][index_l], ppt->selection_num * sizeof(double *), ppt->error_message);
            for (int bin_first = 0; bin_first < ppt->selection_num; bin_first++) {
              class_alloc(pgb2->Cl3[index_type_first][index_type_second][index_l][bin_first], ppt->selection_num * sizeof(double), ppt->error_message);
              }
            }
          }
        }


  /*==================================================================
  ====================================================================
  =============================   SOURCES  ===========================
  ====================================================================
  ===================================================================*/

  /* In this section we prepare the transfer functions (e.g. the density and velocity transfer functions) and write them into
  first_order_sources, we then use these Fourier modes and pair them with prefactors (e.g. k/a/H) and an associated spherical Bessel
  function (e.g. j''_l) to write angular transfer functions (\Delta_l(k,tau)) into the first_order_sources_integ array which has an
  additional dependence on multipole. The angular transfer functions are then multiplied and integrated over in the next sections
  to calculate power- and bi-spectra.
   */
  printf("Starting sources..\n" );

  double tau_window;
  int last_index_window;
  last_index_window = 0;
  double normal;
  double z;
  double width;
  double mean;
  double x;
  width = ppt->selection_width[0];
  mean = ppt->selection_mean[0];
  double test_sum;
  test_sum = 0.0;


  class_call(background_tau_of_z(
                          pba,
                          ppt->selection_mean[0],
                          &selection_mean),
                          ppt->error_message,
                          pgb2->error_message);

  double t1,t2,k1;
  int index_k1, index_t1, index_t2;
  int index_k_bessel;
  printf("pgb2->k_bessel has %d points that span (%g,%g)\n", pgb2->k_size_bessel, pgb2->k_bessel[0], pgb2->k_bessel[pgb2->k_size_bessel-1]);
  printf("ppt->tau_sampling_quadsources has %d points sampled between (%g,%g)\n", ppt->tau_size_quadsources, ppt->tau_sampling_quadsources[0], ppt->tau_sampling_quadsources[ppt->tau_size_quadsources-1] );
  printf("pgb2->tau_sampling_cls has %d points sampled between (%g,%g)\n", pgb2->tau_size_cls, pgb2->tau_sampling_cls[0], pgb2->tau_sampling_cls[pgb2->tau_size_cls-1] );
  printf("pgb2->tau_sampling_selection has %d points sampled between (%g,%g)\n", pgb2->tau_size_selection, pgb2->tau_sampling_selection[0], pgb2->tau_sampling_selection[pgb2->tau_size_selection-1]);
  printf("pgb2->k_bessel has %d points sampled between (%g, %g)\n", pgb2->k_size_bessel, pgb2->k_bessel[0], pgb2->k_bessel[pgb2->k_size_bessel-1]);




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
      printf("entered source v\n" );
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

          pgb2->first_order_sources[pgb2->index_source_v][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);

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

        }
      }
    }

    intermediate = 0;
    intermediate_plus = 0;
    int last_index_lens;
    int last_index_k_lens;
    int index_k_lens;
    double ** phi_plus_psi;
    double ** phi;
    double ** phi_prime;
    double ** psi;
    double ** psi_prime;
    double ** phi_plus_psi_prime;
    /* Allocate a temporary array phi_plus_psi[index_tau][index_k] that stores the quadsources transfer function for phi+psi
        purely for brevity*/
    class_alloc(phi, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
    class_alloc(phi_prime, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
    class_alloc(psi, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
    class_alloc(psi_prime, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
    class_alloc(phi_plus_psi, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);
    class_alloc(phi_plus_psi_prime, ppt->tau_size_quadsources * sizeof(double *), ppt->error_message);

      for (int index = 0; index < ppt->tau_size_quadsources; index++) {

        class_alloc(phi[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);

        class_alloc(phi_prime[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);

        class_alloc(psi[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);

        class_alloc(psi_prime[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);

        class_alloc(phi_plus_psi[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);

        class_alloc(phi_plus_psi_prime[index],
                    ppt->k_size[ppt->index_md_scalars] * sizeof(double ),
                    ppt->error_message);
      }


    /* Here we are simply renaming the quadsources arrays (from perturbations) for brevity, they are still defined over their native indices and ranges. They will need to
      be interpolated later.*/
    for (int index = 0; index < ppt->tau_size_quadsources; index++) {
      for (int index_k = 0; index_k < ppt->k_size[ppt->index_md_scalars]; index_k++) {
        phi_plus_psi[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_phi][index * ppt->k_size[ppt->index_md_scalars] + index_k]
          + ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_psi][index * ppt->k_size[ppt->index_md_scalars] + index_k];

        phi_plus_psi_prime[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_phi_prime][index * ppt->k_size[ppt->index_md_scalars] + index_k]
          + ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_psi_prime][index * ppt->k_size[ppt->index_md_scalars] + index_k];

        phi[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_phi][index * ppt->k_size[ppt->index_md_scalars] + index_k];

        phi_prime[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_phi_prime][index * ppt->k_size[ppt->index_md_scalars] + index_k];

        psi[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_psi][index * ppt->k_size[ppt->index_md_scalars] + index_k];

        psi_prime[index][index_k] = ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_psi_prime][index * ppt->k_size[ppt->index_md_scalars] + index_k];

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
          index = 0;

          double tau = pgb2->tau_sampling_cls[index_tau];


          class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);

          for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {


            double k = pgb2->k_bessel[index_k_bessel];

            class_call(index_of_k_old(k,
                                     &index_k,
                                     &last_index_k,
                                     ppt),
                                     pgb2->error_message,
                                     pgb2->error_message);

            f = (tau-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

            intermediate  = f*psi[index+1][index_k]+(1-f)*psi[index][index_k];

            intermediate_plus =  f*psi[index+1][index_k+1]+(1-f)*psi[index][index_k+1];

            g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

            pgb2->first_order_sources[pgb2->index_source_psi][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);
          }
      }
    }



    if (pgb2->index_source_phi_plus_psi != -1) {
      for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){

          last_index_k = 0;
          index_k_lens = 0;
          int last_index;
          int index_tau_lens;
          double lensing_result;

          index = 0;

          double tau = pgb2->tau_sampling_cls[index_tau];

          class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);

          for (int index_k_bessel = 0; index_k_bessel < bessel_boost*(index_tau+1)/*pgb2->k_size_bessel*/; index_k_bessel++) {



            double k = pgb2->k_bessel[index_k_bessel];

            class_call(index_of_k_old(k,
                           &index_k,
                           &last_index_k,
                           ppt),
                           pgb2->error_message,
                           pgb2->error_message);

            f = (tau-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

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
          int last_index;
          int index_tau_lens;
          double lensing_result;
          f = 0.;
          g = 0.;

          for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

            index = 0;

            double tau = pgb2->tau_sampling_cls[index_tau];

            class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);

            double k = pgb2->k_bessel[index_k_bessel];

            class_call(index_of_k_old(k,
                           &index_k,
                           &last_index_k,
                           ppt),
                           pgb2->error_message,
                           pgb2->error_message);

            f = (tau-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

            intermediate  = f*phi[index+1][index_k]+(1-f)*phi[index][index_k];

            intermediate_plus =  f*phi[index+1][index_k+1]+(1-f)*phi[index][index_k+1];

            g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

            pgb2->first_order_sources[pgb2->index_source_phi][index_tau][index_k_bessel] = (g*intermediate_plus +(1-g)*intermediate);
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


        pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi_prime][0][index_k_bessel] = (pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][1][index_k_bessel]-pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][0][index_k_bessel])/
          (pgb2->tau_sampling_cls[1]-pgb2->tau_sampling_cls[0]);

        pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi_prime][pgb2->tau_size_cls-1][index_k_bessel] = (pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][pgb2->tau_size_cls-1][index_k_bessel]-pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][pgb2->tau_size_cls-2][index_k_bessel])/
          (pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]-pgb2->tau_sampling_cls[pgb2->tau_size_cls-2]);
      }


      for (int index_tau = 1; index_tau < pgb2->tau_size_cls-1; index_tau++){
        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
          pgb2->first_order_sources[pgb2->index_source_phi_prime][index_tau][index_k_bessel] = (pgb2->first_order_sources[pgb2->index_source_phi][index_tau+1][index_k_bessel]-pgb2->first_order_sources[pgb2->index_source_phi][index_tau-1][index_k_bessel])
                                                                                               /(pgb2->tau_sampling_cls[index_tau+1]-pgb2->tau_sampling_cls[index_tau-1]);

          pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi_prime][index_tau][index_k_bessel] = (pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][index_tau+1][index_k_bessel]-pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][index_tau-1][index_k_bessel])
                                                                                                                  /(pgb2->tau_sampling_cls[index_tau+1]-pgb2->tau_sampling_cls[index_tau-1]);
        }
      }
    }



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


      printf("Entering angular transfer function preparation.\n");
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

                pgb2->first_order_sources_integ[pgb2->index_type_density][index_l][index_tau][index_k_bessel] = g_bias*pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_bessel] * j;
              }
            }
          }
        }

        /* Density prime */
        if(pgb2->index_type_quad_density_p != -1){
          double * pvecback_quad_density_p;
          class_alloc(pvecback_quad_density_p, pba->bg_size * sizeof(double), pba->error_message);
          int last_index_quad_density_p;
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_quad_density_p,
                                           pvecback_quad_density_p),
                                           pba->error_message,
                                           pgb2->error_message);

              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);
                double bias = 1.3;
                double prefactor_quad_density_p = pgb2->k_bessel[index_k_bessel]
                                                 /pvecback_quad_density_p[pba->index_bg_H]
                                                 /pvecback_quad_density_p[pba->index_bg_a];

                class_call(bessel_at_x_first_deriv(pgb2, pbs, x, index_l, &j_first_deriv), pbs->error_message, pgb2->error_message);

                pgb2->first_order_sources_integ[pgb2->index_type_quad_density_p][index_l][index_tau][index_k_bessel] = bias
                                                                                                                      *prefactor_quad_density_p
                                                                                                                      *pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_bessel]
                                                                                                                      *j_first_deriv;

              }
            }
          }
        }

        printf("pgb2->index_type_quad_v = %d\n", pgb2->index_type_quad_v);
        if(pgb2->index_type_quad_v != -1){
          double * pvecback_quad_v;
          class_alloc(pvecback_quad_v, pba->bg_size * sizeof(double), pba->error_message);
          int last_index_quad_v;
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_quad_v,
                                           pvecback_quad_v),
                                           pba->error_message,
                                           pgb2->error_message);

              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);
                double bias = 1.0;
                double prefactor_quad_v = 1.0
                                         /pgb2->k_bessel[index_k_bessel];

                class_call(bessel_at_x_first_deriv(pgb2, pbs, x, index_l, &j_first_deriv), pbs->error_message, pgb2->error_message);

                pgb2->first_order_sources_integ[pgb2->index_type_quad_v][index_l][index_tau][index_k_bessel] = prefactor_quad_v
                                                                                                               *pgb2->first_order_sources[pgb2->index_source_theta][index_tau][index_k_bessel]
                                                                                                               *j_first_deriv;
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
                double k1 = pgb2->k_bessel[index_k_bessel];
                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);

                class_call(bessel_at_x_second_deriv(pgb2, pbs, x, index_l, &j_second_deriv), pbs->error_message, pgb2->error_message);

                double prefactor_rsd = 1.0
                             /pvecback_rsd[pba->index_bg_H]
                             /pvecback_rsd[pba->index_bg_a];

              /* Write in to the array */
                pgb2->first_order_sources_integ[pgb2->index_type_rsd][index_l][index_tau][index_k_bessel] = prefactor_rsd
                                *pgb2->first_order_sources[pgb2->index_source_theta][index_tau][index_k_bessel]*
                                j_second_deriv;



              }
            }
          }
        }

        printf("completed RSD source integ.\n");

        /* This is for a quad-source; a term that is linear but appears specifically at second-order multiplied with one other quadsource.
        In the full galaxy number-over density at second order (\Delta^{(2)}) the term is  (-2/a/a/H/Hv')*(v'''). We separate the two
        brackets and treat them as quad_v_p (left) and quad_v_ppp (right), with the "p" denoting the derivative. */

        /* Eq. 3.32 in [1812.09297] */

        double * pvecback_quad_v_p;
        class_alloc(pvecback_quad_v_p, pba->bg_size * sizeof(double), pba->error_message);

        if(pgb2->index_type_quad_v_p != -1 ){
          int index_tau_quad;
          int last_index_quad_v_p;
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){

              index_tau_quad = 0;

              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_quad_v_p,
                                           pvecback_quad_v_p),
                                           pba->error_message,
                                           pgb2->error_message);

              double tau = pgb2->tau_sampling_cls[index_tau];

              class_call(index_of_tau_sampling_quadsources(tau, &index_tau_quad, ppt), pgb2->error_message, pgb2->error_message);
              index_k = 0;
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

                double k = pgb2->k_bessel[index_k_bessel];

                class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);

                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);

                class_call(bessel_at_x_second_deriv(pgb2, pbs, x, index_l, &j_first_deriv), pbs->error_message, pgb2->error_message);


                prefactor1 = 1.0
                             /pvecback_quad_v_p[pba->index_bg_H]
                             /pvecback_quad_v_p[pba->index_bg_a];



                pgb2->first_order_sources_integ[pgb2->index_type_quad_v_p][index_l][index_tau][index_k_bessel] = prefactor1*pgb2->first_order_sources[pgb2->index_source_theta][index_tau][index_k_bessel]*j_second_deriv;
              }
            }
          }
        }



        if(pgb2->index_type_quad_v_pp != -1 ){
          double * pvecback_quad_v_pp;
          class_alloc(pvecback_quad_v_pp, pba->bg_size * sizeof(double), pba->error_message);
          double j_third_deriv;
          int index_tau_quad;
          int last_index_quad_v_pp;
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){

              index_tau_quad = 0;

              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_quad_v_pp,
                                           pvecback_quad_v_pp),
                                           pba->error_message,
                                           pgb2->error_message);

              double tau = pgb2->tau_sampling_cls[index_tau];

              class_call(index_of_tau_sampling_quadsources(tau, &index_tau_quad, ppt), pgb2->error_message, pgb2->error_message);
              index_k = 0;
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

                double k = pgb2->k_bessel[index_k_bessel];

                class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);

                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);

                class_call(bessel_at_x_third_deriv(pgb2, pbs, x, index_l, &j_third_deriv), pbs->error_message, pgb2->error_message);


               double prefactor_quad_v_pp = k
                                            /pvecback_quad_v_pp[pba->index_bg_H]
                                            /pvecback_quad_v_pp[pba->index_bg_a]
                                            /pvecback_quad_v_pp[pba->index_bg_H]
                                            /pvecback_quad_v_pp[pba->index_bg_a];



                pgb2->first_order_sources_integ[pgb2->index_type_quad_v_pp][index_l][index_tau][index_k_bessel] = prefactor_quad_v_pp
                                                                                                                  *pgb2->first_order_sources[pgb2->index_source_theta][index_tau][index_k_bessel]
                                                                                                                  *j_third_deriv;

              }
            }
          }
        }





        /* First Type: Doppler1 */
        if(pgb2->index_type_d1 != -1 ){
          int last_index_d1;
          double * pvecback_d1;
          class_alloc(pvecback_d1, pba->bg_size * sizeof(double), pba->error_message);
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
                //warning f_evo set to zero
                f_evo1 = 0;
                         /*2.
                         /pvecback_d1[pba->index_bg_H]
                         /pvecback_d1[pba->index_bg_a]
                         /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau])
                         +pvecback_d1[pba->index_bg_H_prime]
                         /pvecback_d1[pba->index_bg_H]
                         /pvecback_d1[pba->index_bg_H]
                         /pvecback_d1[pba->index_bg_a];*/

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

        /* First Type: Doppler2 */
        int last_index_d2;
        if(pgb2->index_type_d2 != -1){
          double * pvecback_d2;
          class_alloc(pvecback_d2, pba->bg_size * sizeof(double), pba->error_message);
          double j_d2;
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
                class_call(bessel_at_x(pbs, x , index_l, &j_d2), pbs->error_message, pgb2->error_message);

                f_evo2 = 2.
                         /pvecback_d2[pba->index_bg_H]
                         /pvecback_d2[pba->index_bg_a]
                         /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau])
                         +pvecback_d2[pba->index_bg_H_prime]
                         /pvecback_d2[pba->index_bg_H]
                         /pvecback_d2[pba->index_bg_H]
                         /pvecback_d2[pba->index_bg_a];
                         //alert f_evo1 skipped
                //prefactor1 = -3.0*pvecback_d2[pba->index_bg_a]*pvecback_d2[pba->index_bg_H];
                //WARNING f_evo set to zero
                pgb2->first_order_sources_integ[pgb2->index_type_d2][index_l][index_tau][index_k_bessel] = (0.0/*f_evo2*/-3.0)
                        *pvecback_d2[pba->index_bg_a]
                        *pvecback_d2[pba->index_bg_H]
                        *pgb2->first_order_sources[pgb2->index_source_theta][index_tau][index_k_bessel]
                        *j_d2
                        /pgb2->k_bessel[index_k_bessel]
                        /pgb2->k_bessel[index_k_bessel];
              }
            }
          }
        }

        /* First Type: g1 (first of the GR terms) */
        if(pgb2->index_type_g1 != -1){
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            index = 0;
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){


              double tau = pgb2->tau_sampling_cls[index_tau];

              class_call(index_of_tau_sampling_quadsources(tau, &index, ppt), pgb2->error_message, pgb2->error_message);
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {


                double k = pgb2->k_bessel[index_k_bessel];

                //class_call(index_of_k(k, &index_k_lens, ppt), pgb2->error_message, pgb2->error_message);
                class_call(index_of_k_old(k,
                               &index_k,
                               &last_index_k,
                               ppt),
                               pgb2->error_message,
                               pgb2->error_message);

                f = (tau-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

                intermediate  = f*phi[index+1][index_k]+(1-f)*phi[index][index_k];

                intermediate_plus =  f*phi[index+1][index_k+1]+(1-f)*phi[index][index_k+1];

                g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

                double phi = (g*intermediate_plus +(1-g)*intermediate);

                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);

                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);


                /* Warning not the same as in CLASSgal paper 1307.1459 A.20 which has -2+5s prefactor !? */
                double prefactor_g1 = 1.;//-2.0+5.0*ptr->s_bias;

                pgb2->first_order_sources_integ[pgb2->index_type_g1][index_l][index_tau][index_k_bessel] = prefactor_g1
                        /*pgb2->first_order_sources[pgb2->index_source_phi][index_tau][index_k_bessel]*/
                        *phi
                        *j;
              }
            }
          }
        }

        /* First Type: g2 */

        if(pgb2->index_type_g2 != -1 ){
          double * pvecback_g2;
          class_alloc(pvecback_g2, pba->bg_size * sizeof(double), pba->error_message);
          int last_index_g2;
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){

              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_g2,
                                           pvecback_g2),
                                           pba->error_message,
                                           pgb2->error_message);

              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);
                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);

                f_evo1 = 0.;
                         /*2.
                         /pvecback_g2[pba->index_bg_H]
                         /pvecback_g2[pba->index_bg_a]
                         /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau])
                         +pvecback_g2[pba->index_bg_H_prime]
                         /pvecback_g2[pba->index_bg_H]
                         /pvecback_g2[pba->index_bg_H]
                         /pvecback_g2[pba->index_bg_a];*/

                double prefactor_g2 = -(3.0
                             +pvecback_g2[pba->index_bg_H_prime]
                             /pvecback_g2[pba->index_bg_H]
                             /pvecback_g2[pba->index_bg_H]
                             /pvecback_g2[pba->index_bg_a]
                             +(2.0-5.0*ptr->s_bias)
                             /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau])
                             /pvecback_g2[pba->index_bg_H]
                             /pvecback_g2[pba->index_bg_a]
                             -f_evo1);


                pgb2->first_order_sources_integ[pgb2->index_type_g2][index_l][index_tau][index_k_bessel] = prefactor_g2
                        *pgb2->first_order_sources[pgb2->index_source_psi][index_tau][index_k_bessel]
                        *j;
              }
            }
          }
        }

        /* First Type: g3 */

        if(pgb2->index_type_g3 != -1 ){
          double * pvecback_g3;
          class_alloc(pvecback_g3, pba->bg_size * sizeof(double), pba->error_message);
          int index_tau_quad;
          int last_index_g3;
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              index_tau_quad = 0;

              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_g3,
                                           pvecback_g3),
                                           pba->error_message,
                                           pgb2->error_message);

              double tau = pgb2->tau_sampling_cls[index_tau];


              class_call(index_of_tau_sampling_quadsources(tau, &index_tau_quad, ppt), pgb2->error_message, pgb2->error_message);
              index_k = 0;
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {





                double k = pgb2->k_bessel[index_k_bessel];

                class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);
                /*class_call(index_of_k_old(k,
                               &index_k,
                               &last_index_k,
                               ppt),
                               pgb2->error_message,
                               pgb2->error_message);*/

                f = (tau-ppt->tau_sampling_quadsources[index_tau_quad])/(ppt->tau_sampling_quadsources[index_tau_quad+1]-ppt->tau_sampling_quadsources[index_tau_quad]);

                intermediate  = f*phi_prime[index_tau_quad+1][index_k]+(1-f)*phi_prime[index_tau_quad][index_k];

                intermediate_plus =  f*phi_prime[index_tau_quad+1][index_k+1]+(1-f)*phi_prime[index_tau_quad][index_k+1];

                g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

                double phi_prime = (g*intermediate_plus +(1-g)*intermediate);

                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);
                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);


                prefactor1 = 1.0
                             /pvecback_g3[pba->index_bg_H]
                             /pvecback_g3[pba->index_bg_a];

                pgb2->first_order_sources_integ[pgb2->index_type_g3][index_l][index_tau][index_k_bessel] = prefactor1
                        *phi_prime
                        *j;
              }
            }
          }
        }



        double * pvecback_delta;

        class_alloc(pvecback_delta, pba->bg_size * sizeof(double), pba->error_message);
        if(pgb2->index_type_delta != -1 ){
          int index_tau_quad;
          int last_index_delta;
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              index_tau_quad = 0;

              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_delta,
                                           pvecback_delta),
                                           pba->error_message,
                                           pgb2->error_message);

              double tau = pgb2->tau_sampling_cls[index_tau];

              class_call(index_of_tau_sampling_quadsources(tau, &index_tau_quad, ppt), pgb2->error_message, pgb2->error_message);
              index_k = 0;
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

                double k = pgb2->k_bessel[index_k_bessel];

                class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);

                double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]);

                class_call(bessel_at_x_second_deriv(pgb2, pbs, x, index_l, &j_second_deriv), pbs->error_message, pgb2->error_message);
                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);

                prefactor1 = 1.
                             /pvecback_delta[pba->index_bg_H]
                             /pvecback_delta[pba->index_bg_a];



                double term1 = pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_bessel]*j;

                double term2 = pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_bessel]*j_second_deriv;

                pgb2->first_order_sources_integ[pgb2->index_type_delta][index_l][index_tau][index_k_bessel] =
                      pgb2->first_order_sources_integ[pgb2->index_type_rsd][index_l][index_tau][index_k_bessel]
                      + pgb2->first_order_sources_integ[pgb2->index_type_density][index_l][index_tau][index_k_bessel];

                /*double Omega_m0 = pba->Omega0_cdm + pba->Omega0_b;
                /* infer redshift */
                /*double z = pba->a_today/pvecback_delta[pba->index_bg_a]-1.;
                double Omega_l0 = 1 - Omega_m0;
                double Ez = sqrt(Omega_m0*pow(1+z,3) + Omega_l0);
                double Omega_m = Omega_m0*pow(1+z,3)/(Ez*Ez);
                double f = pow(Omega_m,4/7.);
                //herehere
                double paper_delta = pgb2->first_order_sources_integ[pgb2->index_type_density][index_l][index_tau][index_k_bessel]-f*term2;

                printf("%g      %g\n", pgb2->first_order_sources_integ[pgb2->index_type_delta][index_l][index_tau][index_k_bessel], paper_delta);*/

              }
            }
          }
        }







    /* Lensing source term preparation */

    double j_lens;
    int index_k_lens2;


    if (pgb2->index_type_lens != -1) {
      printf("starting pgb2->index_type_lens\n" );

      int index_tau_quad;
      int last_index_k;
      int index_k;
      //for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      for (int index_l = 2; index_l < 3; index_l++) {
      //for (int index_l = ptr->l_size[ppt->index_md_scalars]-1; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      //for (int index_l = ptr->l_size[ppt->index_md_scalars]-2; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
        index_k = 0;
        index_k_lens2 = 0;
        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

          double k = pgb2->k_bessel[index_k_bessel];
          printf("index_k: %d, k = %g\n",index_k_bessel, k);
          class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);

          for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
            double lens_sum = 0.;
            printf("inside lens loop: index_l x index_k_bessel x index_tau = %dx%dx%d\n", index_l, index_k_bessel, index_tau);
            /* Avoid the last index (tau_sampling_bessel2[index_tau][bessel_boost*(index_tau+1)] = nan) to avoid division by zero */
            for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(index_tau+1)-1; index_tau_bessel++) {


              index_tau_quad = 0;

              double tau = pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel];


              class_call(index_of_tau_sampling_quadsources(tau, &index_tau_quad, ppt), pgb2->error_message, pgb2->error_message);

              double k = pgb2->k_bessel[index_k_bessel];

              class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);
              /*class_call(index_of_k_old(k,
                             &index_k,
                             &last_index_k,
                             ppt),
                             pgb2->error_message,
                             pgb2->error_message);*/

              f = (tau-ppt->tau_sampling_quadsources[index_tau_quad])/(ppt->tau_sampling_quadsources[index_tau_quad+1]-ppt->tau_sampling_quadsources[index_tau_quad]);

              intermediate  = f*phi_plus_psi[index_tau_quad+1][index_k]+(1-f)*phi_plus_psi[index_tau_quad][index_k];

              intermediate_plus =  f*phi_plus_psi[index_tau_quad+1][index_k+1]+(1-f)*phi_plus_psi[index_tau_quad][index_k+1];

              g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

              double phi_plus_psi = (g*intermediate_plus +(1-g)*intermediate);
              double chi_cls = tau0-pgb2->tau_sampling_cls[index_tau];
              double chi_lens = tau0-pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel];

              double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age-pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel]);

              class_call(bessel_at_x(pbs,x, index_l, &j), pbs->error_message, pgb2->error_message);


              double chi_fraction = (chi_cls-chi_lens)/chi_cls/chi_lens;

              lens_sum += (1.)*ptr->l[index_l]*(ptr->l[index_l]+1)*(2.-5*ptr->s_bias)*chi_fraction*phi_plus_psi
                          *j*w_trapz_lens_bessel2[index_tau][index_tau_bessel]/2.;


            }
            pgb2->first_order_sources_integ[pgb2->index_type_lens][index_l][index_tau][index_k_bessel] = lens_sum;
            printf("lens_sum = %g\n", lens_sum );

          }
        }
      }
    }

    if (pgb2->index_type_g4 != -1) {


      int index_tau_quad;
      int last_index_k;
      int index_k;
      //for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      for (int index_l = 2; index_l < 3; index_l++) {
      //for (int index_l = ptr->l_size[ppt->index_md_scalars]-1; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      //for (int index_l = ptr->l_size[ppt->index_md_scalars]-2; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
      //for (int index_l = 0; index_l < 1; index_l++) {
        index_k = 0;
        index_k_lens2 = 0;
        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

          double k = pgb2->k_bessel[index_k_bessel];
          printf("index_k: %d, k = %g\n",index_k_bessel, k);
          class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);

          for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
            double g4_sum = 0.;
            printf("inside g4 loop: index_l x index_k_bessel x index_tau = %dx%dx%d\n", index_l, index_k_bessel, index_tau);

            /* Avoid the last index (tau_sampling_bessel2[index_tau][bessel_boost*(index_tau+1)] = nan) to avoid division by zero */
            for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(index_tau+1)-1; index_tau_bessel++) {

              index_tau_quad = 0;

              double tau = pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel];

              class_call(index_of_tau_sampling_quadsources(tau, &index_tau_quad, ppt), pgb2->error_message, pgb2->error_message);

              double k = pgb2->k_bessel[index_k_bessel];

              class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);

              /* Interpolate */
              f = (tau-ppt->tau_sampling_quadsources[index_tau_quad])/(ppt->tau_sampling_quadsources[index_tau_quad+1]-ppt->tau_sampling_quadsources[index_tau_quad]);

              intermediate  = f*phi_plus_psi[index_tau_quad+1][index_k]+(1-f)*phi_plus_psi[index_tau_quad][index_k];

              intermediate_plus =  f*phi_plus_psi[index_tau_quad+1][index_k+1]+(1-f)*phi_plus_psi[index_tau_quad][index_k+1];

              g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

              double phi_plus_psi = (g*intermediate_plus +(1-g)*intermediate);

              double chi_lens = tau0-pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel];

              double chi_cls = tau0 - pgb2->tau_sampling_cls[index_tau];

              double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age-pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel]);

              class_call(bessel_at_x(pbs,x, index_l, &j), pbs->error_message, pgb2->error_message);
              // WARNING magnification edited
              g4_sum += (2.-5*ptr->s_bias)
                        *phi_plus_psi
                        *j
                        *w_trapz_lens_bessel2[index_tau][index_tau_bessel]
                        /chi_cls;

            }
            pgb2->first_order_sources_integ[pgb2->index_type_g4][index_l][index_tau][index_k_bessel] = g4_sum;
            printf("g4_sum = %g\n", g4_sum);

          }
        }
      }
    }

    if (pgb2->index_type_g5 != -1) {


      int index_tau_quad;
      int last_index_k;
      int index_tau_cls;
      int index_k;
      double * pvecbackg;
      int last_index_g;
      double g5_prefactor;
      class_alloc(pvecbackg, pba->bg_size*sizeof(double), pba->error_message);

      for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      //for (int index_l = 3; index_l < 4; index_l++) {
      //for (int index_l = ptr->l_size[ppt->index_md_scalars]-2; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
        index_k = 0;
        index_k_lens2 = 0;
        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

          double k = pgb2->k_bessel[index_k_bessel];

          class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);

          for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
            double g5_sum = 0.;
            printf("inside g5 loop: index_l x index_k_bessel x index_tau = %dx%dx%d\n", index_l, index_k_bessel, index_tau);
            index_tau_cls = 0;

            class_call(background_at_tau(pba,
                                         pgb2->tau_sampling_cls[index_tau],
                                         pba->long_info,
                                         pba->inter_normal,
                                         &last_index_g,
                                         pvecbackg),
                                         pba->error_message,
                                         pgb2->error_message);
           // WARNING f_evo hardcoded=0
           double f_evo = 0;
               /*2.
                /pvecbackg[pba->index_bg_H]
                /pvecbackg[pba->index_bg_a]
                /(pba->conformal_age - pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel])
                +pvecbackg[pba->index_bg_H_prime]
                /pvecbackg[pba->index_bg_H]
                /pvecbackg[pba->index_bg_H]
                /pvecbackg[pba->index_bg_a];*/
            g5_prefactor = (1.0
                           +pvecbackg[pba->index_bg_H_prime]
                           /pvecbackg[pba->index_bg_a]
                           /pvecbackg[pba->index_bg_H]
                           /pvecbackg[pba->index_bg_H]
                           +(2.0-5.0*ptr->s_bias)
                           /(pba->conformal_age - pgb2->tau_sampling_cls[index_tau])
                           /pvecbackg[pba->index_bg_a]
                           /pvecbackg[pba->index_bg_H]
                           +5.0*ptr->s_bias
                           -f_evo);
            /* Avoid the last index (tau_sampling_bessel2[index_tau][bessel_boost*(index_tau+1)] = nan) to avoid division by zero */
            for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(index_tau+1)-1; index_tau_bessel++) {

              //printf("%dx%dx%dx%d\n", index_l, index_k_bessel, index_tau, index_tau_bessel);

              index_tau_quad = 0;

              double tau = pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel];


              class_call(index_of_tau_sampling_quadsources(tau, &index_tau_quad, ppt), pgb2->error_message, pgb2->error_message);

              double k = pgb2->k_bessel[index_k_bessel];

              class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);



              /* Interpolate */
              //f = (tau-ppt->tau_sampling_quadsources[index_tau_quad])/(ppt->tau_sampling_quadsources[index_tau_quad+1]-ppt->tau_sampling_quadsources[index_tau_quad]);
            //  printf("index_tau_cls =%d\n", index_tau_cls);
            //  printf("index_tau =%d\n", index_tau);
            //  printf("index_tau_bessel = %d\n", index_tau_bessel);

              /*if (index_tau_cls == pgb2->tau_size_cls-1) {

                phi_plus_psi_prime = pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi_prime][index_tau_cls][index_k_bessel];

                //printf("entered *first* if, phi_plus_psi_prime = %g\n", phi_plus_psi_prime);
              }

              else{
                //printf("entered *second* if\n");
                f = (tau-pgb2->tau_sampling_cls[index_tau_cls])/(pgb2->tau_sampling_cls[index_tau_cls+1]-pgb2->tau_sampling_cls[index_tau_cls]);
                //printf("f = %g\n", f );

                phi_plus_psi_prime = f*pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi_prime][index_tau_cls+1][index_k_bessel]+(1-f)*pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi_prime][index_tau_cls][index_k_bessel];
              }*/

              /*f = (tau-ppt->tau_sampling_quadsources[index_tau_quad])/(ppt->tau_sampling_quadsources[index_tau_quad+1]-ppt->tau_sampling_quadsources[index_tau_quad]);

              intermediate  = f*phi_plus_psi_prime[index_tau_quad+1][index_k]+(1-f)*phi_plus_psi_prime[index_tau_quad][index_k];

              intermediate_plus =  f*phi_plus_psi_prime[index_tau_quad+1][index_k+1]+(1-f)*phi_plus_psi_prime[index_tau_quad][index_k+1];

              g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

              double phi_plus_psi_prime = (g*intermediate_plus +(1-g)*intermediate);*/

              f = (tau-ppt->tau_sampling_quadsources[index_tau_quad])/(ppt->tau_sampling_quadsources[index_tau_quad+1]-ppt->tau_sampling_quadsources[index_tau_quad]);

              intermediate  = f*phi_plus_psi[index_tau_quad+1][index_k]+(1-f)*phi_plus_psi[index_tau_quad][index_k];

              intermediate_plus =  f*phi_plus_psi[index_tau_quad+1][index_k+1]+(1-f)*phi_plus_psi[index_tau_quad][index_k+1];

              g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

              double phi_plus_psi = (g*intermediate_plus +(1-g)*intermediate);


              double x = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age-pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel]);

              //class_call(bessel_at_x(pbs,x, index_l, &j), pbs->error_message, pgb2->error_message);

              class_call(bessel_at_x_first_deriv(pgb2, pbs, x, index_l, &j_first_deriv), pbs->error_message, pgb2->error_message);

              g5_sum += /*g5_prefactor*/phi_plus_psi*k
                          *j_first_deriv*w_trapz_lens_bessel2[index_tau][index_tau_bessel];

              //printf("g5_prefactor = %g\n", g5_prefactor);
              //printf("x = %g\n", x);
              /*if (j != 0.0) {
                printf("j = %g\n", j);
              }*/

              //printf("phi_plus_psi_prime = %g\n", phi_plus_psi_prime);
              //printf("w_trapz_lens_bessel2 = %g\n", w_trapz_lens_bessel2[index_tau][index_tau_bessel]);

            }
            pgb2->first_order_sources_integ[pgb2->index_type_g5][index_l][index_tau][index_k_bessel] = g5_sum;
            printf("g5_sum = %g\n", g5_sum);

          }
        }
      }
    }



  printf("Sources complete.\n");


 /* Second-order source term preparation */
 //ppt2->sources[ppt2->index_tp2_delta_cdm][index_k1][index_k2][index_tau_ppt2*ppt2->k3_size[index_k1][index_k2]+index_k3];




 /*if (pgb2->index_source_delta_cdm_so != -1 && configuration == equilateral) {
   printf("Preparing density source term..\n");
   for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
     index_k = 0;

     for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
       index = 0;

       double tau = pgb2->tau_sampling_cls[index_tau];

       class_call(index_of_tau_sampling_quadsources(tau, &index, ppt),ppt->error_message,pgb2->error_message);
       class_call(index_of_tau_sampling_ppt2(tau,
                      &index_tau_ppt2,
                      &last_tau_kernel,
                      pgb2), pgb2->error_message, pgb2->error_message);

       double k = pgb2->k_bessel[index_k_bessel];
       // WARNING CHECK WHETHER tau_size_quadsources is the corresponsing time grid to ppt2->sources
       class_call(index_of_k(k, &index_k, ppt), ppt->error_message, pgb2->error_message);

       f = (pgb2->tau_sampling_cls[index_tau]-ppt->tau_sampling_quadsources[index])/(ppt->tau_sampling_quadsources[index+1]-ppt->tau_sampling_quadsources[index]);

       intermediate  = (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k]+
                        f*ppt2->sources[ppt2->index_tp2_delta_cdm][index_k1][index_k2][index_tau_ppt2*ppt2->k3_size[index_k1][index_k2]+index_k3]
           (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k]);

       intermediate_plus =  (f*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][(index+1) * ppt->k_size[ppt->index_md_scalars] + index_k+1]+
           (1-f)*ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index * ppt->k_size[ppt->index_md_scalars] + index_k+1]);

       g = (pgb2->k_bessel[index_k_bessel]-ppt->k[ppt->index_md_scalars][index_k])/(ppt->k[ppt->index_md_scalars][index_k+1]-ppt->k[ppt->index_md_scalars][index_k]);

       pgb2->second_order_sources_eq[pgb2->index_source_delta_cdm_so][index_tau][index_k_bessel] = g*intermediate_plus +(1-g)*intermediate;


     }
   }*/



  /*==================================================================
  ====================================================================
  =======================   INTEGRATIONS ============================
  ====================================================================
  ===================================================================*/


  printf("Starting integrations..\n");

  /* Allocate array for class_xfer[index_type][bin][index_l][index_k] */
  double **** class_xfer;
  class_alloc(class_xfer, pgb2->type_size * sizeof(double ****), pgb2->error_message);

  for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){

    class_alloc(class_xfer[index_type_first], ppt->selection_num * sizeof(double ***), pgb2->error_message);

    for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {

      class_alloc(class_xfer[index_type_first][bin1], ptr->l_size[ppt->index_md_scalars] * sizeof(double **), pgb2->error_message);

      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){

        class_alloc(class_xfer[index_type_first][bin1][index_l], pgb2->k_size_bessel * sizeof(double *), pgb2->error_message);

      }
    }
  }
  FILE * Dl_file;

  // use appropriate location if you are using MacOS or Linux
  Dl_file = fopen("output/Dl_file.dat","w");


  printf("Starting SONG style fixed grid integrations..\n");

  double source_interp1;
  double source_interp2;

  double tau_first;
  double tau_second;
  int index_of_cls1;
  int index_of_cls2;
  double Pk_song;
  for(int index_type_first = 1; index_type_first < pgb2->type_size; index_type_first++){
    // Alert only looping over the first index in type_second
    for(int index_type_second = 0; index_type_second < /*index_type_first+1 */1/*pgb2->type_size*/; index_type_second++){
      fprintf(Dl_file,"### index_type_first, index_type_second = %d, %d \n", index_type_first, index_type_second);
    //Warning
    //for(int index_type_first = 0; index_type_first < 1; index_type_first++){
      //for(int index_type_second = 1; index_type_second < 2; index_type_second++){
      for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
        for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
          //for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
          for (int index_l = 2; index_l < 3; index_l++) {
          //for(int index_l = ptr->l_size[ppt->index_md_scalars]-1; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
          //for(int index_l = ptr->l_size[ppt->index_md_scalars]-1; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
          //for (int index_l = ptr->l_size[ppt->index_md_scalars]-2; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
          //for (int index_l = 3; index_l < 4; index_l++) {
            printf("index_l = %d (l=%d)\n", index_l, ptr->l[index_l]);
            double tau_sum2 = 0.0;
            for (int index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++) {
              tau_second = pgb2->tau_sampling_selection[bin2][index_tau_second];
              index_of_tau_sampling_cls(tau_second, &index_of_cls2, pgb2);
              double tau_sum1 = 0.0;

              for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++) {

                tau_first = pgb2->tau_sampling_selection[bin1][index_tau_first];

                index_of_tau_sampling_cls(tau_first, &index_of_cls1, pgb2);
                double k_sum1 = 0.0;

                for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

                  /* We loop over the tau-selection indices and interpolate on the tau_cls grid */

                    class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, pgb2->k_bessel[index_k_bessel], &Pk_song), ppm->error_message, pgb2->error_message);


                    source_interp1 = pgb2->first_order_sources_integ[index_type_first][index_l][index_of_cls1-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_of_cls1]-tau_first)
                                  + pgb2->first_order_sources_integ[index_type_first][index_l][index_of_cls1][index_k_bessel]*(tau_first-pgb2->tau_sampling_cls[index_of_cls1-1]);
                    source_interp1 /= (pgb2->tau_sampling_cls[index_of_cls1] - pgb2->tau_sampling_cls[index_of_cls1-1]);

                    source_interp2 = pgb2->first_order_sources_integ[index_type_second][index_l][index_of_cls2-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_of_cls2]-tau_second)
                                  + pgb2->first_order_sources_integ[index_type_second][index_l][index_of_cls2][index_k_bessel]*(tau_second-pgb2->tau_sampling_cls[index_of_cls2-1]);
                    source_interp2 /= (pgb2->tau_sampling_cls[index_of_cls2] - pgb2->tau_sampling_cls[index_of_cls2-1]);


                    k_sum1 += 4. * _PI_ *  Pk_song * pow(pgb2->k_bessel[index_k_bessel],-1.0)* source_interp1*source_interp2*pgb2->w_trapz_k[index_k_bessel];

                    }

                    pgb2->Dl[index_type_first][index_type_second][index_l][bin1][bin2][index_tau_first][index_tau_second] = k_sum1;
                    if (index_l == ptr->l_size[ppt->index_md_scalars]-1) {
                      fprintf(Dl_file,"%g       %g        %g \n", tau_first, tau_second,  k_sum1);
                    }


                    //printf("pgb2->Dl[%d][%d][%d][%d][%d][%d][%d] = %g\n", index_type_first, index_type_second, index_l , bin1, bin2, index_tau_first, index_tau_second, k_sum1);

                  tau_sum1 += k_sum1*selection[bin1][index_tau_first]*w_trapz[bin1][index_tau_first];
                  }
                tau_sum2 += tau_sum1*selection[bin2][index_tau_second]*w_trapz[bin2][index_tau_second];
            }

            pgb2->Cl3[index_type_first][index_type_second][index_l][bin1][bin2] = tau_sum2;
            printf("Linear fixed SONG style: pgb2->Cl3[%d][%d][%d][%d][%d] = %g\n",
                    index_type_first,
                    index_type_second,
                    index_l,
                    bin1,
                    bin2,
                    ptr->l[index_l]*(ptr->l[index_l]+1)*pgb2->Cl3[index_type_first][index_type_second][index_l][bin1][bin2]/(2*_PI_));
            }
          }
        }
      }
    }
    fclose(Dl_file);
    printf("###l x Dl\n" );
    printf("#selection size = %d\n", pgb2->tau_size_selection);
    printf("#selection min/max %g/%g\n", pgb2->tau_sampling_selection[0][0], pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1]);
    printf("#cls size = %d\n", pgb2->tau_size_cls);
    printf("#cls min/max %g/%g\n", pgb2->tau_sampling_cls[0], pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]);
    printf("#bessel_boost = %g\n", bessel_boost);
    printf("#k size = %d\n", pgb2->k_size_bessel);
    printf("#k min/max %g/%g\n", pgb2->k_bessel[0], pgb2->k_bessel[pgb2->k_size_bessel-1]);
    for(int index_type_first = 1; index_type_first < pgb2->type_size; index_type_first++){
      // Alert only looping over the first index in type_second
      for(int index_type_second = 0; index_type_second < /*index_type_first+1 */1/*pgb2->type_size*/; index_type_second++){
        printf("#### Dl[%d][%d][%d][%d][%d][%d][%d]\n",  index_type_first, index_type_second, index_l, 0, 0, bin_mean_index_selection[0], bin_mean_index_selection[0]);
        for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){

          printf("%d      %g\n", ptr->l[index_l], ptr->l[index_l]*(ptr->l[index_l]+1)*pgb2->Dl[index_type_first][index_type_second][index_l][0][0][bin_mean_index_selection[0]][bin_mean_index_selection[0]]/(2*_PI_));
        }
      }
    }


    pgb2->flag = 1;

    int flag_bisp = 1;
    if (flag_bisp == 1) {
        /* code */

      int index_k_limber_test1;
      int index_k_limber_test2;
      double transfer_interp_dens_fixed;
      double transfer_interp_phipluspsi_fixed;
      double transfer_interp_dens_free;
      double transfer_interp_phipluspsi_free;
      double Pk_fixed;
      double Pk_free;
      double one_minus_z;
      double z_minus_one;
      double one_minus_one;
      int last_index_limblens;
      double bispectrum_intDellkDellPhi;
      double Zll11z;
      double Zll11z_integrated;
      double Zllz11;
      double Zllz11_integrated;
      double Zll1z1;
      double sum_z11;
      double sum_1z1;
      double sum_11z;
      double sumA, sumB;
      double sumA_test, sumB_test;
      double chi_z11;
      double chi_1z1;
      double chi_11z;
      double chi_tilde;
      double heaviside_fixed;
      double heaviside_free;
      double DlLDA;
      double Dlg4DA;
      double DlLDB1, DlLDB2, Dlg4DB1, Dlg4DB2;
      int index_of_selection;

      double * limber_int_z11;
      double * limber_int_1z1;
      double * limber_int_11z;

      int tau_size_crop = pgb2->tau_size_selection - (bin_mean_index_selection[0]+1);


      double * tau_sampling_selection_crop;
      double * w_trapz_crop;
      class_alloc(tau_sampling_selection_crop, tau_size_crop * sizeof(double), pba->error_message);
      class_alloc(w_trapz_crop, tau_size_crop * sizeof(double), pba->error_message);
      /* INSTEAD OF DEFINING SELECTION_CROP BY THE FORMULA, SIMPLY MAP THE TAU POINTS FROM THE SELECTION GRID ! */
      for (int index_tau = 0; index_tau < tau_size_crop; index_tau++) {
        tau_sampling_selection_crop[index_tau] = pgb2->tau_sampling_selection[0][bin_mean_index_selection[0]]+index_tau*(pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1] -pgb2->tau_sampling_selection[0][bin_mean_index_selection[0]])/(tau_size_crop-1);
        //tau_sampling_selection_crop[index_tau] =pgb2->tau_sampling_selection[0][1+bin_mean_index_selection[0]+index_tau];
        //printf("-------\n");
        //printf("tau_sampling_selection_crop[%d] = %g\n", index_tau, tau_sampling_selection_crop[index_tau] );
        //printf("tau_sampling_selection[%d] = %g\n", index_tau+1+bin_mean_index_selection[0], pgb2->tau_sampling_selection[0][index_tau+1+bin_mean_index_selection[0]] );


      }


      class_call(array_trapezoidal_weights(tau_sampling_selection_crop,
                                            tau_size_crop,
                                            w_trapz_crop,
                                            pgb2->error_message),
                                            ppt2->error_message,
                                            ppt2->error_message);
      double crop_sum = 0.0;
      double crop_sum2 = 0.0;
      for (int index_tau_crop = 0; index_tau_crop < tau_size_crop; index_tau_crop++) {

        double chi_crop = (tau0-pgb2->tau_sampling_selection[0][1+bin_mean_index_selection[0]+index_tau_crop]);
        double chi_crop2 = (tau0-tau_sampling_selection_crop[index_tau_crop]);
        //printf("chi_crop = %g, chi_crop2 = %g\n", chi_crop, chi_crop2);
        crop_sum+= chi_crop*w_trapz_crop[index_tau_crop];
        crop_sum2+= chi_crop2*w_trapz_crop[index_tau_crop];
        //printf("w_trapz_crop[%d] = %g\n", index_tau_crop, w_trapz_crop[index_tau_crop]);


      }
      //printf("crop_sum = %g *(%g)*\n", crop_sum, 0.5*pow(tau0-pgb2->tau_sampling_selection[0][bin_mean_index_selection[0]],2)-0.5*pow(tau0-pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1],2));
      //printf("crop_sum2 = %g *(%g)*\n", crop_sum2, 0.5*pow(tau0-pgb2->tau_sampling_selection[0][bin_mean_index_selection[0]],2)-0.5*pow(tau0-pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1],2));


      //herehere


      class_alloc(limber_int_z11, pgb2->tau_size_cls * sizeof(double), pba->error_message);
      class_alloc(limber_int_1z1, pgb2->tau_size_cls * sizeof(double), pba->error_message);
      class_alloc(limber_int_11z, pgb2->tau_size_cls * sizeof(double), pba->error_message);
      double * pvecback_limberlens;
      double * pvecback_A94;
      class_alloc(pvecback_limberlens, pba->bg_size * sizeof(double), pba->error_message);
      class_alloc(pvecback_A94, pba->bg_size * sizeof(double), pba->error_message);
      /*double * w_trapz_to_selection_mean;
      class_alloc(w_trapz_to_selection_mean, bin_mean_index_selection[0] * sizeof(double), pba->error_message);*/

      //for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      for (int index_l = 2; index_l < 3; index_l++) {
      //for (int index_l = ptr->l_size[ppt->index_md_scalars]-1; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      //for (int index_l = ptr->l_size[ppt->index_md_scalars]-2; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
        printf("#############index_l = %d (l=%d)##################\n", index_l, ptr->l[index_l]);
        for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
          printf("index_tau = %d\n", index_tau);
          class_call(background_at_tau(pba,
                                       pgb2->tau_sampling_selection[0][index_tau],
                                       pba->long_info,
                                       pba->inter_normal,
                                       &last_index_limblens,
                                       pvecback_A94),
                                       pba->error_message,
                                       pgb2->error_message);

          double z = pba->a_today/pvecback_A94[pba->index_bg_a]-1.;
          sumA = 0.0;
          sumA_test = 0.0;
          //printf("*** pgb2->tau_sampling_selection[0][%d] = %g\n", index_tau, pgb2->tau_sampling_selection[0][index_tau]);
          //for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(index_tau+1)-1; index_tau_bessel++) {
          for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*((pgb2->tau_size_selection-1)-index_tau+1); index_tau_bessel++) {
            /*if( index_tau_bessel == pgb2->tau_size_selection-1) {
              /* skip the iteration */

              //continue;
            //}
            //printf("A: index_tau_bessel = %d\n", index_tau_bessel);

            //chi_tilde = tau0-pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel];
            //double tau_tilde = pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel];

            chi_tilde = tau0-pgb2->tau_sampling_selection_hires[0][index_tau][index_tau_bessel];
            double tau_tilde = pgb2->tau_sampling_selection_hires[0][index_tau][index_tau_bessel];
            //printf("pgb2->tau_sampling_selection_hires%d][%d] = %g\n", index_tau, index_tau_bessel, tau_tilde);

            index_of_tau_sampling_selection(tau_tilde,
                            0,
                            &index_of_selection,
                            pgb2);

            if (index_of_selection == 0) {
              DlLDA = pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][0][bin_mean_index_selection[0]];
              Dlg4DA = pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][0][bin_mean_index_selection[0]];
            }

            else{
              /* Interpolate the Dirac- angular power spectra for the specific tau found by looping the tau_sampling_bessel2 grid */
              DlLDA = pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][index_of_selection-1][bin_mean_index_selection[0]]*(pgb2->tau_sampling_selection[0][index_of_selection]-tau_tilde)
                             +pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][index_of_selection][bin_mean_index_selection[0]]*(tau_tilde-pgb2->tau_sampling_selection[0][index_of_selection-1]);
              DlLDA /= (pgb2->tau_sampling_selection[0][index_of_selection] - pgb2->tau_sampling_selection[0][index_of_selection-1]);

              Dlg4DA = pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][index_of_selection-1][bin_mean_index_selection[0]]*(pgb2->tau_sampling_selection[0][index_of_selection]-tau_tilde)
                            +pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][index_of_selection][bin_mean_index_selection[0]]*(tau_tilde-pgb2->tau_sampling_selection[0][index_of_selection-1]);
              Dlg4DA /= (pgb2->tau_sampling_selection[0][index_of_selection] - pgb2->tau_sampling_selection[0][index_of_selection-1]);
            }



            sumA +=(2./chi_tilde)*DlLDA*Dlg4DA*w_trapz_lens_bessel2[index_tau][index_tau_bessel];
            //sumA_test += (2./chi_tilde)*w_trapz_lens_bessel2[index_tau][index_tau_bessel];
            sumA_test += (2./chi_tilde)*w_trapz_selection_hires[0][index_tau][index_tau_bessel];

          }
          double exact = 2.*log(tau0-pgb2->tau_sampling_selection[0][index_tau])-2.*log(tau0-pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1]);
          printf("sumA_test = %g, chi = %g, tau = %g, should be = %g, ERROR = *%g %%*\n", sumA_test, tau0-pgb2->tau_sampling_selection[0][index_tau], pgb2->tau_sampling_selection[0][index_tau], exact, (sumA_test - exact)*100/exact );
          double weightB = ((tau0-pgb2->tau_sampling_selection[0][bin_mean_index_selection[0]])-(tau0-pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1]))/((pgb2->tau_size_selection-bin_mean_index_cls[0]+1)-1);
          double sum_test = 0.0;
        //  printf("bin_mean_index_selection[0] = %d\n", bin_mean_index_selection[0] );

          /*for (int i = 0; i < pgb2->tau_size_selection; i++) {
            printf("tau_sampling_selection[0][%i] = %g\n", i, pgb2->tau_sampling_selection[0][i]);
          }*/

          sumB = 0.0;
          /*sumB = (2./tau0-pgb2->tau_sampling_bessel2[index_tau][0])*(pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][pgb2->tau_size_selection-1][index_tau]
                                *pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][pgb2->tau_size_selection-1][bin_mean_index_selection[0]]
                                +pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][pgb2->tau_size_selection-1][bin_mean_index_selection[0]]
                                *pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][pgb2->tau_size_selection-1][index_tau])
                                *0.5*weightB;*/
          sumB_test = 0.0;


          //sumB_test = 0.5*weightB*(tau0-pgb2->tau_sampling_selection[0][bin_mean_index_selection[0]]);
          for (int index_tau_crop = 0; index_tau_crop < tau_size_crop; index_tau_crop++) {
            //printf("index_tau_crop = %d\n", index_tau_crop);
            chi_tilde = (tau0-tau_sampling_selection_crop[index_tau_crop]);
            double tau_tilde = tau_sampling_selection_crop[index_tau_crop];
            double chi_tilde_interp;
            index_of_tau_sampling_selection(tau_tilde,
                            0,
                            &index_of_selection,
                            pgb2);

            if (index_of_selection == 0) {
              DlLDB1 = pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][0][index_tau];
              DlLDB2 = pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][0][bin_mean_index_selection[0]];
              Dlg4DB1 = pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][0][index_tau];
              Dlg4DB2 = pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][0][bin_mean_index_selection[0]];
            }

            else{

              chi_tilde_interp = (tau0-pgb2->tau_sampling_selection[0][index_of_selection-1])*(pgb2->tau_sampling_selection[0][index_of_selection]-tau_tilde)
                             +(tau0-pgb2->tau_sampling_selection[0][index_of_selection])*(tau_tilde-pgb2->tau_sampling_selection[0][index_of_selection-1]);
              chi_tilde_interp /= (pgb2->tau_sampling_selection[0][index_of_selection] - pgb2->tau_sampling_selection[0][index_of_selection-1]);

              DlLDB1 = pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][index_of_selection-1][index_tau]*(pgb2->tau_sampling_selection[0][index_of_selection]-tau_tilde)
                             +pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][index_of_selection][index_tau]*(tau_tilde-pgb2->tau_sampling_selection[0][index_of_selection-1]);
              DlLDB1 /= (pgb2->tau_sampling_selection[0][index_of_selection] - pgb2->tau_sampling_selection[0][index_of_selection-1]);

              DlLDB2 = pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][index_of_selection-1][bin_mean_index_selection[0]]*(pgb2->tau_sampling_selection[0][index_of_selection]-tau_tilde)
                             +pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][index_of_selection][bin_mean_index_selection[0]]*(tau_tilde-pgb2->tau_sampling_selection[0][index_of_selection-1]);
              DlLDB2 /= (pgb2->tau_sampling_selection[0][index_of_selection] - pgb2->tau_sampling_selection[0][index_of_selection-1]);

              Dlg4DB1 = pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][index_of_selection-1][index_tau]*(pgb2->tau_sampling_selection[0][index_of_selection]-tau_tilde)
                             +pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][index_of_selection][index_tau]*(tau_tilde-pgb2->tau_sampling_selection[0][index_of_selection-1]);
              Dlg4DB1 /= (pgb2->tau_sampling_selection[0][index_of_selection] - pgb2->tau_sampling_selection[0][index_of_selection-1]);

              Dlg4DB2 = pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][index_of_selection-1][bin_mean_index_selection[0]]*(pgb2->tau_sampling_selection[0][index_of_selection]-tau_tilde)
                             +pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][index_of_selection][bin_mean_index_selection[0]]*(tau_tilde-pgb2->tau_sampling_selection[0][index_of_selection-1]);
              Dlg4DB2 /= (pgb2->tau_sampling_selection[0][index_of_selection] - pgb2->tau_sampling_selection[0][index_of_selection-1]);
            }


            sumB += (2./chi_tilde)*(DlLDB1*Dlg4DB2+DlLDB2*Dlg4DB1)*w_trapz_crop[index_tau_crop];
            /*sumB +=(2./chi_tilde)*(pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][1+bin_mean_index_selection[0]+index_tau_crop][index_tau]
                                  *pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][1+bin_mean_index_selection[0]+index_tau_crop][bin_mean_index_selection[0]]
                                  +pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][1+bin_mean_index_selection[0]+index_tau_crop][bin_mean_index_selection[0]]
                                  *pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][1+bin_mean_index_selection[0]+index_tau_crop][index_tau])
                  *w_trapz_crop[index_tau_crop];*/
            sumB_test += (2./chi_tilde)*w_trapz_crop[index_tau_crop];
          }

          /*for (int index_tau_bessel = bin_mean_index_selection[0]+1; index_tau_bessel < pgb2->tau_size_selection-1; index_tau_bessel++) {
            sumB_test += weightB*(tau0-pgb2->tau_sampling_selection[0][index_tau_bessel]);

            //printf("sumB_test = %g\n", sumB_test );
          }
          sumB_test+= 0.5*weightB*(tau0-pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1]);

          sumB += (2./chi_tilde)*(pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][bin_mean_index_selection[0]][index_tau]
                                *pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][bin_mean_index_selection[0]][bin_mean_index_selection[0]]
                                +pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][bin_mean_index_selection[0]][bin_mean_index_selection[0]]
                                *pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][bin_mean_index_selection[0]][index_tau])
                *0.5*weightB;*/




          //printf("sumB_test = %g *(%g)*\n", sumB_test, 0.5*pow(tau0-pgb2->tau_sampling_selection[0][bin_mean_index_selection[0]],2)-0.5*pow(tau0-pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1],2));
          printf("sumB_test = %g *(%g)*\n", sumB_test, 2.*log(tau0-pgb2->tau_sampling_selection[0][bin_mean_index_selection[0]])-2.*log(tau0-pgb2->tau_sampling_selection[0][pgb2->tau_size_selection-1]));
          //printf("tau_size_crop = %d\n", tau_size_crop );
          //printf("sumB_test final = %g (*%g*)\n", sumB_test, 0.5*(pow(tau0-pgb2->tau_sampling_selection[0][bin_mean_index_selection[0]],2)));
        //  for (int index_tau_bessel = 1; index_tau_bessel < bin_mean_index_selection[0]; index_tau_bessel++) {
            /*if( index_tau_bessel == index_tau) {
              /* skip the iteration */

          //    continue;
           //}
            //chi_tilde = tau0-pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel];
          //  chi_tilde = tau0-pgb2->tau_sampling_selection[0][index_tau_bessel];

            //sumB +=(2./chi_tilde)*(pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][index_tau_bessel][index_tau]
                      //            *pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][index_tau_bessel][bin_mean_index_selection[0]]
                        //          +pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][index_tau_bessel][bin_mean_index_selection[0]]
                          //        *pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][index_tau_bessel][index_tau])
              //    *weightB;

            /*if (z >5.14) {
              printf("----------\n");
              printf("sumB = %g\n", sumB);
              printf("2./chi_tilde = %g\n", 2./chi_tilde );
              //printf("tau_tilde = %g\n", );
              printf("index_tau_bessel = %d\n", index_tau_bessel);
              /*printf("B: Dl1 = %g\n", pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][index_tau_bessel][index_tau]);
              printf("B: Dl2 = %g\n", pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][index_tau_bessel][bin_mean_index_selection[0]]);
              printf("B: Dl3 = %g\n", pgb2->Dl[pgb2->index_type_lens][pgb2->index_type_density][index_l][0][0][index_tau_bessel][bin_mean_index_selection[0]]);
              printf("B: Dl4 = %g\n", pgb2->Dl[pgb2->index_type_g4][pgb2->index_type_density][index_l][0][0][index_tau_bessel][index_tau]);
              printf("pgb2->index_type_lens = %d\n", pgb2->index_type_lens);
              printf("pgb2->index_type_g4 = %d\n", pgb2->index_type_g4);
              printf("pgb2->index_type_density = %d\n", pgb2->index_type_density);*/
            /*  printf("w_trapz[%d][%d] = %g\n", index_tau, index_tau_bessel, w_trapz_lens_bessel2[index_tau][index_tau_bessel]);
              printf("weightB =%g\n", weightB );
            }*/
        //  }




          double min_limberI, min_limberII, min_limberIII;
          double max_limberI, max_limberII, max_limberIII;
          double A, B, C;

          class_call(drc3jj (ptr->l[index_l],
                             ptr->l[index_l],
                             1,
                             -1,
                             &min_limberI,
                             &max_limberI,
                             &A,
                             1000,
                             pgb2->error_message),
                             pgb2->error_message,
                             pgb2->error_message);

          class_call(drc3jj (ptr->l[index_l],
                             ptr->l[index_l],
                             -1,
                             1,
                             &min_limberII,
                             &max_limberII,
                             &B,
                             1000,
                             pgb2->error_message),
                             pgb2->error_message,
                             pgb2->error_message);

          class_call(drc3jj (ptr->l[index_l],
                             ptr->l[index_l],
                             0,
                             0,
                             &min_limberIII,
                             &max_limberIII,
                             &C,
                             1000,
                             pgb2->error_message),
                             pgb2->error_message,
                             pgb2->error_message);

          double A_test;

          /*Al1l2l3(l,
                  l,
                  &A_test,
                  pgb2);*/




          double A_LLL = 0.5*(A+B)/C;
          /*Al1l2l3(ptr->l[index_l],
                  ptr->l[index_l],
                  &A_LLL,
                  pgb2);*/
          printf("%g      %g\n", z, A_LLL*(ptr->l[index_l]*(ptr->l[index_l]+1))*(sumA+sumB));
        }
      }
      exit(0);

      for (int l = 382; l < 383; l++) {
        printf("#############\n");
        index_k_limber_test1 = 0;
        index_k_limber_test2 = 0;
        for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {
          //printf("tau_sampling_cls[%d] = %g\n", index_tau, pgb2->tau_sampling_cls[index_tau] );
          //printf("tau_sampling_selection[0][%d] = %g\n", index_tau, pgb2->tau_sampling_selection[0][index_tau] );
          class_call(background_at_tau(pba,
                                       pgb2->tau_sampling_cls[index_tau],
                                       pba->long_info,
                                       pba->inter_normal,
                                       &last_index_limblens,
                                       pvecback_limberlens),
                                       pba->error_message,
                                       pgb2->error_message);

         /*class_call(background_at_tau(pba,
                                      pgb2->tau_sampling_cls[index_tau],
                                      pba->long_info,
                                      pba->inter_normal,
                                      &last_index_g,
                                      pvecbackg),
                                      pba->error_message,
                                      pgb2->error_message);*/

          z = pba->a_today/pvecback_limberlens[pba->index_bg_a]-1.;

          double k_fixed = (l+0.5)/(tau0-pgb2->tau_sampling_cls[bin_mean_index_cls[0]]);
          double chi_fixed = (tau0-pgb2->tau_sampling_cls[bin_mean_index_cls[0]]);

          double k_free = (l+0.5)/(tau0-pgb2->tau_sampling_cls[index_tau]);
          double chi_free = (tau0-pgb2->tau_sampling_cls[index_tau]);

          class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, k_fixed, &Pk_fixed), ppm->error_message, pgb2->error_message);

          //printf("Pk_fixed b4 = %g\n", Pk_fixed );
          /* Apply normalisation to power spectra */
          double Pk_fixed_norm = 2*_PI_*_PI_/k_fixed/k_fixed/k_fixed;
          //printf("Pk_fixed norm = %g\n", Pk_fixed_norm );
          class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, k_free, &Pk_free), ppm->error_message, pgb2->error_message);
          //printf("Pk_free  = %g\n", Pk_free );
          double Pk_free_norm= 2*_PI_*_PI_/k_free/k_free/k_free;
          //double Pk_free_norm= 0.;
          //printf("Pk_free norm = %g\n", Pk_free_norm );

          class_call(index_of_k_bessel(k_fixed, &index_k_limber_test1, pgb2), pgb2->error_message, pgb2->error_message);
          class_call(index_of_k_bessel(k_free, &index_k_limber_test2, pgb2), pgb2->error_message, pgb2->error_message);

          double k_found_minus1 = pgb2->k_bessel[index_k_limber_test1];
          double k_found_plus1 = pgb2->k_bessel[index_k_limber_test1+1];
          double k_found_minus2 = pgb2->k_bessel[index_k_limber_test2];
          double k_found_plus2 = pgb2->k_bessel[index_k_limber_test2+1];

          //printf("k_minus_fixed = %g, k_fixed = %g, k_plus_fixed = %g\n", k_found_minus1, k_fixed, k_found_plus1);
          //printf("k_minus_free = %g, k_free = %g, k_plus_free = %g\n", k_found_minus2, k_free, k_found_plus2);

          transfer_interp_dens_fixed = pgb2->first_order_sources[pgb2->index_source_delta_cdm][bin_mean_index_cls[0]][index_k_limber_test1]*(pgb2->k_bessel[index_k_limber_test1+1]-k_fixed)
                                  + pgb2->first_order_sources[pgb2->index_source_delta_cdm][bin_mean_index_cls[0]][index_k_limber_test1+1]*(k_fixed-pgb2->k_bessel[index_k_limber_test1]);
          transfer_interp_dens_fixed /= (pgb2->k_bessel[index_k_limber_test1+1] - pgb2->k_bessel[index_k_limber_test1]);

          transfer_interp_phipluspsi_fixed = pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][bin_mean_index_cls[0]][index_k_limber_test1]*(pgb2->k_bessel[index_k_limber_test1+1]-k_fixed)
                                  + pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][bin_mean_index_cls[0]][index_k_limber_test1+1]*(k_fixed-pgb2->k_bessel[index_k_limber_test1]);
          transfer_interp_phipluspsi_fixed /= (pgb2->k_bessel[index_k_limber_test1+1] - pgb2->k_bessel[index_k_limber_test1]);

          transfer_interp_dens_free = pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_limber_test2]*(pgb2->k_bessel[index_k_limber_test2+1]-k_free)
                                  + pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_limber_test2+1]*(k_free-pgb2->k_bessel[index_k_limber_test2]);
          transfer_interp_dens_free /= (pgb2->k_bessel[index_k_limber_test2+1] - pgb2->k_bessel[index_k_limber_test2]);

          transfer_interp_phipluspsi_free = pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][index_tau][index_k_limber_test2]*(pgb2->k_bessel[index_k_limber_test2+1]-k_free)
                                  + pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][index_tau][index_k_limber_test2+1]*(k_free-pgb2->k_bessel[index_k_limber_test2]);
          transfer_interp_phipluspsi_free /= (pgb2->k_bessel[index_k_limber_test2+1] - pgb2->k_bessel[index_k_limber_test2]);

          heaviside(1.-z, &one_minus_z);
          heaviside(z-1., &z_minus_one);
          heaviside(1.-1.,&one_minus_one);

          sum_z11 = 0.0;
          sum_1z1 = 0.0;
          sum_11z = 0.0;

          /*interpolate_Dl_for_l(l2,
                               &type_A,
                               &type_C,
                               &bin1,
                               &bin2,
                               &last_index_l_DAC12,
                               &index_tau_first,
                               &index_tau_second,
                               &DAC12,
                               pgb2,
                               ptr,
                               ppt);*/


          for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(index_tau+1)-1; index_tau_bessel++) {
            chi_z11 = tau0-pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel];
            heaviside(chi_z11-chi_fixed, &heaviside_fixed);
            //heaviside(chi_z11-chi_fixed, heaviside1);




            sum_z11 += (chi_z11-chi_fixed)*heaviside_fixed*heaviside_fixed*w_trapz_lens_bessel2[index_tau][index_tau_bessel]/chi_z11/chi_z11/chi_z11/chi_fixed;
          }
          limber_int_z11[index_tau] = sum_z11;
          //printf("limber_int_z11[%d] =%g\n", index_tau,limber_int_z11[index_tau]);
          for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(bin_mean_index_cls[0]+1)-1; index_tau_bessel++) {

            chi_1z1 = tau0-pgb2->tau_sampling_bessel2[bin_mean_index_cls[0]][index_tau_bessel];
            chi_11z = tau0-pgb2->tau_sampling_bessel2[bin_mean_index_cls[0]][index_tau_bessel];


            heaviside(chi_11z-chi_fixed, &heaviside_fixed);
            heaviside(chi_11z-chi_free, &heaviside_free);



            double step = (chi_1z1-chi_free)*heaviside_fixed*heaviside_free*w_trapz_lens_bessel2[bin_mean_index_cls[0]][index_tau_bessel]/chi_1z1/chi_1z1/chi_1z1/chi_free;
            //if (chi_free < chi_fixed && step != 0.) {

            sum_1z1 += (chi_1z1-chi_free)*heaviside_fixed*heaviside_free*w_trapz_lens_bessel2[bin_mean_index_cls[0]][index_tau_bessel]/chi_1z1/chi_1z1/chi_1z1/chi_free;
            sum_11z += (chi_11z-chi_fixed)*heaviside_fixed*heaviside_free*w_trapz_lens_bessel2[bin_mean_index_cls[0]][index_tau_bessel]/chi_11z/chi_11z/chi_11z/chi_fixed;
            if (index_tau == 481) {
              printf("index_tau_bessel = %d\n", index_tau_bessel);
              printf("**chi_1z1-chi_free) = %g\n", chi_1z1-chi_free );

              printf("chi_1z1 = %g\n", chi_1z1 );
              printf("2. heaviside_fixed =%g\n", heaviside_fixed );
              printf("2. heaviside_free =%g\n", heaviside_free );
              printf("step = %g\n", step);
              printf("sum_1z1 = %g\n", sum_1z1);
            }
          }
          limber_int_1z1[index_tau] = sum_1z1;
          limber_int_11z[index_tau] = sum_11z;
          /*printf("chi_free = %g\n", chi_free);
          printf("limber_int_1z1[%d] =%g\n", index_tau,limber_int_1z1[index_tau]);
          printf("limber_int_11z[%d] =%g\n", index_tau,limber_int_11z[index_tau]);*/

          Zllz11_integrated = -_PI_*_PI_
                              /4.
                              /pow(chi_fixed,4)
                              *Pk_fixed
                              *Pk_fixed_norm
                              *Pk_fixed
                              *Pk_fixed_norm
                              *transfer_interp_dens_fixed
                              *transfer_interp_dens_fixed
                              *transfer_interp_phipluspsi_fixed
                              *transfer_interp_phipluspsi_fixed
                              *2
                              *limber_int_z11[index_tau]
                              *sqrt(l*(l+1))
                              *pow(l*(l+1),3/2);

          Zll11z_integrated = -_PI_*_PI_
                              /4.
                              /chi_fixed
                              /chi_fixed
                              /chi_free
                              /chi_free
                              *Pk_fixed
                              *Pk_fixed_norm
                              *Pk_free
                              *Pk_free_norm
                              *transfer_interp_dens_fixed
                              *transfer_interp_dens_free
                              *transfer_interp_phipluspsi_fixed
                              *transfer_interp_phipluspsi_free
                              *(limber_int_1z1[index_tau]+limber_int_11z[index_tau])
                              *sqrt(l*(l+1))
                              *pow(l*(l+1),3/2);



          // r1 is free, r2=r3= fixed at z=1
          //one_minus_z);
          //heaviside(z-1., &z_minus_one);
          //heaviside(1.-1.,&one_minus_one);
          double factor_z11A = (chi_free-chi_fixed)/(2*pow(chi_free,2)*pow(chi_fixed,2));
          double bracket_z11A1 = chi_free-chi_fixed;
          double bracket_z11A2 = (2*chi_free*chi_fixed-chi_free*chi_fixed-chi_fixed*chi_fixed)/chi_fixed;
          double subsetA_z11 = z_minus_one*one_minus_one*factor_z11A*(bracket_z11A1+bracket_z11A2);
          double subsetA_z11_old = z_minus_one*one_minus_one*((chi_free-chi_fixed)/(2*pow(chi_free,2)*pow(chi_fixed,2))*(((chi_free-chi_fixed))+(2*chi_free*chi_fixed-chi_free*chi_fixed-chi_fixed*chi_fixed)/chi_fixed));
          double subsetB_z11 = z_minus_one*one_minus_one*((chi_free-chi_fixed)/(2*pow(chi_free,2)*pow(chi_fixed,2))*(((chi_free-chi_fixed))+(2*chi_free*chi_fixed-chi_free*chi_fixed-chi_fixed*chi_fixed)/chi_fixed));

          Zllz11 = _PI_*_PI_
                   /4.
                   /pow(chi_fixed, 4)
                   *Pk_fixed_norm*Pk_fixed
                   *Pk_fixed_norm*Pk_fixed
                   *transfer_interp_phipluspsi_fixed
                   *transfer_interp_phipluspsi_fixed
                   *transfer_interp_dens_fixed
                   *transfer_interp_dens_fixed
                   *sqrt(l*(l+1))
                   *pow(l*(l+1),3/2.)
                   *(-subsetA_z11-subsetB_z11);

         double factor_1z1A = (chi_fixed - chi_free)/(2*pow(chi_free,2)*pow(chi_fixed,2));
         double bracket_1z1A1 = chi_fixed - chi_free;
         double bracket_1z1A2 = (2*chi_free*chi_fixed - chi_free*chi_fixed - chi_fixed*chi_fixed)/chi_fixed;
         double subsetA_1z1 = one_minus_z*z_minus_one*factor_1z1A*(bracket_1z1A1+bracket_1z1A2);
         //double subsetA_1z1 = one_minus_z*z_minus_one*((chi_fixed - chi_free)/(2*pow(chi_free,2)*pow(chi_fixed,2))*(((chi_fixed - chi_free))+(2*chi_free*chi_fixed - chi_free*chi_fixed - chi_fixed*chi_fixed)/chi_fixed));
         double subsetB_1z1 = one_minus_one*one_minus_z*(((chi_fixed-chi_fixed)/(2*pow(chi_fixed,4)))*(((chi_fixed-chi_fixed))+(2*chi_fixed*chi_fixed-chi_free*chi_fixed-chi_free*chi_fixed)/chi_free));

         Zll1z1 = _PI_*_PI_
                  /4.
                  /chi_free
                  /chi_free
                  /chi_fixed
                  /chi_fixed
                  *Pk_free_norm*Pk_free
                  *Pk_fixed_norm*Pk_fixed
                  *transfer_interp_phipluspsi_fixed
                  *transfer_interp_phipluspsi_free
                  *transfer_interp_dens_free
                  *transfer_interp_dens_fixed
                  *sqrt(l*(l+1))
                  *pow(l*(l+1),3/2.)
                  *(-subsetA_1z1-subsetB_1z1);

         /*printf("1/chi_free = %g\n", 1./chi_free);
         printf("1/chi_fixed = %g\n", 1./chi_fixed);
         printf("Pk_free = %g\n", Pk_free);
         printf("Pk_fixed = %g\n", Pk_fixed);
         printf("phi_plus_psi_free = %g\n", transfer_interp_phipluspsi_free);
         printf("phi_plus_psi_fixed = %g\n", transfer_interp_phipluspsi_fixed);
         printf("dens_fixed = %g\n", transfer_interp_dens_fixed);
         printf("dens_free = %g\n", transfer_interp_dens_free);
         printf("sqrt(l*(l+1)) = %g\n",sqrt(l*(l+1)));
         printf("*pow(l*(l+1),3/2.) = %g", pow(l*(l+1),3/2.));*/


         double subsetA_11z = one_minus_z*z_minus_one*((chi_fixed-chi_free)/(2*pow(chi_free,2)*pow(chi_fixed,2))*(((chi_fixed-chi_free))+(2*chi_free*chi_fixed-chi_free*chi_fixed-chi_fixed*chi_fixed)/chi_fixed));
         double subsetB_11z = one_minus_one*one_minus_z*(((chi_fixed-chi_fixed)/(2*pow(chi_fixed,4)))*(((chi_fixed-chi_fixed))+(2*chi_fixed*chi_fixed-chi_free*chi_fixed-chi_free*chi_fixed)/chi_free));

         Zll11z = _PI_*_PI_
                  /4.
                  /chi_free
                  /chi_free
                  /chi_fixed
                  /chi_fixed
                  *Pk_free_norm*Pk_free
                  *Pk_fixed_norm*Pk_fixed
                  *transfer_interp_phipluspsi_fixed
                  *transfer_interp_phipluspsi_free
                  *transfer_interp_dens_free
                  *transfer_interp_dens_fixed
                  *sqrt(l*(l+1))
                  *pow(l*(l+1),3/2.)
                  *(-subsetA_1z1-subsetB_1z1);

          /*printf("subsetA_1z1 = %g\n", subsetA_1z1);
          printf("subsetB_1z1 = %g\n", subsetB_1z1);
          printf("subsetA_z11 = %g\n", subsetA_z11);
          printf("subsetA_z11_old = %g\n", subsetA_z11_old);
          printf("subsetB_z11 = %g\n", subsetB_z11);
          printf("subsetA_11z = %g\n", subsetA_11z);
          printf("subsetB_11z = %g\n", subsetB_11z);*/





          /*Zll11z = -1.
                    *one_minus_z
                    *z_minus_one
                    *_PI_*_PI_
                    /4
                    /(chi_fixed*chi_fixed)
                    /(chi_free*chi_free)
                    *Pk_free*Pk_free_norm
                    *Pk_fixed*Pk_fixed_norm
                    *transfer_interp_dens_fixed
                    *transfer_interp_dens_free
                    *transfer_interp_phipluspsi_fixed
                    *transfer_interp_phipluspsi_free
                    *(chi_fixed-chi_free);
                    /*(2*chi_fixed*chi_fixed*chi_free*chi_free)
                    *pow(l*(l+1),3/2.)
                    *sqrt(l*(l+1))
                    *((2*chi_free*chi_fixed-chi_fixed*(chi_free+chi_fixed))/chi_fixed)+(chi_fixed-chi_free);*/


          /*Zllz11 = -2.
                    *one_minus_one
                    *z_minus_one
                    *_PI_*_PI_
                    /4
                    /(chi_fixed*chi_fixed)
                    /(chi_fixed*chi_fixed)
                    *Pk_fixed*Pk_fixed_norm
                    *Pk_fixed*Pk_free_norm
                    *transfer_interp_dens_fixed
                    *transfer_interp_dens_fixed
                    *transfer_interp_phipluspsi_fixed
                    *transfer_interp_phipluspsi_fixed
                    *(chi_free-chi_fixed);
                    /*(2*chi_fixed*chi_fixed*chi_free*chi_free)
                    *pow(l*(l+1),3/2.)
                    *sqrt(l*(l+1))
                    *((2*chi_free*chi_fixed-chi_fixed*(chi_free+chi_fixed))/chi_fixed)+(chi_free-chi_fixed);*/

          //printf("one_minus_z = %g\n", one_minus_z);
          //printf("z_minus_one = %g\n", z_minus_one);
          //printf("chi_fixed = %g\n", chi_fixed );
          //printf("chi_free = %g\n", chi_free );
          //printf("Pk_fixed = %g\n", Pk_fixed );
          //printf("Pk_free = %g\n", Pk_free );
          //printf("%g, %g, %g, %g\n", transfer_interp_dens_fixed, transfer_interp_phipluspsi_fixed, transfer_interp_dens_free, transfer_interp_phipluspsi_free);


          double min_limberI, min_limberII, min_limberIII;
          double max_limberI, max_limberII, max_limberIII;
          double A_limberI, A_limberII, A_limberIII;

          class_call(drc3jj (l,
                             l,
                             1,
                             -1,
                             &min_limberI,
                             &max_limberI,
                             &A_limberI,
                             1000,
                             pgb2->error_message),
                             pgb2->error_message,
                             pgb2->error_message);

          class_call(drc3jj (l,
                             l,
                             -1,
                             1,
                             &min_limberII,
                             &max_limberII,
                             &A_limberII,
                             1000,
                             pgb2->error_message),
                             pgb2->error_message,
                             pgb2->error_message);

          class_call(drc3jj (l,
                             l,
                             0,
                             0,
                             &min_limberIII,
                             &max_limberIII,
                             &A_limberIII,
                             1000,
                             pgb2->error_message),
                             pgb2->error_message,
                             pgb2->error_message);

          double A_test;

          /*Al1l2l3(l,
                  l,
                  &A_test,
                  pgb2);*/




          double A_limber = 0.5*(A_limberI+A_limberII)/A_limberIII;

          printf("A_limber, A_test = %g, %g\n", A_limber, A_test);

          bispectrum_intDellkDellPhi = A_limber*(2*Zll11z+Zllz11);
          double bispectrum_intDellkDellPhi_integrated = A_limber*(2*Zll11z_integrated+Zllz11_integrated);

          double bispectrum_limber = 4.*A_limber*(Zllz11+Zll1z1+Zll11z)/_PI_/_PI_;

          /*printf("Zllz11 = %g\n",Zllz11 );
          printf("Zll1z1 = %g\n",Zll1z1 );
          printf("Zll11z = %g\n",Zll11z );
          printf("A_limber = %g\n", A_limber );*/

          //printf("original: %g      %d      %g\n", z, l, bispectrum_intDellkDellPhi);
          printf(" %g      %d      %g\n", z, l, bispectrum_intDellkDellPhi_integrated);
          //printf("%g      %g\n", Zll11z, Zllz11);
        }
      }
    }
    exit(0);






  /* CLASS STYLE ANGULAR POWER SPECTRUM INTEGRATING TIME FIRST */
  /*printf("Entering CLASS style -transfer- preparation.. \n");
  double source_interp;
  int index_of_cls;
  double tau_one_class;
  for(int index_type = 0; index_type < pgb2->type_size; index_type++){
    for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
      for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
        printf("index_l = %d (l=%d)\n", index_l, ptr->l[index_l]);
        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
          double tau_sum = 0.0;
          /* We loop over the tau-selection indices and interpolate on the tau_cls grid */
        /*  for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {

            tau_one_class = pgb2->tau_sampling_selection[bin1][index_tau];

            index_of_tau_sampling_cls(tau_one_class, &index_of_cls, pgb2);


            source_interp = pgb2->first_order_sources_integ[index_type][index_l][index_of_cls-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_of_cls]-tau_one_class)
                          + pgb2->first_order_sources_integ[index_type][index_l][index_of_cls][index_k_bessel]*(tau_one_class-pgb2->tau_sampling_cls[index_of_cls-1]);
            source_interp /= (pgb2->tau_sampling_cls[index_of_cls] - pgb2->tau_sampling_cls[index_of_cls-1]);



            tau_sum += source_interp*selection[bin1][index_tau]*w_trapz[bin1][index_tau];
            }

          class_xfer[index_type][bin1][index_l][index_k_bessel] = tau_sum;
          }
        }
      }
    }
    /*printf("Entering CLASS style k-integration... \n");
    double Pkk;
    for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
      printf("index_type_first = %d\n", index_type_first );
      for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
        printf("index_type_second = %d\n", index_type_second );
        for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
          printf("bin1 = %d\n", bin1 );
          for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
            printf("bin2 = %d\n", bin2 );
            for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
              printf("index_l = %d (l=%d)\n", index_l, ptr->l[index_l]);
              double k_sum = 0.;
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
                class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, pgb2->k_bessel[index_k_bessel], &Pkk), ppm->error_message, pgb2->error_message);

                k_sum += 4. * _PI_ *  Pkk * pow(pgb2->k_bessel[index_k_bessel],-1.0)* class_xfer[index_type_first][bin1][index_l][index_k_bessel]
                  * class_xfer[index_type_second][bin2][index_l][index_k_bessel]*pgb2->w_trapz_k[index_k_bessel];
              }
              pgb2->Cl3[index_type_first][index_type_second][index_l][bin1][bin2] = k_sum;
              printf("CLASS Style: pgb2->Cl3[%d][%d][%d][%d][%d] = %g\n",
                      index_type_first,
                      index_type_second,
                      index_l,
                      bin1,
                      bin2,
                      ptr->l[index_l]*(ptr->l[index_l]+1)*pgb2->Cl3[index_type_first][index_type_second][index_l][bin1][bin2]/(2*_PI_));
            }
          }
        }
      }
    }*/
    printf("# CLASS STYLE ANGPOWSPEC COMPUTATION (pgb2->Dl3, Cl3)\n" );
    if (k_sampling == 0) {
      printf("#k is linearly sampled with %d points\n", pgb2->k_size_bessel );
    }
    if (k_sampling == 1) {
      printf("#k is log-sampled with %d points\n", pgb2->k_size_bessel );
      //printf("#k_logbase = %g\n",k_logbase);
    }

    printf("#k_max = %g\n",pgb2->k_bessel[pgb2->k_size_bessel-1]);
    printf("#tau_size_selection = %d\n", pgb2->tau_size_selection);
    printf("#tau_size_cls= %d\n", pgb2->tau_size_cls);
    printf("#bessel_boost = %d\n", bessel_boost );
    printf("#z1=%g, z2=%g (%g) Gaussian\n", ppt->selection_mean[0],  ppt->selection_mean[0],ppt->selection_width[0]);
    printf("#l    l(l+1)Cl3/2pi   \n");
    //for(int index_type = 0; index_type < pgb2->type_size; index_type++){
    for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
      for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
        for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
        //for (int index_l = ptr->l_size[ppt->index_md_scalars]-2; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
        //for (int index_l = 3; index_l < 4; index_l++) {
          printf("%d    %d    %g\n",index_l, ptr->l[index_l],
                              ptr->l[index_l]*(ptr->l[index_l]+1)*(pgb2->Cl3[0][0][index_l][0][0])/(2*_PI_));
        }
      }
    }


  //  }
    int index_l_found;
    int last_index_l = 0;
    int last_index_l2 = 0;
    int type_A;
    int type_B;
    int type_C;
    int type_D;

    type_A = pgb2->index_type_lens;
    type_B = pgb2->index_type_g4;
    type_C = pgb2->index_type_density;
    type_D = pgb2->index_type_density;
    double spline;
    /* Codes for Dl DAC12 = Dl for type A and C on bins 1 and 2 */
    double DAC12, DBD13, DBC12, DAD13, DAC21, DBD23, DBC21, DAD23, DAC32, DBD31, DBC32, DAD31;

    int last_index_l_DAC12;
    int last_index_l_DBD13;
    int last_index_l_DBC12;
    int last_index_l_DAD13;

    int last_index_l_DAC21;
    int last_index_l_DBD23;
    int last_index_l_DBC21;
    int last_index_l_DAD23;

    int last_index_l_DAC32;
    int last_index_l_DBD31;
    int last_index_l_DBC32;
    int last_index_l_DAD31;
    printf("#Printing quadratic Bispectrum..\n");
    printf("#type_A = %s", getName(type_A));
    printf("#type_B = %s", getName(type_B));
    printf("#type_C = %s", getName(type_C));
    printf("#type_D = %s", getName(type_D));
    printf("#g_bias = %g\n", g_bias );

    int l1,l2,l3;
    double * pvecback_z3;
    class_alloc(pvecback_z3, pba->bg_size * sizeof(double), pba->error_message);
    for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      printf("##### l = ptr->l[%d] = %d\n", index_l, ptr->l[index_l]);
      l1=ptr->l[index_l];
      l2=ptr->l[index_l];
      l3=ptr->l[index_l];
      for (int bin1 = 0; bin1 < 1; bin1++) {
        int index_tau_first_min = bin_mean_index_selection[bin1];
        for (int bin2 = 0; bin2 < 1; bin2++) {
          int index_tau_second_min = bin_mean_index_selection[bin2];
          for (int bin3 = 0; bin3 < 1; bin3++) {
            for (int index_tau_first = index_tau_first_min; index_tau_first < index_tau_first_min+1; index_tau_first++) {
              for (int index_tau_second = index_tau_second_min; index_tau_second < index_tau_second_min+1; index_tau_second++) {
                  for (int index_tau_third = 0; index_tau_third < pgb2->tau_size_selection; index_tau_third++) {
                    double z3;
                    class_call(background_at_tau(pba,
                                                 pgb2->tau_sampling_selection[bin3][index_tau_third],
                                                 pba->long_info,
                                                 pba->inter_normal,
                                                 &last_index,
                                                 pvecback_z3),
                                                 pba->error_message,
                                                 pgb2->error_message);

                      /*infer redhsift */
                    z3 = pba->a_today/pvecback_z3[pba->index_bg_a]-1.;
                    last_index_l_DAC12 = 0;
                    last_index_l_DBD13 = 0;
                    last_index_l_DBC12 = 0;
                    last_index_l_DAD13 = 0;

                    last_index_l_DAC21 = 0;
                    last_index_l_DBD23 = 0;
                    last_index_l_DBC21 = 0;
                    last_index_l_DAD23 = 0;

                    last_index_l_DAC32 = 0;
                    last_index_l_DBD31 = 0;
                    last_index_l_DBC32 = 0;
                    last_index_l_DAD31 = 0;
                  /*for(int l = 2; l < 500; l++){
                    l1=l;
                    l2=l;
                    l3=l;*/
                  //  printf("last_index_l_DBC13, last_index_l_DAD12, last_index_l_DBD13, last_index_l_DAC12, last_index_l_DBC23, last_index_l_DAD21 = %d, %d, %d, %d, %d, %d\n", last_index_l_DBC13, last_index_l_DAD12, last_index_l_DBD13, last_index_l_DAC12, last_index_l_DBC23, last_index_l_DAD21);
                    //printf("last_index_l_DBD23,last_index_l_DAC21,last_index_l_DBC31,last_index_l_DAD32,last_index_l_DBD31,last_index_l_DAC32 = = %d, %d, %d, %d, %d, %d\n", last_index_l_DBD23,last_index_l_DAC21,last_index_l_DBC31,last_index_l_DAD32,last_index_l_DBD31,last_index_l_DAC32);
                    //printf("l = %d\n", l);
                    //double l_squared = l*l;
                    /*index_of_l(l,
                               &index_l_found,
                               &last_index_l2,
                               ppt,
                               ptr);*/
                    //printf("l = %d closest = %d\n", l, ptr->l[index_l_found]);
                    double sum_first = 0.0;
                    //for (int index_tau_bessel_first = 0; index_tau_bessel_first < bessel_boost*(index_tau_first+1)-1; index_tau_bessel++) {

                      interpolate_Dl_for_l(l2,
                                           &type_A,
                                           &type_C,
                                           &bin1,
                                           &bin2,
                                           &last_index_l_DAC12,
                                           &index_tau_first,
                                           &index_tau_second,
                                           &DAC12,
                                           pgb2,
                                           ptr,
                                           ppt);


                    //  printf("DAC12 = %g \n", DAC12);
                      //printf("last_index_l_DBD13 = %d\n", last_index_l_DBD13);

                       interpolate_Dl_for_l(l3,
                                            &type_B,
                                            &type_D,
                                            &bin1,
                                            &bin3,
                                            &last_index_l_DBD13,
                                            &index_tau_first,
                                            &index_tau_third,
                                            &DBD13,
                                            pgb2,
                                            ptr,
                                            ppt);

                        //printf("DBD13 = %g \n", DBD13);


                        interpolate_Dl_for_l(l2,
                                             &type_B,
                                             &type_C,
                                             &bin1,
                                             &bin2,
                                             &last_index_l_DBC12,
                                             &index_tau_first,
                                             &index_tau_second,
                                             &DBC12,
                                             pgb2,
                                             ptr,
                                             ppt);

                       //printf("last_index_l_DBC13 = %d\n", last_index_l_DBC13);
                       interpolate_Dl_for_l(l3,
                                            &type_A,
                                            &type_D,
                                            &bin1,
                                            &bin3,
                                            &last_index_l_DAD13,
                                            &index_tau_first,
                                            &index_tau_third,
                                            &DAD13,
                                            pgb2,
                                            ptr,
                                            ppt);

                      //  sum_first += (DAC12*DBD13+DBC12*DAD13)*w_trapz_lens_bessel2[index_tau][index_tau_bessel];
                      //}

                      interpolate_Dl_for_l(l1,
                                           &type_A,
                                           &type_C,
                                           &bin2,
                                           &bin1,
                                           &last_index_l_DAC21,
                                           &index_tau_second,
                                           &index_tau_first,
                                           &DAC21,
                                           pgb2,
                                           ptr,
                                           ppt);

                       interpolate_Dl_for_l(l3,
                                            &type_B,
                                            &type_D,
                                            &bin2,
                                            &bin3,
                                            &last_index_l_DBD23,
                                            &index_tau_second,
                                            &index_tau_third,
                                            &DBD23,
                                            pgb2,
                                            ptr,
                                            ppt);


                       interpolate_Dl_for_l(l1,
                                            &type_B,
                                            &type_C,
                                            &bin2,
                                            &bin1,
                                            &last_index_l_DBC21,
                                            &index_tau_second,
                                            &index_tau_first,
                                            &DBC21,
                                            pgb2,
                                            ptr,
                                            ppt);

                       interpolate_Dl_for_l(l3,
                                            &type_A,
                                            &type_D,
                                            &bin2,
                                            &bin3,
                                            &last_index_l_DAD23,
                                            &index_tau_second,
                                            &index_tau_third,
                                            &DAD23,
                                            pgb2,
                                            ptr,
                                            ppt);


                      interpolate_Dl_for_l(l2,
                                           &type_A,
                                           &type_C,
                                           &bin3,
                                           &bin2,
                                           &last_index_l_DAC32,
                                           &index_tau_third,
                                           &index_tau_second,
                                           &DAC32,
                                           pgb2,
                                           ptr,
                                           ppt);

                       interpolate_Dl_for_l(l1,
                                            &type_B,
                                            &type_D,
                                            &bin3,
                                            &bin1,
                                            &last_index_l_DBD31,
                                            &index_tau_third,
                                            &index_tau_first,
                                            &DBD31,
                                            pgb2,
                                            ptr,
                                            ppt);


                        interpolate_Dl_for_l(l2,
                                             &type_B,
                                             &type_C,
                                             &bin3,
                                             &bin2,
                                             &last_index_l_DBC32,
                                             &index_tau_third,
                                             &index_tau_second,
                                             &DBC32,
                                             pgb2,
                                             ptr,
                                             ppt);

                       interpolate_Dl_for_l(l1,
                                            &type_A,
                                            &type_D,
                                            &bin3,
                                            &bin1,
                                            &last_index_l_DAD31,
                                            &index_tau_third,
                                            &index_tau_first,
                                            &DAD31,
                                            pgb2,
                                            ptr,
                                            ppt);




                        double minI, minII, minIII;
                        double maxI, maxII, maxIII;
                        double AI, AII, AIII;
                        double min2I, min2II, min2III;
                        double max2I, max2II, max2III;
                        double A2I, A2II, A2III;
                        double min3I, min3II, min3III;
                        double max3I, max3II, max3III;
                        double A3I, A3II, A3III;

                        class_call(drc3jj (l2,
                                           l3,
                                           1,
                                           -1,
                                           &minI,
                                           &maxI,
                                           &AI,
                                           1000,
                                           pgb2->error_message),
                                           pgb2->error_message,
                                           pgb2->error_message);

                        class_call(drc3jj (l2,
                                           l3,
                                           -1,
                                           1,
                                           &minII,
                                           &maxII,
                                           &AII,
                                           1000,
                                           pgb2->error_message),
                                           pgb2->error_message,
                                           pgb2->error_message);

                        class_call(drc3jj (l2,
                                           l3,
                                           0,
                                           0,
                                           &minIII,
                                           &maxIII,
                                           &AIII,
                                           1000,
                                           pgb2->error_message),
                                           pgb2->error_message,
                                           pgb2->error_message);


                        double A = 0.5*(AI+AII)/AIII;

                        class_call(drc3jj (l1,
                                           l3,
                                           1,
                                           -1,
                                           &min2I,
                                           &max2I,
                                           &A2I,
                                           1000,
                                           pgb2->error_message),
                                           pgb2->error_message,
                                           pgb2->error_message);

                        class_call(drc3jj (l1,
                                           l3,
                                           -1,
                                           1,
                                           &min2II,
                                           &max2II,
                                           &A2II,
                                           1000,
                                           pgb2->error_message),
                                           pgb2->error_message,
                                           pgb2->error_message);

                        class_call(drc3jj (l1,
                                           l3,
                                           0,
                                           0,
                                           &min2III,
                                           &max2III,
                                           &A2III,
                                           1000,
                                           pgb2->error_message),
                                           pgb2->error_message,
                                           pgb2->error_message);

                        double A2 = 0.5*(A2I+A2II)/A2III;

                        class_call(drc3jj (l1,
                                           l2,
                                           1,
                                           -1,
                                           &min3I,
                                           &max3I,
                                           &A3I,
                                           1000,
                                           pgb2->error_message),
                                           pgb2->error_message,
                                           pgb2->error_message);

                        class_call(drc3jj (l1,
                                           l2,
                                           -1,
                                           1,
                                           &min3II,
                                           &max3II,
                                           &A3II,
                                           1000,
                                           pgb2->error_message),
                                           pgb2->error_message,
                                           pgb2->error_message);

                        class_call(drc3jj (l1,
                                           l2,
                                           0,
                                           0,
                                           &min3III,
                                           &max3III,
                                           &A3III,
                                           1000,
                                           pgb2->error_message),
                                           pgb2->error_message,
                                           pgb2->error_message);

                        double A3 = 0.5*(A3I+A3II)/A3III;



                    //printf("%g    %g\n", z3, A*(DAC12*DBD13+DBC12*DAD13+DAC21*DBD23+DBC21*DAD23+DAC32*DBD31+DBC32*DAD31));
                    printf("%g    %g\n", z3, A*sqrt(l2*(l2+1))*sqrt(l3*(l3+1))*(DAC12*DBD13+DBC12*DAD13)
                                             +A2*sqrt(l1*(l1+1))*sqrt(l3*(l3+1))*(DAC21*DBD23+DBC21*DAD23)
                                             +A3*sqrt(l2*(l2+1))*sqrt(l1*(l1+1))*(DAC32*DBD31+DBC32*DAD31));
                    /*printf("DAC12 = %g\n", DAC12);
                    printf("DBD13 = %g\n", DBD13 );
                    printf("DBC12 = %g\n", DBC12 );
                    printf("DAD13 = %g\n", DAD13 );
                    printf("DAC21 = %g\n", DAC21 );
                    printf("DBD23 = %g\n", DBD23 );
                    printf("DBC21 = %g\n", DBC21 );
                    printf("DAD23 = %g\n", DAD23 );
                    printf("DAC32 = %g\n", DAC32 );
                    printf("DBD31 = %g\n", DBD31 );
                    printf("DBC32 = %g\n", DBC32 );
                    printf("DAD31 = %g\n", DAD31 );*/

                    //printf("Dl_plus = %g\n",pgb2->Dl[0][0][index_l_found][bin1][bin2][index_tau_first][index_tau_second] );
                }
              }
            }
          }
        }
      }
    }

  printf("End of galbispectra2!\n");
  return _SUCCESS_;

} // end of galbispectra2_init()
