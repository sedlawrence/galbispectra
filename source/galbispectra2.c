/* Module to compute galaxy bispectra */

#include <math.h>
#include "galbispectra2.h"
#include "perturbations2.h"

/*==================================================================
====================================================================
=============================   TOOLS ==============================
====================================================================
===================================================================*/





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
  printf("inside iol: last_index_l = %d\n", *last_index_l);
  if (*last_index_l > 0){

    l_start = ptr->l[*last_index_l - 1];
    l_end = ptr->l[*last_index_l];

    if(l_end == l){
      *index_l=*last_index_l;
      printf("1. setting *index_l = *last_index_l\n" );
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
        printf("1. setting *index_l = *index\n" );
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
       printf("2. setting *index_l = *index\n" );
       *index_l=index;
       *last_index_l=index;
       break;
     }
    }
  }

  if (l == ptr->l[*last_index_l]) {
    printf("setting *index_l = *last_index_l\n" );
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

  if(k<k_start){printf("inside index_of_k: THING THAT SHOULD NOT HAPPEN HAPPENED\n");exit(2);}
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
  int bessel_boost = 2;
  double e = 2.71828;
  tau0 = pba->conformal_age;

  printf("tau0 = %g\n",tau0);


  /* Please ensure the pgb2->tau_size_selection grid is of a higher resolution than pgb2->r_size * pgb2->alpha_size = 13*/
  pgb2->tau_size_cls = 500; //prev on 500
  pgb2->tau_size_selection = 1; //prev on 601
  pgb2->r_size = 150; //previously on 1501 1515

  /* alpha_size must be an odd positive-integer lyin order to correctly fill the pgb2->r and pgb2->alpha grids using symmetries */
   //previously on 65
  double res_factor =0; //prev 4.25

  pgb2->alpha_log_size = 1; //prev on 65
  //pgb2->alpha_log_size = ceil(pow(2,res_factor)*(pgb2->alpha_log_size-1)+1);
  /* a total of pgb2->alpha_size alpha rays are distributed across the log_alpha ray grid, which always samples the diagonal of tau1-tau2
  space and alpha_over_window, which samples the region of tau1-tau2 where the pair of window functions (for each bin pair) is sampled.
  log_to_window_ratio defines how to distrubute alpha_size across the two grids. */
  //double log_alpha_fraction = 0.5;
  //pgb2->alpha_size = ceil(pgb2->alpha_log_size/log_alpha_fraction);
  //pgb2->alpha_window_size = pgb2->alpha_size-pgb2->alpha_log_size;
  if (pgb2->alpha_log_size % 2 == 0) {
    pgb2->alpha_log_size = pgb2->alpha_log_size + 1;
    printf("ERROR! alpha_size must be an ODD positive-integer\n" );
    exit(0);
  }
  pgb2->alpha_window_size = 65; //1025
  //pgb2->alpha_window_size = ceil(pow(2,res_factor)*(pgb2->alpha_window_size-1)+1);
  pgb2->alpha_size = pgb2->alpha_log_size+pgb2->alpha_window_size;

  /* The epsilon parameter defines how the logarithmically sampled alpha points are sampled. The lower this number, the more
    conentrated they are about PI/4. */
  double epsilon = 1e-5;

  int N_step = (pgb2->alpha_log_size-1)/2-1;
  printf("N_step = %d\n", N_step);
  printf("pgb2->alpha_log_size = %d\n", pgb2->alpha_log_size);


  double alpha_ratio = pow(_PI_/4/epsilon,1./N_step);
  printf("**epsilon*(ratio)^(N_step) = %g**\n", epsilon*(pow(alpha_ratio,N_step)));
  printf("alpha_ratio = %g\n", alpha_ratio);



  /* make the following 1 if alpha is to be sampled logarithmically */
  int alpha_log_sampling = 1;
  printf("alpha_size = %d\n",pgb2->alpha_size);
  printf("alpha_log_size = %d\n",pgb2->alpha_log_size);
  printf("alpha_window_size = %d", pgb2->alpha_window_size);



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

  clock_t begin = clock();








/************************************************************
 ********* Second-order matter kernel computation *********
************************************************************/
/* Uncomment the following lines to print the second-order CDM kernel */
  /* Fix k2 = 10e-5, find what this value is in both the ppt and ppt2 grids fpectively */
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



  class_alloc(pgb2->r,
              ppt->selection_num * sizeof(double**),
              pgb2->error_message);

  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    class_alloc(pgb2->r[bin1],
                ppt->selection_num * sizeof(double*),
                pgb2->error_message);

    for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
      class_alloc(pgb2->r[bin1][bin2],
                  pgb2->r_size * sizeof(double),
                  pgb2->error_message);
    }
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

    for (index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
      /* Warning Dirac-case */
      pgb2->tau_sampling_selection[bin][index_tau] = selection_mean;
      //pgb2->tau_sampling_selection[bin][index_tau] = tau_min + index_tau*(tau_max-tau_min)/(pgb2->tau_size_selection-1);

    }
  }


  // ALERT HARDCODED!!
  for (index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {
    /* Warning Dirac-case */

    pgb2->tau_sampling_cls[index_tau] = overall_tau_min + index_tau*(overall_tau_max-overall_tau_min)/(pgb2->tau_size_cls-1);
    //pgb2->tau_sampling_cls[index_tau] = selection_mean;
    //printf("pgb2->tau_sampling_cls[%d] = %g\n" , index_tau, pgb2->tau_sampling_cls[index_tau]);

  }



  /* Find the ranges of alpha to span such that we sample within tau_sampling_cls[0] to tau_sampling_cls[pgb2->tau_size_selection-1]*/
  double tau_max_cls = pgb2->tau_sampling_cls[pgb2->tau_size_cls-1];
  double tau_min_cls = pgb2->tau_sampling_cls[0];
  double r_max;
  double r_min = sqrt((1.-(11750./pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))*(1.-(11750./pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))+(1.-(11750./pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))*(1.-(11750./pgb2->tau_sampling_cls[pgb2->tau_size_cls-1])));
  double r_min2 = 0.2411267469;
  printf("1. r_min = %g\n", r_min);
  printf("r_min2 = %g\n", r_min2);





  /* Under new parameterisation tau1 = tau0(1-rsin(alpha)). tau2 = tau0(1-rcos(alpha)). Max r value can be found by
  rearranging tau_min = tau1(r_max,pi/2)=tau2(r_max,0). This is assuming alpha runs over 0 to pi/2. */
  // NOTE this has been edited to be hardcoded



  printf("r spans (%g,%g) with %d points\n",r_min, r_max, pgb2->r_size);


  printf("tau_sampling_cls has %d points that span (%g, %g)\n", pgb2->tau_size_cls, pgb2->tau_sampling_cls[0], pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]);
  //printf("tau_sampling_bessel has %d points that span (%g, %g)\n", pgb2->tau_size_bessel, pgb2->tau_sampling_bessel[index_tau_bessel], pgb2->tau_sampling_cls[pgb2->tau_size_bessel-1]);

  printf("conformal_age = %g\n",pba->conformal_age);

  // Time-grid values for tau1 = tau0(1-r*sin(alpha)), tau2 = tau0(1-r*cos(alpha)). The motivation for this parameterisation is to
  // effectively sample highly oscillatory regions of upcoming integrands. We sample one dimension with fine precision and the other
  // with coarse precision. This way when we integrate over alpha and r we have fewer iterations since alpha_size * r_size < tau_size_selection^2.



  // NOTE: THIS HAS BEEN CHANGED
  //double alpha_min = 0.6071841396;
  //double  = 0.0;

  //double alpha_max = _PI_/2.0-alpha_min;
  //double alpha_max = _PI_/2.0;

  int middle_index_alpha = (pgb2->alpha_size-1)/2;
  printf("middle_index_alpha = %d\n", middle_index_alpha );

  //printf(" =%g,  %g\n", alpha_min, (tau0-tau_max_cls)/(r_min2*tau_max_cls));
  //printf("alpha_max =%g\n", alpha_max );



  /*for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    // Linear spacing
    pgb2->alpha[index_alpha] =  alpha_min+index_alpha*(alpha_max-alpha_min)/(pgb2->alpha_size-1);


  }*/

  printf("alpha_size = %d\n",pgb2->alpha_size);
  if (res_factor<0) {
    printf("ERROR! logarithmic resolution factor < 1!\n");
    exit(2);
  }


  double alpha_over_window_base = 1.34;  //prev on 1.34 2.52
  double base = pow(alpha_over_window_base, 1/pow(2,res_factor)); // previously on 1.34


  //double alpha_middle = alpha_min + (alpha_max-alpha_min)/2.;

  //for (int index_alpha = middle_index_alpha-1; index_alpha > -1; index_alpha--) {

    //pgb2->alpha[index_alpha] =  alpha_middle-pow(base,middle_index_alpha-index_alpha)*(alpha_middle-alpha_min)/(pow(base,middle_index_alpha));


  //}
  //pgb2->alpha[middle_index_alpha] = alpha_middle;
  // log spacing in z1=z2=1 case




  printf("tau_max_cls = %g\n", pgb2->tau_sampling_cls[pgb2->tau_size_selection-1] );
  /* Allocate and define the pgb2->tau_sampling_bessel2[index_tau_cls][index_tau_selection] array which is used to integrate the innermost (time/chi) integral of the lensing term.
    There is bessel_boost times more resolution in this grid compared with tau_sampling cl. The bessel function is called using the bessel2 grid,
    so it is not necessary to have as high resolution in tau_cls compared with most other first order terms (as long as the bessel_boost makes the
    resolution in the tau_sampling_bessel2 grid sufficient.). pgb2->tau_sampling_cls[index_tau] is
    equivalent to pgb2->tau_sampling_bessel2[index_tau][bessel_boost*index_tau]*/

  double ** w_trapz_lens_bessel2;
  class_alloc(pgb2->tau_sampling_bessel2,
              pgb2->tau_size_cls * sizeof(double*),
              ppt->error_message);
  class_alloc(w_trapz_lens_bessel2,
              pgb2->tau_size_cls * sizeof(double*),
              ppt->error_message);

  class_alloc(pgb2->tau_sampling_bessel,
              pgb2->tau_size_bessel * sizeof(double),
              ppt->error_message);
  for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {

    class_alloc(pgb2->tau_sampling_bessel2[index_tau],
                bessel_boost*(index_tau+1) * sizeof(double),
                ppt->error_message);
    class_alloc(w_trapz_lens_bessel2[index_tau],
                bessel_boost*(index_tau+1) * sizeof(double),
                ppt->error_message);


  }
  //printf("tau_cls[0] = %g\n", pgb2->tau_sampling_cls[0]);
  for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(0+1); index_tau_bessel++) {
    //printf("bessel_boost*(index_tau+1) = %g\n", bessel_boost*(0+1));
    pgb2->tau_sampling_bessel2[0][index_tau_bessel] = pgb2->tau_sampling_cls[0];
    //printf(" pgb2->tau_sampling_bessel2[0][%d] = %g\n", index_tau_bessel, pgb2->tau_sampling_bessel2[0][index_tau_bessel]);
  }


  for (int index_tau = 1; index_tau < pgb2->tau_size_cls; index_tau++) {
    //printf("tau_cls[%d] = %g\n", index_tau, pgb2->tau_sampling_cls[index_tau] );
    for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(index_tau+1); index_tau_bessel++) {
      double no_of_wedges = bessel_boost*(index_tau+1.)-1.;

      pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel] = pgb2->tau_sampling_cls[index_tau] + index_tau_bessel*(tau0-pgb2->tau_sampling_cls[index_tau])/(no_of_wedges);
      //printf("pgb2->tau_sampling_bessel2[%d][%d] =%g\n", index_tau, index_tau_bessel, pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel]);
    }

  }


  for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {
    class_call(array_trapezoidal_weights(pgb2->tau_sampling_bessel2[index_tau],
                                         bessel_boost*(index_tau+1),
                                         w_trapz_lens_bessel2[index_tau],
                                         pgb2->error_message),
                                         pgb2->error_message,
                                         pgb2->error_message);
    /*for (int index_tau_bessel = 0; index_tau_bessel < bessel_boost*(index_tau+1); index_tau_bessel++) {
      printf("w_trapz_lens_bessel2[%d][%d] = %g\n", index_tau, index_tau_bessel, w_trapz_lens_bessel2[index_tau][index_tau_bessel] );
    }*/
  }











  for (int index_tau_bessel = 0; index_tau_bessel < pgb2->tau_size_bessel; index_tau_bessel++) {
    pgb2->tau_sampling_bessel[index_tau_bessel] = overall_tau_min + index_tau_bessel*(pba->conformal_age-overall_tau_min)/(pgb2->tau_size_bessel-1);
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

  /*int index_test_l1, index_test_l2, index_test_l3;
  int index_test_k1, index_test_k2, index_test_k3;

  /* Bessels Debug */
  /*
  index_test_l1 = 0;
  index_test_l2 = 1;
  index_test_l3 = 2;
  index_test_k1 = 0;
  index_test_k2 = 20;
  index_test_k3 = 60;
  double j_test1,j_test2,j_test3,j_test4,j_test5,j_test6,j_test7,j_test8,j_test9;
  printf("#k_test: %g    %g    %g \n",pgb2->k_bessel[index_test_k1], pgb2->k_bessel[index_test_k2], pgb2->k_bessel[index_test_k3]);
  printf("#l_test: %d   %d    %d\n", ptr->l[index_test_l1], ptr->l[index_test_l2], ptr->l[index_test_l3]);
  //printf("#tau \n");
  for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {


    double tau_test = pgb2->tau_sampling_cls[index_tau];
    double test_k1 = pgb2->k_bessel[index_test_k1];
    double test_k2 = pgb2->k_bessel[index_test_k2];
    double test_k3 = pgb2->k_bessel[index_test_k3];

    class_call(bessel_at_x_second_deriv(pgb2, pbs, test_k1*(tau0-tau_test), index_test_l1, &j_test1), pbs->error_message, pgb2->error_message);
    class_call(bessel_at_x_second_deriv(pgb2, pbs, test_k1*(tau0-tau_test), index_test_l2, &j_test2), pbs->error_message, pgb2->error_message);
    class_call(bessel_at_x_second_deriv(pgb2, pbs, test_k1*(tau0-tau_test), index_test_l3, &j_test3), pbs->error_message, pgb2->error_message);

    class_call(bessel_at_x_second_deriv(pgb2, pbs, test_k2*(tau0-tau_test), index_test_l1, &j_test4), pbs->error_message, pgb2->error_message);
    class_call(bessel_at_x_second_deriv(pgb2, pbs, test_k2*(tau0-tau_test), index_test_l2, &j_test5), pbs->error_message, pgb2->error_message);
    class_call(bessel_at_x_second_deriv(pgb2, pbs, test_k2*(tau0-tau_test), index_test_l3, &j_test6), pbs->error_message, pgb2->error_message);

    class_call(bessel_at_x_second_deriv(pgb2,pbs, test_k3*(tau0-tau_test), index_test_l1, &j_test7), pbs->error_message, pgb2->error_message);
    class_call(bessel_at_x_second_deriv(pgb2, pbs, test_k3*(tau0-tau_test), index_test_l2, &j_test8), pbs->error_message, pgb2->error_message);
    class_call(bessel_at_x_second_deriv(pgb2, pbs, test_k3*(tau0-tau_test), index_test_l3, &j_test9), pbs->error_message, pgb2->error_message);

    printf("%g   %g   %g    %g    %g    %g    %g    %g    %g    %g\n",tau_test, j_test1,j_test2,j_test3,j_test4,j_test5,j_test6,j_test7,j_test8,j_test9 );


  }*/





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




  /* Absolute maximum value r can have is when alpha =pi/4 and tau is set to tau_min. If r1 or r2 are above this value, this is because of a sing-
  ularity in the denominator */
  //double tau_one;
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

  double **** r_window;
  double *** alpha_over_window;


  class_alloc(r_window,
              ppt->selection_num * sizeof(double***),
              pgb2->error_message);

  class_alloc(alpha_over_window,
              ppt->selection_num * sizeof(double**),
              pgb2->error_message);



  for (int i = 0; i < ppt->selection_num; i++) {
    class_alloc(r_window[i],
                ppt->selection_num * sizeof(double**),
                pgb2->error_message);
    class_alloc(alpha_over_window[i],
                ppt->selection_num * sizeof(double*),
                pgb2->error_message);


    for (int j = 0; j < ppt->selection_num; j++) {
      class_alloc(r_window[i][j],
                  pgb2->alpha_window_size * sizeof(double*),
                  pgb2->error_message);
      class_alloc(alpha_over_window[i][j],
                  pgb2->alpha_window_size * sizeof(double),
                  pgb2->error_message);

      for (int k = 0; k < pgb2->alpha_window_size; k++) {

        class_alloc(r_window[i][j][k],
                    pgb2->r_size * sizeof(double),
                    pgb2->error_message);


      }
    }
  }



  double *** alpha2;
  double **** r_bins;
  double **** w_trapz_r_bins;
  double ***  w_trapz_alpha2;

  class_alloc(alpha2,
              ppt->selection_num * sizeof(double**),
              pgb2->error_message);
  class_alloc(r_bins,
              ppt->selection_num * sizeof(double***),
              pgb2->error_message);

  class_alloc(w_trapz_alpha2,
              ppt->selection_num * sizeof(double**),
              pgb2->error_message);

  class_alloc(w_trapz_r_bins,
              ppt->selection_num * sizeof(double***),
              pgb2->error_message);

  for (int i = 0; i < ppt->selection_num; i++) {
    class_alloc(alpha2[i],
                ppt->selection_num * sizeof(double*),
                pgb2->error_message);
    class_alloc(r_bins[i],
                ppt->selection_num * sizeof(double**),
                pgb2->error_message);

    class_alloc(w_trapz_alpha2[i],
                ppt->selection_num * sizeof(double*),
                pgb2->error_message);
    class_alloc(w_trapz_r_bins[i],
                ppt->selection_num * sizeof(double**),
                pgb2->error_message);

    for (int j = 0; j < ppt->selection_num; j++) {

      class_alloc(alpha2[i][j],
                  pgb2->alpha_size * sizeof(double),
                  pgb2->error_message);

      class_alloc(r_bins[i][j],
                  pgb2->alpha_size * sizeof(double*),
                  pgb2->error_message);

      class_alloc(w_trapz_alpha2[i][j],
                  pgb2->alpha_size * sizeof(double),
                  pgb2->error_message);

      class_alloc(w_trapz_r_bins[i][j],
                  pgb2->alpha_size * sizeof(double*),
                  pgb2->error_message);

      for (int k = 0; k < pgb2->alpha_size; k++) {
        class_alloc(r_bins[i][j][k],
                    pgb2->r_size * sizeof(double),
                    pgb2->error_message);
        class_alloc(w_trapz_r_bins[i][j][k],
                    pgb2->r_size * sizeof(double),
                    pgb2->error_message);
      }
    }
  }
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
  pgb2->index_type_vp = -1; /*p denotes prime. The quad terms are linear but only appear at second order with another quad source */
  pgb2->index_type_quad_v_ppp = -1;

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
    pgb2->index_type_delta = index_type; /* Eq. 3.28 in [1812.09297] */
    index_type++;
    //pgb2->index_source_phi = index_source;
    //index_source++;
    //pgb2->index_source_v = index_type;
    //index_source++;
    //pgb2->index_source_phi = index_type;
    //index_source++;*/
  }
  // turn on for rsd

 if (k5 == 5.0){
    //pgb2->index_source_v = index_source;
    //index_source++;
    //pgb2->index_type_vp = index_source;
    //index_type++;
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
    index_source++;



    //pgb2->index_type_g1 = index_type;
    //index_type++;
    //pgb2->index_type_g2 = index_type;
    //index_type++;
    //pgb2->index_type_g3 = index_type;
    //index_type++;
    //pgb2->index_type_g4 = index_type;
    //index_type++;
    pgb2->index_type_g5 = index_type;
    index_type++;
  }*/



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

  /* Define an array of values of first order transfer functions:
              pgb2->second_order_sources_eq[index_type][index_tau][index_k_bessel] */
  class_alloc(pgb2->second_order_sources_eq, pgb2->source_size * sizeof(double **), ppt->error_message);
    for (int index_type = 0; index_type < pgb2->source_size; index_type++) {
      /* Allocate memory for pgb2->second_order_sources_eq[index_type][index_tau] */
      class_alloc(pgb2->second_order_sources_eq[index_type],
                  pgb2->tau_size_cls * sizeof(double *),
                  ppt->error_message);
      /* Allocate memory for pgb2->second_order_sources_eq[index_type] */
      for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++) {
          /* Loop over type and tau. For each of them, allocate memory
           for pgb2->second_order_sources_eq[index_type][index_tau][index_k_bessel]  */
        class_alloc(pgb2->second_order_sources_eq[index_type][index_tau],
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
    double **** window_pair;
    class_alloc(window_pair,
                ppt->selection_num * sizeof(double***),
                pgb2->error_message);

    for (int i = 0; i < ppt->selection_num; i++) {

      class_alloc(window_pair[i],
                  ppt->selection_num * sizeof(double**),
                  pgb2->error_message);

      for (int j = 0; j < ppt->selection_num; j++) {

        class_alloc(window_pair[i][j],
                    pgb2->alpha_size * sizeof(double*),
                    pgb2->error_message);

        for (int k = 0; k < pgb2->alpha_size; k++) {

          class_alloc(window_pair[i][j][k],
                      pgb2->r_size * sizeof(double),
                      pgb2->error_message);
        }
      }
    }
    printf("window_pair array allocated.\n" );
    printf("First order sources_integ allocated\n" );

    //pgb2->Dl[index_type_first][index_type_second][index_l][bin_first][bin_second][index_tau_first][index_tau_second];
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
    class_alloc(pgb2->Dl3, pgb2->type_size * sizeof(double ******), pgb2->error_message);
    for (int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
      class_alloc(pgb2->Dl3[index_type_first], pgb2->type_size * sizeof(double *****), pgb2->error_message);
      for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
        class_alloc(pgb2->Dl3[index_type_first][index_type_second], ptr->l_size[ppt->index_md_scalars] * sizeof(double ****), pgb2->error_message);
        for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
          class_alloc(pgb2->Dl3[index_type_first][index_type_second][index_l], ppt->selection_num * sizeof(double ***), ppt->error_message);
          for (int bin_first = 0; bin_first < ppt->selection_num; bin_first++) {
            class_alloc(pgb2->Dl3[index_type_first][index_type_second][index_l][bin_first], ppt->selection_num * sizeof(double**), ppt->error_message);
            for  (int bin_second = 0; bin_second < ppt->selection_num; bin_second++) {
              class_alloc(pgb2->Dl3[index_type_first][index_type_second][index_l][bin_first][bin_second], pgb2->alpha_size * sizeof(double *), pgb2->error_message);
              for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
                class_alloc(pgb2->Dl3[index_type_first][index_type_second][index_l][bin_first][bin_second][index_alpha], pgb2->r_size * sizeof(double), pgb2->error_message);
              }
            }
          }
        }
      }
    }

    /*  double ***** Dl2_integrand;

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
      }*/
      printf("Allocating size %ix%ix%ix%ix%i bytes \n", pgb2->type_size, pgb2->type_size, ptr->l_size[ppt->index_md_scalars], pgb2->tau_size_cls, pgb2->tau_size_cls);

      /* Allocate array for Cl[index_type_first][index_type_second][index_l][bin_first][bin_second] */
      /*class_alloc(pgb2->Cl, pgb2->type_size * sizeof(double ****), pgb2->error_message);
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
      }*/


      int index1=0;
      int index2=0;



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



  /* Next we cut-off the window function by integrating from the tau value corresponding to the mean redshift to a lower time value, when the
  value of the integration exceeds some percentage (e.g. 49%) we take the lower_cutoff index and do the same for the higher time values to
  find the upper_cutoff. We then go on to parameterise the r,alpha polar time grid using these bounds. */

  /* The bin cut-off corresponds to the smaller and larger index in the pgb2->tau_sampling_selection[bin] grid, which incapsulates
  a chosen proportion of the integrated weight of the Gaussian window function. */


  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    double left_sum = 0.0;
    double right_sum = 0.0;
    class_call(background_tau_of_z(
                            pba,
                            ppt->selection_mean[bin1],
                            &selection_mean_tau),
                            ppt->error_message,
                            pgb2->error_message);

    printf("bin = %d, selection mean z = %g, selection_mean_tau = %g\n", bin1, ppt->selection_mean[bin1],selection_mean_tau);

    index_of_tau_sampling_selection(selection_mean_tau,
                    bin1,
                    &selection_mean_tau_index,
                    pgb2);


    for (int index_tau_left = selection_mean_tau_index; index_tau_left > -1; index_tau_left--) {

      left_sum += selection[bin1][index_tau_left]*w_trapz[bin1][index_tau_left];

      if (left_sum > 0.499){
        bin_cut_off_lower[bin1] = index_tau_left;
        lower_cutoff = index_tau_left;

        break;
      }

      if (index_tau_left == 0){
        bin_cut_off_lower[bin1] = index_tau_left;
        lower_cutoff = index_tau_left;

        break;
      }
    }
    for (int index_tau_right = selection_mean_tau_index; index_tau_right < pgb2->tau_size_selection; index_tau_right++) {

      right_sum += selection[bin1][index_tau_right]*w_trapz[bin1][index_tau_right];


      if (right_sum > 0.499){
        bin_cut_off_upper[bin1] = index_tau_right;
        upper_cutoff = index_tau_right;

        break;
      }

      if (index_tau_right == pgb2->tau_size_selection-1){
        if (pgb2->tau_sampling_selection[bin1][pgb2->tau_size_selection-1] == tau0) {
          bin_cut_off_upper[bin1] = index_tau_right-1;
          upper_cutoff = index_tau_right-1;
        }
        else{
          bin_cut_off_upper[bin1] = index_tau_right;
          upper_cutoff = index_tau_right;
        }

        break;
      }

    }

  }


  /* Create a log-alpha grid that log samples the diagonal region (largest weight) of the tau1-tau2 space. It does this independently
  of bins, multipole and type of terms. */

  double * alpha_log;

  class_alloc(alpha_log,
            pgb2->alpha_log_size * sizeof(double),
            pgb2->error_message);

  double * test_write;

  class_alloc(test_write,
            pgb2->alpha_log_size * sizeof(double),
            pgb2->error_message);

  for (int i = 0; i < pgb2->alpha_log_size; i++) {
    test_write[i] = 5.;
  }


  int log_index_middle = ceil(pgb2->alpha_log_size/2);
  printf("log_index_middle = %d\n",log_index_middle);
  printf("pgb2->alpha_log_size = %d\n",pgb2->alpha_log_size);
  int half_size = ((pgb2->alpha_log_size+1)/2)-1;
  printf("half_size = %d\n", half_size);
  // WARNING
  double alpha_log_min = 0.; //0.; // prev on 0.7766   0.488854
  double alpha_log_max = _PI_/2.; // prev on   0.7942   0.7942 1.08194

  double alpha_log_base = 1.34; //prev 1.02 1.34 2.52
  double log_b = alpha_log_base;
  alpha_log_base = pow(alpha_log_base,1/pow(2,res_factor));


  //alpha_log[0] = 0.;
  alpha_log[log_index_middle] = _PI_/4.;
  for (int index_log_alpha = log_index_middle-1; index_log_alpha > -1; index_log_alpha--) {
    alpha_log[index_log_alpha] =  _PI_/4.-pow(alpha_log_base, half_size-index_log_alpha)*((_PI_/4.)-alpha_log_min)/(pow(alpha_log_base,half_size));
    //alpha_log[index_log_alpha] = _PI_/4.-epsilon*pow(alpha_ratio, (log_index_middle-1)-index_log_alpha);

    //printf("alpha_log[%d] = %g (difference = %g)\n", index_log_alpha, alpha_log[index_log_alpha],alpha_log[index_log_alpha]-alpha_log[index_log_alpha+1] );
  }


  /* ALLOCATING HERE IS FINE */
  for (int index_log_alpha = log_index_middle+1; index_log_alpha < pgb2->alpha_log_size; index_log_alpha++) {
    alpha_log[index_log_alpha] =  _PI_/4.+pow(alpha_log_base, index_log_alpha)*(alpha_log_max-(_PI_/4.))/(pow(alpha_log_base,pgb2->alpha_log_size-1));
    //alpha_log[index_log_alpha] = _PI_/4.+epsilon*pow(alpha_ratio, index_log_alpha-(log_index_middle+1));
    //printf("alpha_ratio = %g\n", alpha_ratio);
    //printf("alpha_log[%d] = %g (alpha-pi/4 = %g)\n", index_log_alpha, alpha_log[index_log_alpha],alpha_log[index_log_alpha]-_PI_/4. );
  }
  /*for (int index_log_alpha = 0; index_log_alpha < pgb2->alpha_log_size; index_log_alpha++) {
    //printf("alpha_log[%d] = %g (difference = %g)\n", index_log_alpha, alpha_log[index_log_alpha], alpha_log[index_log_alpha]-alpha_log[index_log_alpha-1] );
    printf("%g \n", alpha_log[index_log_alpha]);
  }*/






  /* ALLOCATING HERE IS BAD */



  /* We are going to loop over the window-bins and create our polar grid. There will be a unique grid between any two pairs of window bins,
    this is to ensure that the time-space is adequately sampled within the window-regions. We first find the full-range that each window-regions function spans */
  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {

      double tau1_max = pgb2->tau_sampling_selection[bin1][pgb2->tau_size_selection-1];
      double tau2_max = pgb2->tau_sampling_selection[bin2][pgb2->tau_size_selection-1];
      double tau1_min = pgb2->tau_sampling_selection[bin1][0];
      double tau2_min = pgb2->tau_sampling_selection[bin2][0];

      printf("5. pgb2->tau_sampling_selection[%d][%d] = %g\n",bin1,upper_cutoff,pgb2->tau_sampling_selection[bin1][bin_cut_off_upper[bin1]] );
      double tau1_cutoff_upper = pgb2->tau_sampling_selection[bin1][bin_cut_off_upper[bin1]];
      double tau1_cutoff_lower = pgb2->tau_sampling_selection[bin1][0];
      double tau2_cutoff_upper = pgb2->tau_sampling_selection[bin2][bin_cut_off_upper[bin2]];
      double tau2_cutoff_lower = pgb2->tau_sampling_selection[bin2][0];


      double r_min = sqrt((tau0-tau1_cutoff_upper)*(tau0-tau1_cutoff_upper)+(tau0-tau2_cutoff_upper)*(tau0-tau2_cutoff_upper))/tau0;


      double alpha_of_r_min = atan((tau0-tau1_cutoff_upper)/(tau0-tau2_cutoff_upper));

      double alpha_of_r_max = atan((tau0-tau1_cutoff_lower)/(tau0-tau2_cutoff_lower));

      double alpha_min = atan((tau0-tau1_cutoff_upper)/(tau0-tau2_cutoff_lower));

      double alpha_max = atan((tau0-tau1_cutoff_lower)/(tau0-tau2_cutoff_upper));


      printf("alpha_min[%d][%d] = %g\n",bin1,bin2, alpha_min);
      printf("alpha_of_r_max[%d][%d] = %g\n", bin1, bin2, alpha_of_r_max);
      //printf("alpha_min_test  = %g\n",alpha_min_test);
      printf("alpha_max[%d][%d] = %g**\n", bin1,bin2, alpha_max);
      //printf("alpha_max_test  = %g\n",alpha_max_test);
      printf("2. r_min = %g\n", r_min);



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

      double r_abs_max = (1.-tau1_cutoff_lower/tau0)/sin(alpha_of_r_max);
      double tau_check = tau0*(1.-r_min*cos(alpha_of_r_min));

      /* There is a rectangular (square) box that is cut from time-space for cross(auto)-bin correlations using the tau_cutoff
      quantities. This is the region we wish to integrate which contains the weight of the window function. We want to distribute
      the alpha "sun-rays" such that there is a proportionate number of them on the longer/shorter side of the rectangle */

      // NOTE: alpha_over_window needs to be defined. Must incorporate alpha_of_r_max!
      double ratio =(tau1_cutoff_upper-tau1_cutoff_lower)/(tau2_cutoff_upper-tau2_cutoff_lower);
      double auto_tau1, auto_tau2;
      printf("ratio =%g\n", ratio);
      //printf("tau1_cutoff_lower = %g\n", tau1_cutoff_lower);
      //printf("tau1_cutoff_upper = %g\n", tau1_cutoff_upper);

      /* Here we define the "sun-rays" that span the rectangular box in time-space, the first indices span the horizontal side of the
      rectangle up to index_of_r_max which is the index associated with the ray that goes to the bottom left corner of the rectangle. This
      is identically the longest sun ray and also the ray which separates the rays spanning the horizontal and vertical sides.*/
      if (ratio >= 1.0){
        printf("entered first if\n");
        printf("ratio >= 1\n");
        int horizontal_size = ceil((pgb2->alpha_window_size-1)*(1./(ratio+1.0)));
        printf("horizontal_size = %d\n", horizontal_size );
        //for (int index_alpha = 0; index_alpha < pgb2->alpha_window_size; index_alpha++) {
        //  alpha_over_window[bin1][bin2][index_alpha] =  alpha_min+index_alpha*(alpha_max-alpha_min)/(pgb2->alpha_window_size-1);
        //}
        for (int index_alpha = horizontal_size-1; index_alpha > -1; index_alpha--) {

          double increment_h = (alpha_of_r_max-alpha_min)/(horizontal_size-1);
          /* alpha log-spacing */
          if (alpha_log_sampling == 1) {
            //warning

            alpha_over_window[bin1][bin2][index_alpha] =  alpha_of_r_max-pow(base,horizontal_size-index_alpha)*(alpha_of_r_max-alpha_min)/(pow(base,horizontal_size));
            //alpha_over_window[bin1][bin2][index_alpha] =  alpha_of_r_max-epsilon*pow(alpha_log_base, horizontal_size-index_alpha)*((_PI_/4.)-alpha_min)/(pow(alpha_log_base,horizontal_size));
            //printf("alpha_over_window[%d][%d][%d] = %g\n", bin1, bin2, index_alpha, alpha_over_window[bin1][bin2][index_alpha] );

            //alpha_over_window[bin1][bin2][index_alpha] = _alpha_of_r_max-epsilon*pow(alpha_ratio, horizontal_size-index_alpha);
            //
          }

          if (alpha_log_sampling != 1) {
            alpha_over_window[bin1][bin2][index_alpha] = alpha_of_r_max-(horizontal_size-index_alpha)*(alpha_of_r_max-alpha_min)/(pgb2->alpha_window_size-1);

            //alpha_over_window[bin1][bin2][index_alpha] = 0.7766-(horizontal_size-index_alpha)*(0.7766-0.0)/(horizontal_size);
            //printf("alpha_over_window[%d][%d][%d] = %g\n", bin1, bin2, index_alpha, alpha_over_window[bin1][bin2][index_alpha]);
          }

          double r_max_per_alpha = (1.0-tau2_cutoff_lower/tau0)/(cos(alpha_over_window[bin1][bin2][index_alpha]));

          for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            r_window[bin1][bin2][index_alpha][index_r] = r_min+index_r*(r_max_per_alpha-r_min)/(pgb2->r_size-1);

            auto_tau1 = tau0*(1.-r_window[bin1][bin2][index_alpha][index_r]*sin(alpha_over_window[bin1][bin2][index_alpha]));
            auto_tau2 = tau0*(1.-r_window[bin1][bin2][index_alpha][index_r]*cos(alpha_over_window[bin1][bin2][index_alpha]));
          }
        }
        /* We want a sun-ray to go directly to (tau1,tau2)=(tau1_cutoff_lower,tau2_cutoff_lower) (along the r_max line for this pair of
      bins. Note that pgb2->alpha_window_size = horizontal+1+vertical. The 1 referes to this line. r_abs_max refers to the largest r value along this
      given sunray, alpha_of_r_max is the corresponding angle (alpha).*/
      // WARNING consider the following:
        alpha_over_window[bin1][bin2][horizontal_size] = alpha_of_r_max;
        //alpha_over_window[bin1][bin2][horizontal_size] = 0;

        for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
          r_window[bin1][bin2][horizontal_size][index_r] = r_min+index_r*(r_abs_max-r_min)/(pgb2->r_size-1);
        }

        int vertical_size = pgb2->alpha_window_size-1-horizontal_size;




        for (int index_alpha = horizontal_size+1; index_alpha < pgb2->alpha_window_size; index_alpha++) {
          double increment_v = (alpha_max-alpha_min)/(vertical_size-1);
          //alpha_over_window[bin1][bin2][index_alpha] = alpha_of_r_max + increment_v+index_alpha*(alpha_max-(alpha_of_r_max+increment_v))/(vertical_size-1);
          if (alpha_log_sampling == 1) {
            alpha_over_window[bin1][bin2][index_alpha] =  alpha_of_r_max+pow(base,index_alpha)*(alpha_max-alpha_of_r_max)/(pow(base,pgb2->alpha_window_size-1));
            //alpha_over_window[bin1][bin2][index_alpha] = alpha_of_r_max+epsilon*pow(alpha_ratio, index_alpha-(log_index_middle+1));
          }

          if (alpha_log_sampling != 1) {
            alpha_over_window[bin1][bin2][index_alpha] = alpha_of_r_max+(index_alpha-horizontal_size)*(alpha_max-alpha_of_r_max)/(pgb2->alpha_window_size-1);

            //alpha_over_window[bin1][bin2][index_alpha] = 0.7942+(index_alpha-horizontal_size)*(_PI_/2.-0.7942)/(pgb2->alpha_window_size-1-horizontal_size);
            //printf("alpha_over_window[%d][%d][%d] = %g\n", bin1, bin2, index_alpha, alpha_over_window[bin1][bin2][index_alpha]);
          }
          double r_max_per_alpha = (1.0-tau1_cutoff_lower/tau0)/(sin(alpha_over_window[bin1][bin2][index_alpha]));

          //printf("r_max_per_alpha[%d][%d][%d] = %g\n",bin1,bin2,index_alpha, r_max_per_alpha);
          for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            r_window[bin1][bin2][index_alpha][index_r] = r_min+index_r*(r_max_per_alpha-r_min)/(pgb2->r_size-1);
            //printf("r_window = %g\n", r_window[bin1][bin2][index_alpha][index_r]);
            auto_tau1 = tau0*(1.-r_window[bin1][bin2][index_alpha][index_r]*sin(alpha_over_window[bin1][bin2][index_alpha]));
            auto_tau2 = tau0*(1.-r_window[bin1][bin2][index_alpha][index_r]*cos(alpha_over_window[bin1][bin2][index_alpha]));
          }
        }
      }


      //NOTE  ALERT this needs to be changed to resemble the above
      if (ratio < 1.0 && ratio >0.) {
        printf("entered second if\n");
        printf("ratio should be less than 1 but greater than zero *%g*\n", ratio );
        //int horizontal_size = ceil((pgb2->alpha_window_size-1)*(1.0-(ratio/(ratio+1.0))));

        ratio =1./ratio;
        printf("ratio = %d\n", ratio);
        int vertical_size = ceil((pgb2->alpha_window_size-1)*(1./(ratio+1.0)));
        int horizontal_size = pgb2->alpha_window_size-vertical_size-1;

        for (int index_alpha = horizontal_size-1; index_alpha > -1; index_alpha--) {

          double increment_h = (alpha_of_r_max-alpha_min)/(horizontal_size-1);
          if (alpha_log_sampling != 1) {
            alpha_over_window[bin1][bin2][index_alpha] = alpha_of_r_max-(horizontal_size-index_alpha)*(alpha_of_r_max-alpha_min)/(pgb2->alpha_window_size-1);
          }
          if (alpha_log_sampling == 1) {
            alpha_over_window[bin1][bin2][index_alpha] =  alpha_of_r_max-pow(base,horizontal_size-index_alpha)*(alpha_of_r_max-alpha_min)/(pow(base,horizontal_size));
          }
          double r_max_per_alpha = (1.0-tau2_cutoff_lower/tau0)/(cos(alpha_over_window[bin1][bin2][index_alpha]));

          for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            r_window[bin1][bin2][index_alpha][index_r] = r_min+index_r*(r_max_per_alpha-r_min)/(pgb2->r_size-1);

            auto_tau1 = tau0*(1.-r_window[bin1][bin2][index_alpha][index_r]*sin(alpha_over_window[bin1][bin2][index_alpha]));
            auto_tau2 = tau0*(1.-r_window[bin1][bin2][index_alpha][index_r]*cos(alpha_over_window[bin1][bin2][index_alpha]));
            }
          }
        /* We want a sun-ray to go directly to (tau1,tau2)=(tau1_cutoff_lower,tau2_cutoff_lower) (along the r_max line for this pair of
      bins. Note that pgb2->alpha_window_size = horizontal+1+vertical. The 1 referes to this line. r_abs_max refers to the largest r value along this
      given sunray, alpha_of_r_max is the corresponding angle (alpha).*/

        alpha_over_window[bin1][bin2][horizontal_size] = alpha_of_r_max;
        for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
          r_window[bin1][bin2][horizontal_size][index_r] = r_min+index_r*(r_abs_max-r_min)/(pgb2->r_size-1);
        }

        //int vertical_size = pgb2->alpha_window_size-1-horizontal_size;

        printf("vertical_size =%d\n", vertical_size );

        for (int index_alpha = horizontal_size+1; index_alpha < pgb2->alpha_window_size; index_alpha++) {
          double increment_v = (alpha_max-alpha_min)/(vertical_size-1);
          //alpha_over_window[bin1][bin2][index_alpha] = alpha_of_r_max + increment_v+index_alpha*(alpha_max-(alpha_of_r_max+increment_v))/(vertical_size-1);
          if (alpha_log_sampling == 1) {
            alpha_over_window[bin1][bin2][index_alpha] =  alpha_of_r_max+pow(base,index_alpha)*(alpha_max-alpha_of_r_max)/(pow(base,pgb2->alpha_window_size-1));
          }

          if (alpha_log_sampling != 1) {
            alpha_over_window[bin1][bin2][index_alpha] = alpha_of_r_max+(index_alpha-horizontal_size)*(alpha_max-alpha_of_r_max)/(pgb2->alpha_window_size-1);
          }
          double r_max_per_alpha = (1.0-tau1_cutoff_lower/tau0)/(sin(alpha_over_window[bin1][bin2][index_alpha]));

          //printf("r_max_per_alpha[%d][%d][%d] = %g\n",bin1,bin2,index_alpha, r_max_per_alpha);
          for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            r_window[bin1][bin2][index_alpha][index_r] = r_min+index_r*(r_max_per_alpha-r_min)/(pgb2->r_size-1);
            //printf("r_window = %g\n", r_window[bin1][bin2][index_alpha][index_r]);
            auto_tau1 = tau0*(1.-r_window[bin1][bin2][index_alpha][index_r]*sin(alpha_over_window[bin1][bin2][index_alpha]));
            auto_tau2 = tau0*(1.-r_window[bin1][bin2][index_alpha][index_r]*cos(alpha_over_window[bin1][bin2][index_alpha]));
          }
        }
      }

      if (ratio == 0.0 || ratio < 0.0){
        printf("ERROR ratio of horizontal to vertical side of rectangle cut-out is %g. bin = %d with bin = %d\n", ratio,bin1, bin2);
        exit(2);
      }
    }
  }

  /*for (int index_alpha = 0; index_alpha < pgb2->alpha_; index_alpha++) {
    printf("alpha_over_window[0][0][%d] = %g\n", index_alpha, alpha_over_window[0][0][index_alpha]);
  }*/


  // NOTE: the else part doesn't have the increment, the middle alpha index is not correct!



  printf("tau0 = %g\n",tau0 );
  printf("tau_max_cls = %g\n", tau_max_cls);
  printf("tau_min_cls = %g\n", tau_min_cls);



  r_abs_max = (1.0-(tau_min_cls/tau0))/sin(_PI_/4);

  printf("r_abs_max = %G\n", r_abs_max);


  printf("r_abs_max = %g\n", r_abs_max);
  printf("r_abs_max = %g\n", (1-pgb2->tau_sampling_cls[0]/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1])/cos(_PI_/4));

  /*double * r_max_per_alpha;

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
  }*/



  printf("pgb2->alpha_size b4 = %d\n", pgb2->alpha_size);


  printf("pgb2->alpha_size after = %d\n", pgb2->alpha_size);
  printf("selection_num = %d\n", ppt->selection_num );

  double r_diag_min;
  double r_diag_max;

  printf("r_diag_min = %g\n",r_diag_min );
  printf("r_diag_max = %g\n",r_diag_max );

  /* We now need to merge the alpha_over_window grid with alpha_log to form alpha2-grid, that has increasing alpha values with index (so
  that trapezoidal weights are configured correctly). */

  printf("Writing the alpha2- and r_bins-grid..\n");
  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
      for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
        alpha2[bin1][bin2][index_alpha] =0.0;
      }
    }
  }
  /* We are going to form an ordered alpha2 list from the log-sampled diagonal "alpha_log" array and the window specific "alpha_over_window"
    array. We do this by looping over indices of alpha2 and iteratively comparing the two smaller arrays to find the next largest alpha-
    value to add to the final grid. Once all elements in one list are used up, we simply add all of the elements in the remaining list. */
  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
      int index_A = 0;
      int index_B = 0;
      double value_A;
      double value_B;

      /*WARNING*/
      if (bin1 == bin2) {
        //r_diag_max = (1.-12000/*pgb2->tau_sampling_selection[bin1][0]*//tau0)/sin(_PI_/4.);
        r_diag_max = (1.-pgb2->tau_sampling_selection[bin1][0]/tau0)/sin(_PI_/4.);
        //r_diag_max = 0.15;
      }
      // compare redshifts
      double bin1_mean = ppt->selection_mean[bin1];
      double bin2_mean = ppt->selection_mean[bin2];

      if (bin1_mean < bin2_mean) {
        r_diag_max = (1.-pgb2->tau_sampling_selection[bin1][0]/tau0)/sin(_PI_/4.);
      }

      if (bin1_mean > bin2_mean) {
        r_diag_max = (1.-pgb2->tau_sampling_selection[bin2][0]/tau0)/cos(_PI_/4.);
      }

      r_diag_min =  r_window[bin1][bin2][0][0];

      double random_fraction = 0.4;

      double rand_index_multiple = ceil(1./random_fraction);
      double n =1;
      for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
        pgb2->r[bin1][bin2][index_r] = r_diag_min+(r_diag_max-r_diag_min)*(index_r)/(pgb2->r_size-1);

        /* Replace the above line with the following if random r-sampling is required */

        /*if (index_r == n*rand_index_multiple) {
          double lower = pgb2->r[bin1][bin2][index_r-1];
          double upper =  r_diag_min+(r_diag_max-r_diag_min)*(index_r)/(pgb2->r_size-1);

          /* Define a random number between 1 and 10 */
        /*  int random = rand() % (10 - 1 + 1) + 1;

          pgb2->r[bin1][bin2][index_r] = lower+(1./random)*(upper-lower);
          //printf("index_r = %d, lower = %g, upper = %g, random = %d, r =%g \n", index_r, lower, upper, random, pgb2->r[bin1][bin2][index_r]);
          n+=1;
        }
        else{
          pgb2->r[bin1][bin2][index_r] = r_diag_min+(r_diag_max-r_diag_min)*(index_r)/(pgb2->r_size-1);
        }*/

      }



      for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {

        if((index_A >= pgb2->alpha_log_size) && (index_B >= pgb2->alpha_window_size)){
          printf("breaking loop.\n" );
          break;
        }

        /* Check alpha hasn't already been written for this index, this occurs when we have the same alpha value in the two
        smaller grids. */
        double a = alpha2[bin1][bin2][index_alpha];

        if (a != 0){

          continue;
        }

        if( index_A >= pgb2->alpha_log_size){


          value_B = alpha_over_window[bin1][bin2][index_B];
          for (int i = index_B; i < pgb2->alpha_window_size; i++) {
            alpha2[bin1][bin2][index_alpha] = alpha_over_window[bin1][bin2][i];

            //for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            //  r_bins[bin1][bin2][index_alpha][index_r] = r_diag_min+(r_diag_max-r_diag_min)*(index_r)/(pgb2->r_size-1);/*r_window[bin1][bin2][i][index_r]*/;
            //  r_bins[bin1][bin2][index_alpha][index_r] = pgb2->r[bin1][bin2][index_r];
            //}
            index_alpha +=1;
          }
          break;
        }

        if( index_B >= pgb2->alpha_window_size){
          value_A = alpha_log[index_A];
          for (int i = index_A; i < pgb2->alpha_log_size; i++) {
            alpha2[bin1][bin2][index_alpha] = alpha_log[i];

            //for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            //  r_bins[bin1][bin2][index_alpha][index_r] = r_diag_min+(r_diag_max-r_diag_min)*(index_r)/(pgb2->r_size-1);
              //r_bins[bin1][bin2][index_alpha][index_r] = pgb2->r[bin1][bin2][index_r];
            //}
            index_alpha +=1;
          }
          break;
        }

        value_A = alpha_log[index_A];
        value_B = alpha_over_window[bin1][bin2][index_B];


        if (value_A < value_B) {
          alpha2[bin1][bin2][index_alpha] = value_A;

          //for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            //printf("%dx%dx%dx%d\n", bin1, bin2, index_alpha, index_r);
          //  r_bins[bin1][bin2][index_alpha][index_r] = r_diag_min+(r_diag_max-r_diag_min)*(index_r)/(pgb2->r_size-1);
            //r_bins[bin1][bin2][index_alpha][index_r] = pgb2->r[bin1][bin2][index_r];
          //}
          index_A +=1;
        }

        if (value_B < value_A) {
          alpha2[bin1][bin2][index_alpha] = value_B;

          //for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
          //  r_bins[bin1][bin2][index_alpha][index_r] = r_diag_min+(r_diag_max-r_diag_min)*(index_r)/(pgb2->r_size-1);//r_window[bin1][bin2][index_B][index_r];
            //r_bins[bin1][bin2][index_alpha][index_r] = pgb2->r[bin1][bin2][index_r];
          //}
          index_B +=1;
        }

        if (value_A == value_B) {

          alpha2[bin1][bin2][index_alpha] = value_A;

          alpha2[bin1][bin2][index_alpha+1] = value_B;

          //for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
          //  r_bins[bin1][bin2][index_alpha][index_r] = r_diag_min+(r_diag_max-r_diag_min)*(index_r)/(pgb2->r_size-1);//r_window[bin1][bin2][index_B][index_r];
            //r_bins[bin1][bin2][index_alpha][index_r] = pgb2->r[bin1][bin2][index_r];
          //}
          //for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
          //  r_bins[bin1][bin2][index_alpha+1][index_r] = r_diag_min+(r_diag_max-r_diag_min)*(index_r)/(pgb2->r_size-1);//r_window[bin1][bin2][index_B][index_r];
            //r_bins[bin1][bin2][index_alpha][index_r] = pgb2->r[bin1][bin2][index_r];
          //}
          index_A+=1;
          index_B+=1;
        }
      }
    }
  }
  double number;
  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
      for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
        for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
            number = pgb2->r[bin1][bin2][index_r];
            r_bins[bin1][bin2][index_alpha][index_r] = number;
            //printf("r_bins[0][0][%d][%d] = %g, pgb2->r = %g (step_size: %g)\n", index_alpha, index_r, r_bins[bin1][bin2][index_alpha][index_r], pgb2->r[bin1][bin2][index_r], 1000.*(r_bins[bin1][bin2][0][index_r]-r_bins[bin1][bin2][0][index_r-1]));
        }
      }
    }
  }







  /* Time debug PRINT SUN RAY time-sampling*/
  /*double tau100, tau200, tau101, tau201, tau110, tau210, tau111, tau211;
  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
      tau100 = tau0*(1.-r_bins[0][0][index_alpha][index_r]*sin(alpha2[0][0][index_alpha]));
      tau200 = tau0*(1.-r_bins[0][0][index_alpha][index_r]*cos(alpha2[0][0][index_alpha]));
      tau101 = tau0*(1.-r_bins[0][1][index_alpha][index_r]*sin(alpha2[0][1][index_alpha]));
      tau201 = tau0*(1.-r_bins[0][1][index_alpha][index_r]*cos(alpha2[0][1][index_alpha]));
      tau110 = tau0*(1.-r_bins[1][0][index_alpha][index_r]*sin(alpha2[1][0][index_alpha]));
      tau210 = tau0*(1.-r_bins[1][0][index_alpha][index_r]*cos(alpha2[1][0][index_alpha]));
      tau111 = tau0*(1.-r_bins[1][1][index_alpha][index_r]*sin(alpha2[1][1][index_alpha]));
      tau211 = tau0*(1.-r_bins[1][1][index_alpha][index_r]*cos(alpha2[1][1][index_alpha]));

      //printf("%g      %g      %g      %g      %g      %g      %g      %g\n", tau100, tau200, tau101, tau201, tau110, tau210, tau111,tau211 );
    }
  }*/

  /* this only works for z1=z2, we are interpolating to find the exact dirac result of Dl between two r values for the selection
    redshift */
  double selection_mean_tau_bin1;

  class_call(background_tau_of_z(
                        pba,
                        ppt->selection_mean[0],
                        &selection_mean_tau_bin1),
                        ppt->error_message,
                        pgb2->error_message);
  printf("1. selection_mean_tau_bin1 = %g\n", selection_mean_tau_bin1);

  //printf("bin = %d, selection mean z = %g, selection_mean_tau = %g\n", bin1, ppt->selection_mean[0], selection_mean_tau_bin1);

  index_of_tau_sampling_selection(selection_mean_tau_bin1,
                  0,
                  &selection_mean_tau_index,
                  pgb2);
  double r_find;

  r_find = (1.-selection_mean_tau_bin1/tau0)/sin(alpha2[0][0][middle_index_alpha]);

  int closest_r_index_to_mean_bin;

  double double_index_r = (r_find - r_bins[0][0][middle_index_alpha][0])/(r_bins[0][0][middle_index_alpha][pgb2->r_size-1] - r_bins[0][0][middle_index_alpha][0]) * (pgb2->r_size);
  printf("double_index_r = %g\n", double_index_r);
  int index_r_find = (int) floor(double_index_r);

  if (index_r_find < 1) {
    index_r_find = 1;
  }

  else if(index_r_find > pgb2->r_size -1){
    index_r_find = pgb2->r_size -1;
  }
  double tau_found = tau0*(1.-r_bins[0][0][middle_index_alpha][index_r_find]*alpha2[0][0][middle_index_alpha]);
  printf("tau value found is %g\n", tau_found);




  /* Initialise and fill an array window_pair, parameterised by the polar r- and alpha-parameters which pick the exact (as opposed to
  interpolated value of the product of the window function product */

  double tau_one_wind;
  double tau_two_wind;
  double * pvecback_wind1;
  double * pvecback_wind2;
  double dNdz1;
  double dln_dNdz_dz1;
  double dNdz2;
  double dln_dNdz_dz2;
  double selection1;
  double selection2;
  class_alloc(pvecback_wind2, pba->bg_size * sizeof(double), pba->error_message);
  class_alloc(pvecback_wind1, pba->bg_size * sizeof(double), pba->error_message);


  printf("Completed sampling.\n");

  /*for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    printf("alpha2[0][0][%d] = %g\n", index_alpha, alpha2[0][0][index_alpha] );
  }*/
  /*for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
    printf("r_bins[0][0][0][%d] = %g\n", index_r, r_bins[0][0][0][index_r] );
  }*/
  //exit(0);


















  /*==================================================================
  ====================================================================
  =============================   SOURCES  ===========================
  ====================================================================
  ===================================================================*/
  printf("Starting sources..\n" );

  /*double pq = pgb2->k_bessel[0];
  int testy_index = 0;
  printf("pq = %g\n",pq);
  class_call(index_of_k(pq, &testy_index, ppt), pgb2->error_message, pgb2->error_message);
  printf("exited search function.\n");
  exit(0);*/
  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
      for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
        for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
          tau_one_wind = tau0*(1.0-r_bins[bin1][bin2][index_alpha][index_r]*sin(alpha2[bin1][bin2][index_alpha]));
          tau_two_wind = tau0*(1.0-r_bins[bin1][bin2][index_alpha][index_r]*cos(alpha2[bin1][bin2][index_alpha]));

          class_call(background_at_tau(pba,
                                       tau_one_wind,
                                       pba->long_info,
                                       pba->inter_normal,
                                       &last_index,
                                       pvecback_wind1),
                     pba->error_message,
                     ptr->error_message);

          /* infer redshift */


          class_call(background_at_tau(pba,
                                       tau_two_wind,
                                       pba->long_info,
                                       pba->inter_normal,
                                       &last_index,
                                       pvecback_wind2),
                     pba->error_message,
                     ptr->error_message);

          /* infer redshift */
          double z1 = pba->a_today/pvecback_wind1[pba->index_bg_a]-1.;
          double z2 = pba->a_today/pvecback_wind2[pba->index_bg_a]-1.;



          class_call(transfer_selection_function(ppr,
                                      ppt,
                                      ptr,
                                      bin1,
                                      z1,
                                      &selection1),
                                      ptr->error_message,
                                      pgb2->error_message);

          class_call(transfer_selection_function(ppr,
                                      ppt,
                                      ptr,
                                      bin2,
                                      z2,
                                      &selection2),
                                      ptr->error_message,
                                      pgb2->error_message);
          /* get corresponding dN/dtau = dN/dz * dz/dtau = dN/dz * H */
          selection1 *= pvecback_wind1[pba->index_bg_H];
          selection2 *= pvecback_wind1[pba->index_bg_H];


          window_pair[bin1][bin2][index_alpha][index_r] = selection1*selection2; //window_bin1*pvecback_wind1[pba->index_bg_H];//*window_bin2*pvecback_wind2[pba->index_bg_H];
          //printf("window_pair[%d][%d][%d][%d] = %g\n", bin1, bin2, index_alpha, index_r, window_pair[bin1][bin2][index_alpha][index_r]);
          //printf("z1 = %g\n", z1 );
          //printf("z2 = %g\n", z2 );


        }
      }
    }
  }

  /* Create trapezoidal weights for r- and alpha- integration */
  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
      class_call(array_trapezoidal_weights(alpha2[bin1][bin2],
                                           pgb2->alpha_size,
                                           w_trapz_alpha2[bin1][bin2],
                                           pgb2->error_message),
                                           pgb2->error_message,
                                           pgb2->error_message);

      for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
        class_call(array_trapezoidal_weights(r_bins[bin1][bin2][index_alpha],
                                             pgb2->r_size,
                                             w_trapz_r_bins[bin1][bin2][index_alpha],
                                             pgb2->error_message),
                                             pgb2->error_message,
                                             pgb2->error_message);
      }
    }
  }
  /* Now we normalise the window function pair such that  tau0*tau0*r*window1*window2 dadr =1. */
  double inner_window_sum;
  double outer_window_sum;

  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
      outer_window_sum = 0.0;
      for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
        inner_window_sum = 0.0;
        for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
          inner_window_sum += tau0*tau0*r_bins[bin1][bin2][index_alpha][index_r]*window_pair[bin1][bin2][index_alpha][index_r]*w_trapz_r_bins[bin1][bin2][index_alpha][index_r];

        }

        outer_window_sum += inner_window_sum*w_trapz_alpha2[bin1][bin2][index_alpha];
      }
      printf("window_pair renorm factor = %g\n", 1./outer_window_sum);
      for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
        for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
          window_pair[bin1][bin2][index_alpha][index_r]/=outer_window_sum;
        }
      }
    }
  }

  double sum_test_out = 0.0;
  for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
    double sum_test_in = 0.0;
    for (int index_r = 0; index_r < pgb2->r_size; index_r++) {
      sum_test_in += tau0*tau0*r_bins[0][0][index_alpha][index_r]*window_pair[0][0][index_alpha][index_r]*w_trapz_r_bins[0][0][index_alpha][index_r];
    }
    sum_test_out+= sum_test_in*w_trapz_alpha2[0][0][index_alpha];
  }
  printf("sum_test_out (should be 1)= %g\n", sum_test_out);






  double * w_trapz_r;
  double ** w_trapz_r2;





  /* Because we are manipulating the window function to emphasise regions of interest, we need to renomalise it such that when
  the product of the two window functions is integrated over r and alpha, the result is one. We will do this by working out a
  renorm factor for any pair of bins and using this later in the double time (r and alpha) integration over Dl.*/

  double ** renorm_factor;

  class_alloc(renorm_factor,
              ppt->selection_num * sizeof(double*),
              pgb2->error_message);

  for (int i = 0; i < ppt->selection_num; i++) {
    class_alloc(renorm_factor[i],
                ppt->selection_num * sizeof(double),
                pgb2->error_message);
  }
  double selection_sum = 0.;
  for (int index_tau = 0; index_tau < pgb2->tau_size_selection; index_tau++) {
    selection_sum += selection[0][index_tau]*w_trapz[0][index_tau];
  }
  printf("selection_sum = %g\n",selection_sum);
  double inner_sum, outer_sum;
  double window1_renorm,window2_renorm,tau_one_renorm,tau_two_renorm;
  int index_tau_first_polar_rn, index_tau_second_polar_rn;
  for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
    for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {

      outer_sum = 0.0;
      for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
        inner_sum = 0.0;

        for (int index_r = 0; index_r < pgb2->r_size; index_r++) {

          tau_one_renorm = tau0*(1.-r_bins[bin1][bin2][index_alpha][index_r]*sin(alpha2[bin1][bin2][index_alpha]));
          tau_two_renorm = tau0*(1.-r_bins[bin1][bin2][index_alpha][index_r]*cos(alpha2[bin1][bin2][index_alpha]));
          //printf("tau1 = %g, tau2 = %g\n",tau_one_renorm, tau_two_renorm);

          index_of_tau_sampling_selection(tau_one_renorm, bin1, &index_tau_first_polar_rn, pgb2);
          index_of_tau_sampling_selection(tau_two_renorm, bin2, &index_tau_second_polar_rn, pgb2);


          /* Interpolate between the two indices found on the tau_sampling_selection grid to find the exact W(tau_one_renorm) and W(tau_two_renorm) */


          window1_renorm = selection[bin1][index_tau_first_polar_rn-1]*(pgb2->tau_sampling_selection[bin1][index_tau_first_polar_rn]-tau_one_renorm)
                        + selection[bin1][index_tau_first_polar_rn]*(tau_one_renorm-pgb2->tau_sampling_selection[bin1][index_tau_first_polar_rn-1]);
          window1_renorm /= (pgb2->tau_sampling_selection[bin1][index_tau_first_polar_rn] - pgb2->tau_sampling_selection[bin1][index_tau_first_polar_rn-1]);

          window2_renorm = selection[bin2][index_tau_second_polar_rn-1]*(pgb2->tau_sampling_selection[bin2][index_tau_second_polar_rn]-tau_two_renorm )
                        + selection[bin2][index_tau_second_polar_rn]*(tau_two_renorm -pgb2->tau_sampling_selection[bin2][index_tau_second_polar_rn-1]);
          window2_renorm /= (pgb2->tau_sampling_selection[bin2][index_tau_second_polar_rn] - pgb2->tau_sampling_selection[bin2][index_tau_second_polar_rn-1]);


          //printf("window1 = %g, window 2 = %g\n", window1_renorm, window2_renorm );
          inner_sum += tau0*tau0*r_bins[bin1][bin2][index_alpha][index_r]*window1_renorm*window2_renorm*w_trapz_r_bins[bin1][bin2][index_alpha][index_r];
          //printf("inner_sum = %g\n", inner_sum);
        }
        outer_sum += inner_sum*w_trapz_alpha2[bin1][bin2][index_alpha];
        //printf("outer_sum = %g\n", outer_sum);
      }
      printf("outer_sum = %g\n", outer_sum);
      renorm_factor[bin1][bin2] = 1./outer_sum;
      printf("renorm_factor[%d][%d] = %g\n", bin1, bin2, renorm_factor[bin1][bin2]);
    }
  }






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

            //printf("1. pgb2->first_order_sources[pgb2->index_source_psi][%d][%d] = %g\n", index_tau, index_k_bessel, pgb2->first_order_sources[pgb2->index_source_psi][index_tau][index_k_bessel] );
          }
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
          //f = 0.;
          //g = 0;

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
          //index_k_lens = 0;
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

            //class_call(index_of_k(k, &index_k_lens, ppt), pgb2->error_message, pgb2->error_message);
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

    /*for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
      for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){

        printf("pgb2->tau_sampling_cls[%d] = %g\n", index_tau, pgb2->tau_sampling_cls[index_tau]);
        printf("pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][%d][%d] = *%g*\n", index_tau, index_k_bessel, pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][index_tau][index_k_bessel] );
        printf("pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi_prime][%d][%d] = %g\n", index_tau, index_k_bessel, pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi_prime][index_tau][index_k_bessel]);

      }
    }*/




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
                double bias = 1.3;
                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);

                pgb2->first_order_sources_integ[pgb2->index_type_density][index_l][index_tau][index_k_bessel] = bias*pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_bessel] * j;
              }
            }
          }
        }



        /* RSD */
        int index_tau_bin1;
        index_of_tau_sampling_cls(selection_mean_tau_bin1, &index_tau_bin1, pgb2);
        printf("index_tau_bin1 = %d\n", index_tau_bin1 );
        printf("pgb2->tau_sampling_cls[%d] = %g\n", index_tau_bin1,pgb2->tau_sampling_cls[index_tau_bin1]);
        printf("2. selection_mean_tau_bin1 = %g\n", selection_mean_tau_bin1);

        if (pgb2->index_type_rsd != -1) {
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              printf("indexing index_l x index_tau = %dx%d\n", index_l, index_tau);
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
                                *pgb2->first_order_sources[pgb2->index_source_theta][index_tau][index_k_bessel]*
                                j_second_deriv;

              }
            }
          }
        }
        printf("completed RSD source integ.\n");


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
        double * pvecback_d2;
        class_alloc(pvecback_d2, pba->bg_size * sizeof(double), pba->error_message);
        /* First Type: Doppler2 */
        int last_index_d2;
        if(pgb2->index_type_d2 != -1){
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
        //free(pvecback_d2);
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
        double * pvecback_g2;
        class_alloc(pvecback_g2, pba->bg_size * sizeof(double), pba->error_message);
        if(pgb2->index_type_g2 != -1 ){
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
                //printf("prefactor_g2 = %g\n", prefactor_g2);

                //printf("pgb2->first_order_sources_integ[pgb2->index_type_g2][%d][%d][%d] = %g\n", index_l, index_tau, index_k_bessel, pgb2->first_order_sources_integ[pgb2->index_type_g2][index_l][index_tau][index_k_bessel] );
                //printf("2. pgb2->first_order_sources[pgb2->index_source_psi] = %g\n", pgb2->first_order_sources[pgb2->index_source_psi][index_tau][index_k_bessel] );
              }
            }
          }
        }

        /* First Type: g3 */
        double * pvecback_g3;
        class_alloc(pvecback_g3, pba->bg_size * sizeof(double), pba->error_message);
        if(pgb2->index_type_g3 != -1 ){
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
              //printf("inside lens loop, tau = %g\n", pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel]);

              class_call(index_of_tau_sampling_quadsources(tau, &index_tau_quad, ppt), pgb2->error_message, pgb2->error_message);
              index_k = 0;
              for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {





                double k = pgb2->k_bessel[index_k_bessel];
              //  printf("k = %g\n",k);
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

        /* This is for a quad-source; a term that is linear but appears specifically at second-order multiplied with one other quadsource.
        In the full galaxy number-over density at second order (\Delta^{(2)}) the term is  (-2/a/a/H/Hv')*(v'''). We separate the two
        brackets and treat them as quad_v_p (left) and quad_v_ppp (right), with the "p" denoting the derivative. */

        /* Eq. 3.32 in [1812.09297] */

        double * pvecback_vp;
        class_alloc(pvecback_vp, pba->bg_size * sizeof(double), pba->error_message);
        if(pgb2->index_type_vp != -1 ){
          int index_tau_quad;
          int last_index_vp;
          for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
            for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
              printf("index_tau = %d\n", index_tau );
              index_tau_quad = 0;

              class_call(background_at_tau(pba,
                                           pgb2->tau_sampling_cls[index_tau],
                                           pba->long_info,
                                           pba->inter_normal,
                                           &last_index_vp,
                                           pvecback_vp),
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


                prefactor1 = k
                             /pvecback_vp[pba->index_bg_H]
                             /pvecback_vp[pba->index_bg_a];



                pgb2->first_order_sources_integ[pgb2->index_type_vp][index_l][index_tau][index_k_bessel] = prefactor1*pgb2->first_order_sources[pgb2->index_source_v][index_tau][index_k_bessel]*j_second_deriv;
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

                class_call(bessel_at_x_second_deriv(pgb2, pbs, x, index_l, &j_first_deriv), pbs->error_message, pgb2->error_message);
                class_call(bessel_at_x(pbs, x , index_l, &j), pbs->error_message, pgb2->error_message);

                prefactor1 = k
                             /pvecback_delta[pba->index_bg_H]
                             /pvecback_delta[pba->index_bg_a];



                double term1 = pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_bessel]*j;

                double term2 = prefactor1*pgb2->first_order_sources[pgb2->index_source_v][index_tau][index_k_bessel]*j_second_deriv;

                pgb2->first_order_sources_integ[pgb2->index_type_delta][index_l][index_tau][index_k_bessel] =
                      pgb2->first_order_sources_integ[pgb2->index_type_rsd][index_l][index_tau][index_k_bessel]
                      + pgb2->first_order_sources_integ[pgb2->index_type_density][index_l][index_tau][index_k_bessel];
              }
            }
          }
        }






    /* Lensing source term preparation */







    /* Lensing integral, at present this is the slowest part of the code. We check that this computation is required by the first if
      statement. */

    /*if (pgb2->index_type_lens != -1) {
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

            ///for (index_tau_lens = (index_tau * bessel_boost); index_tau_lens < pgb2->tau_size_bessel-1; index_tau_lens++) {
            for (index_tau_lens = (index_tau * bessel_boost); index_tau_lens < pgb2->tau_size_bessel-1; index_tau_lens++) {
              // NOTE that we skip the last index "pgb2->tau_size_bessel-1" to avoid the singularity. This would analytically drop out anyways
              // for l not 1.
              printf("index_k_bessel x index_tau_lens x index_tau %dx%dx%d\n", index_k_bessel, index_tau_lens, index_tau );
              tau_lens = pgb2->tau_sampling_bessel[index_tau_lens];

              if(x>pbs->x_max ){
              printf("ALERT! x= %g \n",x);continue;}

              class_call(bessel_at_x(pbs,pgb2->k_bessel[index_k_bessel]*(pba->conformal_age -tau_lens), index_l, &j), pbs->error_message, pgb2->error_message);

              index_of_tau_sampling_cls(tau_lens, &index_tau_cls, pgb2);

              first_index_tau_in_lens = index_tau * bessel_boost;

              /* If we are taking the first or last trapezoidal weight then take the first weight, else take the other weight */
        /*      if (index_tau_lens == first_index_tau_in_lens || index_tau_lens == (pgb2->tau_size_bessel-1) ){
                weight = (pgb2->tau_sampling_bessel[first_index_tau_in_lens]-pgb2->tau_sampling_bessel[pgb2->tau_size_bessel-1])/(2.0*((pgb2->tau_size_bessel-1)-first_index_tau_in_lens));
              }

              else{
                weight = (pgb2->tau_sampling_bessel[first_index_tau_in_lens]-pgb2->tau_sampling_bessel[pgb2->tau_size_bessel-1])/((pgb2->tau_size_bessel-1)-first_index_tau_in_lens);
              }

              lensing_result += ((pba->conformal_age - pgb2->tau_sampling_bessel[index_tau_lens])-(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]))
                * j * weight * pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][index_tau_lens][index_k_bessel]/((pba->conformal_age - pgb2->tau_sampling_bessel[index_tau_lens])*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]));
              //  lensing_result += ((pba->conformal_age - pgb2->tau_sampling_bessel2[][index_tau_lens])-(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]))
                //  * j * weight * pgb2->first_order_sources_integrand[pgb2->index_source_phi_plus_psi][index_tau_lens][index_k_bessel]/((pba->conformal_age - pgb2->tau_sampling_bessel[index_tau_lens])*(pba->conformal_age - pgb2->tau_sampling_cls[index_tau]));


            }

            pgb2->first_order_sources_integ[pgb2->index_type_lens][index_l][index_tau][index_k_bessel] = lensing_result;

          }
        }
      }
    }*/
    double j_lens;
    int index_k_lens2;


    if (pgb2->index_type_lens != -1) {


      int index_tau_quad;
      int last_index_k;
      int index_k;
      for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
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

              lens_sum += (-1.)*ptr->l[index_l]*(ptr->l[index_l]+1)*(2.-5*ptr->s_bias)*chi_fraction*phi_plus_psi
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
      for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++) {
      //for (int index_l = ptr->l_size[ppt->index_md_scalars]-2; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
        index_k = 0;
        index_k_lens2 = 0;
        for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {

          double k = pgb2->k_bessel[index_k_bessel];
          printf("index_k: %d, k = %g\n",index_k_bessel, k);
          class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);

          for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
            double g4_sum = 0.;
            printf("inside lens loop: index_l x index_k_bessel x index_tau = %dx%dx%d\n", index_l, index_k_bessel, index_tau);

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
      //double phi_plus_psi_prime;
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
          //printf("index_k: %d, k = %g\n",index_k_bessel, k);
          class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);

          for (int index_tau = 0; index_tau < pgb2->tau_size_cls; index_tau++){
            double g5_sum = 0.;
            printf("inside lens loop: index_l x index_k_bessel x index_tau = %dx%dx%d\n", index_l, index_k_bessel, index_tau);
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
              //printf("inside lens loop, tau = %g\n", pgb2->tau_sampling_bessel2[index_tau][index_tau_bessel]);

              class_call(index_of_tau_sampling_quadsources(tau, &index_tau_quad, ppt), pgb2->error_message, pgb2->error_message);

              double k = pgb2->k_bessel[index_k_bessel];
            //  printf("k = %g\n",k);
              class_call(index_of_k(k, &index_k, ppt), pgb2->error_message, pgb2->error_message);



              //printf("pgb2->tau_sampling_bessel2[%d][%d] = %g\n", index_tau, index_tau_bessel, tau);

              //index_of_tau_sampling_cls(tau, &index_tau_cls, pgb2);










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


    /*double * pvecbackg;
    int last_index_g = 0;
    class_alloc(pvecbackg, pba->bg_size*sizeof(double), pba->error_message);

    /* Prepare the two integrated GR terms (g4 and g5) for the integral function */
    /*if ((pgb2->index_type_g5 != -1) || (pgb2->index_type_g4 != -1)) {
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
            /*  if (index_tau_lens == first_index_tau_in_lens || index_tau_lens == (pgb2->tau_size_bessel-1) ){
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
    }*/
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
  printf("Starting SONG style fixed grid integrations..\n");
  double source_interp1;
  double source_interp2;
  double tau_first;
  double tau_second;
  int index_of_cls1;
  int index_of_cls2;
  double Pk_song;
  for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    for(int index_type_second = 0; index_type_second < pgb2->type_size; index_type_second++){
    //Warning
    //for(int index_type_first = 0; index_type_first < 1; index_type_first++){
      //for(int index_type_second = 1; index_type_second < 2; index_type_second++){
      for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
        for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
          for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
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

    /* We now wish to intepolate the Cl and Dl arrays across a range of multipoles using a spline */
 // herehere
    double * restrict_array;
    double * splined_array;
    int splined_array_size = 300;
    class_alloc(restrict_array, ptr->l_size[ppt->index_md_scalars] * sizeof(double), pba->error_message);
    class_alloc(splined_array, splined_array_size * sizeof(double), pba->error_message);
    for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
        restrict_array[index_l] = (ptr->l[index_l])*(ptr->l[index_l]);
    }
    printf("restrict_array[splined_array_size-1] = %g\n", restrict_array[splined_array_size-1]);




     /*for(int index_l = 0; index_l < splined_array_size; index_l++){
       double l = index_l+2.0;
       printf("l = %g\n", l );

       class_call(array_interpolate_spline(ptr->l,
                                ptr->l_size[ppt->index_md_scalars],
                                ptr->l,
                                ptr->l,
                                1,
                                l,
                                &last_index_l,
                                ptr->l,//splined_array,
                                splined_array_size, /** from 1 to n_columns */
                /*                pgb2->error_message),
                                pgb2->error_message,
                                pgb2->error_message);
       /*class_call(array_interpolate_two(ptr->l,
                                        1,
                                        index_l,   /** from 0 to (n_columns_x-1) */
                        /*                restrict_array,
                                        1,
                                        ptr->l_size[ppt->index_md_scalars],  /** must be the same for array_x and array_y */
                                /*        l,
                                        splined_array,
                                        splined_array_size, /** from 1 to n_columns_y */
                                  /*      pgb2->error_message),
                                        pgb2->error_message,
                                        pgb2->error_message);

        /*class_call(array_interpolate_two_arrays_one_column(
        					    ptr->l, /* assumed to be a vector (i.e. one column array) */
        	/*				    restrict_array,
        					    1,
        					    index_l, /* between 0 and (n_columns_y-1) */
        		/*			    ptr->l_size[ppt->index_md_scalars],  /** must be the same for array_x and array_y */
        			/*		    l,
        					    splined_array,
                      pgb2->error_message),
                      pgb2->error_message,
                      pgb2->error_message);*/
      //  printf("l = %g,   splined_array[%d] = *%g*      (%g) \n", l, index_l, splined_array, l*l);
    // }
    // exit(0);



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
    printf("#alpha size = %d\n", pgb2->alpha_size);
    printf("#k_max = %g\n",pgb2->k_bessel[pgb2->k_size_bessel-1]);
    if (alpha_log_sampling == 1) {
      printf("#alpha is log sampled: log spacing base = %g\n",alpha_log_base );
    }
    if (alpha_log_sampling != 1) {
      printf("#alpha is linearly sampled\n");
    }
    printf("#r size = %d\n", pgb2->r_size );
    printf("#log base: %g\n",log_b );
    printf("#epsilon: %g\n", epsilon);
    printf("#alpha_window_size = %d\n", pgb2->alpha_window_size);
    printf("#alpha_log_size = %d\n",pgb2->alpha_log_size);
    printf("#tau_size_selection = %d\n", pgb2->tau_size_selection);
    printf("#tau_size_cls= %d\n", pgb2->tau_size_cls);
    printf("#bessel_boost = %d\n", bessel_boost );
    printf("#z1=%g, z2=%g (%g) Gaussian\n", ppt->selection_mean[0],  ppt->selection_mean[0],ppt->selection_width[0]);
    printf("#renorm_factor[%d][%d] =%g\n",0,0,renorm_factor[0][0] );
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

    type_A = pgb2->index_type_density;
    type_B = pgb2->index_type_rsd;
    type_C = pgb2->index_type_delta;
    type_D = pgb2->index_type_delta;
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
    for (int index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++) {
      for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_selection; index_tau_first++) {
        for (int index_tau_third = 0; index_tau_third < pgb2->tau_size_selection; index_tau_third++) {
          for (int bin1 = 0; bin1 < 1; bin1++) {
            for (int bin2 = 0; bin2 < 1; bin2++) {
              for (int bin3 = 0; bin3 < 1; bin3++) {
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
                for(int l = 2; l < 200; l++){
                //  printf("last_index_l_DBC13, last_index_l_DAD12, last_index_l_DBD13, last_index_l_DAC12, last_index_l_DBC23, last_index_l_DAD21 = %d, %d, %d, %d, %d, %d\n", last_index_l_DBC13, last_index_l_DAD12, last_index_l_DBD13, last_index_l_DAC12, last_index_l_DBC23, last_index_l_DAD21);
                  //printf("last_index_l_DBD23,last_index_l_DAC21,last_index_l_DBC31,last_index_l_DAD32,last_index_l_DBD31,last_index_l_DAC32 = = %d, %d, %d, %d, %d, %d\n", last_index_l_DBD23,last_index_l_DAC21,last_index_l_DBC31,last_index_l_DAD32,last_index_l_DBD31,last_index_l_DAC32);
                  //printf("l = %d\n", l);
                  double l_squared = l*l;
                  /*index_of_l(l,
                             &index_l_found,
                             &last_index_l2,
                             ppt,
                             ptr);*/
                  //printf("l = %d closest = %d\n", l, ptr->l[index_l_found]);
                  interpolate_Dl_for_l(l,
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

                   interpolate_Dl_for_l(2*l,
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


                    interpolate_Dl_for_l(l,
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
                   interpolate_Dl_for_l(2*l,
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

                    interpolate_Dl_for_l(l,
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

                     interpolate_Dl_for_l(2*l,
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


                      interpolate_Dl_for_l(l,
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

                     interpolate_Dl_for_l(2*l,
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


                    interpolate_Dl_for_l(l,
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

                     interpolate_Dl_for_l(l,
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


                      interpolate_Dl_for_l(l,
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

                     interpolate_Dl_for_l(l,
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


                  /*double spline2 = pgb2->Dl[0][0][index_l_found-1][bin1][bin2][index_tau_first][index_tau_second]
                    *(ptr->l[index_l_found]-l)+pgb2->Dl[0][0][index_l_found][bin1][bin2][index_tau_first][index_tau_second]
                      *(l-ptr->l[index_l_found-1]);*/
                  /*double spline = restrict_array[index_l_found-1]*(ptr->l[index_l_found]-l)
                                + restrict_array[index_l_found]*(l-ptr->l[index_l_found-1]);*/

                  //spline2 /= (ptr->l[index_l_found] - ptr->l[index_l_found-1]);

                  //double error = (l_squared-spline)*100/l_squared;
                //  printf("Dl_minus = %g\n",pgb2->Dl[0][0][index_l_found-1][bin1][bin2][index_tau_first][index_tau_second] );
                  //printf("%d    %g\n", l, l*(l+1)*spline/2/_PI_);
                  printf("%d    %g\n", l, l*l*(DAC12*DBD13+DBC12*DAD13+DAC21*DBD23+DBC21*DAD23+DAC32*DBD31+DBC32*DAD31));
                  //printf("Dl_plus = %g\n",pgb2->Dl[0][0][index_l_found][bin1][bin2][index_tau_first][index_tau_second] );
                }
              }
            }
          }
        }
      }
    }
    exit(0);

    for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
    //for (int index_l = ptr->l_size[ppt->index_md_scalars]-2; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
    //for (int index_l = 3; index_l < 4; index_l++) {
      printf("%d    %d    %g\n",index_l, ptr->l[index_l],
                          ptr->l[index_l]*(ptr->l[index_l]+1)*(pgb2->Dl[pgb2->index_type_density][pgb2->index_type_density][index_l][0][0][0][0])/(2*_PI_));
    }





    printf("BULKED Code\n");










  double tau1_auto,tau2_auto;
  double window1_auto;
  double window2_auto;
  int index_tau_first_auto, index_tau_second_auto;

  double tau1_hard,tau2_hard;
  double window1_hard;
  double window2_hard;
  int index_tau_first_hard, index_tau_second_hard;
  double sum1_auto;
  double sum2_auto;
  double sum1_hard;
  double sum2_hard;
  double type1k0, type2k0, type1k5, type2k5, type1k8, type2k8;

  FILE *k2 = fopen("k2.dat", "w");
    if (k2 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

  FILE *k125 = fopen("k125.dat", "w");
    if (k125 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

  FILE *k275 = fopen("k275.dat", "w");
    if (k275 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

  FILE *file3 = fopen("D2.dat", "w");
    if (file3 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
  FILE *file4 = fopen("D125.dat", "w");
    if (file4 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
  FILE *file5 = fopen("D275.dat", "w");
    if (file5 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

  FILE *d2r1 = fopen("D2r1.dat", "w");

    if (d2r1 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
  FILE *d2r2 = fopen("D2r2.dat", "w");
    if (d2r2 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
  FILE *d2r3 = fopen("D2r3.dat", "w");
    if (d2r3 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

  FILE *d125a1 = fopen("D125a1.dat", "w");
    if (d125a1 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
  FILE *d125a2 = fopen("D125a2.dat", "w");
    if (d125a2 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
  FILE *d125a3 = fopen("D125a3.dat", "w");
    if (d125a3 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

  FILE *d275r1 = fopen("D275r1.dat", "w");
    if (d275r1 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
  FILE *d275r2 = fopen("D275r2.dat", "w");
    if (d275r2 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
  FILE *d275r3 = fopen("D275r3.dat", "w");
    if (d275r3 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

  printf("pgb2->tau_sampling_cls[0] = %g\n",pgb2->tau_sampling_cls[0]);
  printf("pgb2->tau_sampling_cls[end] = %g\n",pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]);

    /*Turn this on to print square contour map of tau_sampling_cls grid*/
  /*for (int index_tau_first = 0; index_tau_first < pgb2->tau_size_cls; index_tau_first++) {
    printf("%d/%d\n", index_tau_first, pgb2->tau_sampling_cls-1);
    for (int index_tau_second = 0; index_tau_second < pgb2->tau_size_cls; index_tau_second++) {
      double sumk0 = 0.;
      double sumk5 = 0.;
      double sumk8 = 0;
      for(int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++){
        printf("%dx%d\n",  index_tau_first, index_tau_second );
        /* Call primordial power spectrum (Fourier space)*/
      /*  class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, pgb2->k_bessel[index_k_bessel], &Pk), ppm->error_message, pgb2->error_message);


        type1k0 = pgb2->first_order_sources_integ[0][0][index_tau_first][index_k_bessel];

        type2k0 = pgb2->first_order_sources_integ[0][0][index_tau_second][index_k_bessel];

        type1k5 = pgb2->first_order_sources_integ[0][5][index_tau_first][index_k_bessel];

        type2k5 = pgb2->first_order_sources_integ[0][5][index_tau_second][index_k_bessel];

        type1k8 = pgb2->first_order_sources_integ[0][8][index_tau_first][index_k_bessel];

        type2k8 = pgb2->first_order_sources_integ[0][8][index_tau_second][index_k_bessel];

        //printf("source start =%g, source middle = %g, source end %g\n", pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_firstq-1][index_k_bessel], type1_interpq, type1q);

        /* Using area of trapezoid, we can sum a number of areas of trapezoids to approximate the integral */
      /*  sumk0 += pow(pgb2->k_bessel[index_k_bessel],-1.0) * 4. * _PI_ *  Pk * type1k0 * type2k0  * pgb2->w_trapz_k[index_k_bessel];

        sumk5 += pow(pgb2->k_bessel[index_k_bessel],-1.0) * 4. * _PI_ *  Pk * type1k5 * type2k5  * pgb2->w_trapz_k[index_k_bessel];

        sumk8 += pow(pgb2->k_bessel[index_k_bessel],-1.0) * 4. * _PI_ *  Pk * type1k8 * type2k8  * pgb2->w_trapz_k[index_k_bessel];

      }


      fprintf(file3,"%g       %g       %g\n", pgb2->tau_sampling_cls[index_tau_first], pgb2->tau_sampling_cls[index_tau_second], sumk0);

      fprintf(file4,"%g       %g       %g\n", pgb2->tau_sampling_cls[index_tau_first], pgb2->tau_sampling_cls[index_tau_second], sumk5);

      fprintf(file5,"%g       %g       %g\n", pgb2->tau_sampling_cls[index_tau_first], pgb2->tau_sampling_cls[index_tau_second], sumk8);



    }
  }*/
  int index_r1 = ceil(0.25*(pgb2->r_size+1.)-1.);
  int index_r_middle = ceil(0.5*(pgb2->r_size+1.)-1.);
  int index_r3 = ceil(0.75*(pgb2->r_size+1.)-1.);

  printf("index_r1 = %d (r=%g)\n",index_r1, r_bins[0][0][0][index_r1] );

  printf("index_r_middle = %d (r=%g)\n",index_r_middle,r_bins[0][0][0][index_r_middle] );
  printf("index_r3 = %d (r = %g)\n",index_r3,r_bins[0][0][0][index_r3] );

  /* This bin is for the use of pgb2->Dl3 which is the dirac-Angular power spectrum that is indexed with bins, such that bin1,bin2,alpha and r
  determine tau1 and tau2 of the correlation */

// Dl3
  double type1q, type2q, type1_interpq, type2_interpq, integq, tau_one_polarq, tau_two_polarq;
  int index_tau_first_polarq, index_tau_second_polarq;
  double window_int_check2q;
  double window1q, window2q,window1_interpq, window2_interpq, inner_integrationq, temp2q, outer_integrationq;
  double tau_oneq, tau_twoq;
  int index_tau_firstq, index_tau_secondq;
  double integrand_no_window;
  double testa;
  double k_integrand;
  int last_index_mock1;
  int last_index_mock2;
  double * pvecback_rsd1;
  double * pvecback_rsd2;
  class_alloc(pvecback_rsd1, pba->bg_size * sizeof(double), pba->error_message);
  class_alloc(pvecback_rsd2, pba->bg_size * sizeof(double), pba->error_message);
  double first_order_sources1_interp;
  double first_order_sources2_interp;
  double integrand;
  double j_second_deriv1;
  double j_second_deriv2;
  double prev_outer_integration;
  for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    for (int k = 0; k < index_type_first+1; k++){
      printf("X");
    }
    printf("\n");

    for(int index_type_second = 0; index_type_second < index_type_first+1; index_type_second++){
      for (int bin1 = 0; bin1 < ppt->selection_num; bin1++) {
        for (int bin2 = 0; bin2 < ppt->selection_num; bin2++) {
          printf("(bin1,bin2) = (%d,%d)\n", bin1,bin2);
          //for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
          for (int index_l = 7; index_l < 9; index_l++) {

            printf("index_l = %d (l=%d)\n", index_l, ptr->l[index_l]);
            outer_integrationq = 0.0;
            prev_outer_integration = 0.0;
            window_int_check2q = 0.0;
            double counter = 0.0;
            for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {
            //for (int index_alpha = 0; index_alpha < middle_index_alpha+1; index_alpha++) {
            //for(int index_r = 0; index_r < pgb2->r_size; index_r++){
              printf("index_alpha = %d\n",index_alpha);
              inner_integrationq = 0.0;
              temp2q = 0.0;

              for(int index_r = 0; index_r < pgb2->r_size; index_r++){
              //for (int index_alpha = 0; index_alpha < pgb2->alpha_size; index_alpha++) {

                  /* code */

                tau_oneq = tau0*(1.0-r_bins[bin1][bin2][index_alpha][index_r]*sin(alpha2[bin1][bin2][index_alpha]));
                tau_twoq = tau0*(1.0-r_bins[bin1][bin2][index_alpha][index_r]*cos(alpha2[bin1][bin2][index_alpha]));
                //printf("tau_oneq = %g, tau_twoq = %g\n", tau_oneq, tau_twoq);


                /* Find index in the linear sampling_cls grid for interpolation of arrays defined with this grid. */

                index_of_tau_sampling_cls(tau_oneq, &index_tau_firstq, pgb2);

                index_of_tau_sampling_cls(tau_twoq, &index_tau_secondq, pgb2);

                class_call(background_at_tau(pba,
                                             tau_oneq,
                                             pba->long_info,
                                             pba->inter_normal,
                                             &last_index_mock1,
                                             pvecback_rsd1),
                                             pba->error_message,
                                             pgb2->error_message);



                 class_call(background_at_tau(pba,
                                              tau_twoq,
                                              pba->long_info,
                                              pba->inter_normal,
                                              &last_index_mock2,
                                              pvecback_rsd2),
                                              pba->error_message,
                                              pgb2->error_message);



                integq = 0.0;

                for(int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++){


                  /* Call primordial power spectrum (Fourier space)*/
                  class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, pgb2->k_bessel[index_k_bessel], &Pk), ppm->error_message, pgb2->error_message);

                  double x1 = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - tau_oneq);
                  double x2 = pgb2->k_bessel[index_k_bessel]*(pba->conformal_age - tau_twoq);

                  double prefactor_rsd1 = 1.0
                               /pvecback_rsd1[pba->index_bg_H]
                               /pvecback_rsd1[pba->index_bg_a];

                  double prefactor_rsd2 = 1.0
                               /pvecback_rsd2[pba->index_bg_H]
                               /pvecback_rsd2[pba->index_bg_a];

                  class_call(bessel_at_x_second_deriv(pgb2, pbs, x1, index_l, &j_second_deriv1), pbs->error_message, pgb2->error_message);
                  //printf("called bessel1\n" );
                  class_call(bessel_at_x_second_deriv(pgb2, pbs, x2, index_l, &j_second_deriv2), pbs->error_message, pgb2->error_message);


                  //class_call(bessel_at_x(pbs, x1 , index_l &j1), pbs->error_message, pgb2->error_message);
                  //class_call(bessel_at_x(pbs, x2 , index_l, &j2), pbs->error_message, pgb2->error_message);

                  //pgb2->first_order_sources_integ[pgb2->index_type_density][index_l][index_tau][index_k_bessel] = pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau][index_k_bessel] * j;
                  //printf("called bessel2\n" );
                  type1q = pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_firstq][index_k_bessel];

                  type2q = pgb2->first_order_sources_integ[index_type_second][index_l][index_tau_secondq][index_k_bessel];

                  /* Interpolate between the two indices found on the tau_sampling_cls grid to find the exact tau_one and tau_two */

                  first_order_sources1_interp = pgb2->first_order_sources[pgb2->index_source_theta][index_tau_firstq-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_tau_firstq]-tau_oneq)
                                              + pgb2->first_order_sources[pgb2->index_source_theta][index_tau_firstq][index_k_bessel]*(tau_oneq-pgb2->tau_sampling_cls[index_tau_firstq-1]);
                  first_order_sources1_interp /= (pgb2->tau_sampling_cls[index_tau_firstq] - pgb2->tau_sampling_cls[index_tau_firstq-1]);

                  first_order_sources2_interp = pgb2->first_order_sources[pgb2->index_source_theta][index_tau_secondq-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_tau_secondq]-tau_twoq)
                                              + pgb2->first_order_sources[pgb2->index_source_theta][index_tau_secondq][index_k_bessel]*(tau_twoq-pgb2->tau_sampling_cls[index_tau_secondq-1]);
                  first_order_sources2_interp /= (pgb2->tau_sampling_cls[index_tau_secondq] - pgb2->tau_sampling_cls[index_tau_secondq-1]);

                  /*first_order_sources1_interp = pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau_firstq-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_tau_firstq]-tau_oneq)
                                              + pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau_firstq][index_k_bessel]*(tau_oneq-pgb2->tau_sampling_cls[index_tau_firstq-1]);
                  first_order_sources1_interp /= (pgb2->tau_sampling_cls[index_tau_firstq] - pgb2->tau_sampling_cls[index_tau_firstq-1]);

                  first_order_sources2_interp = pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau_secondq-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_tau_secondq]-tau_twoq)
                                              + pgb2->first_order_sources[pgb2->index_source_delta_cdm][index_tau_secondq][index_k_bessel]*(tau_twoq-pgb2->tau_sampling_cls[index_tau_secondq-1]);
                  first_order_sources2_interp /= (pgb2->tau_sampling_cls[index_tau_secondq] - pgb2->tau_sampling_cls[index_tau_secondq-1]);*/

                  type1_interpq = pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_firstq-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_tau_firstq]-tau_oneq)
                                + pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_firstq][index_k_bessel]*(tau_oneq-pgb2->tau_sampling_cls[index_tau_firstq-1]);
                  type1_interpq /= (pgb2->tau_sampling_cls[index_tau_firstq] - pgb2->tau_sampling_cls[index_tau_firstq-1]);

                  type2_interpq = pgb2->first_order_sources_integ[index_type_second][index_l][index_tau_secondq-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_tau_secondq]-tau_twoq)
                                + pgb2->first_order_sources_integ[index_type_second][index_l][index_tau_secondq][index_k_bessel]*(tau_twoq-pgb2->tau_sampling_cls[index_tau_secondq-1]);

                  type2_interpq /= (pgb2->tau_sampling_cls[index_tau_secondq] - pgb2->tau_sampling_cls[index_tau_secondq-1]);

                  //printf("k = %g\n", pgb2->k_bessel[index_k_bessel]);
                  //printf("tau_minus = %g, tau_exact = %g, tau_plus = %g\n", pgb2->tau_sampling_cls[index_tau_firstq-1], tau_oneq, pgb2->tau_sampling_cls[index_tau_firstq]);
                  //printf("type1_minus = %g, type1_interpq = %g, type1_plus = %g\n", pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_firstq-1][index_k_bessel], type1_interpq, pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_firstq][index_k_bessel]);
                  //printf("source start =%g, source middle = %g, source end %g\n", pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_firstq-1][index_k_bessel], type1_interpq, type1q);

                  /* Using area of trapezoid, we can sum a number of areas of trapezoids to approximate the integral */

                  integq +=  4. * _PI_ *  Pk * pow(pgb2->k_bessel[index_k_bessel],-1.0)*first_order_sources1_interp * first_order_sources2_interp
                           * prefactor_rsd1 * prefactor_rsd2 * j_second_deriv1 * j_second_deriv2 * pgb2->w_trapz_k[index_k_bessel];

                //  integq +=  4. * _PI_ *  Pk * pow(pgb2->k_bessel[index_k_bessel],-1.0)*first_order_sources1_interp * first_order_sources2_interp
                          //    * j1 * j2 * pgb2->w_trapz_k[index_k_bessel];

                  //integq +=  4. * _PI_ *  Pk * pow(pgb2->k_bessel[index_k_bessel],-1.0)*type1_interpq * type2_interpq * pgb2->w_trapz_k[index_k_bessel];

                  k_integrand = 4. * _PI_ *  Pk * pow(pgb2->k_bessel[index_k_bessel],-1.0)*first_order_sources1_interp * first_order_sources2_interp
                                * prefactor_rsd1 * prefactor_rsd2 * j_second_deriv1 * j_second_deriv2 ;

                  /*if (counter == 0.0 && index_k_bessel == pgb2->k_size_bessel-1) {

                    //printf("%g        %g\n", pgb2->k_bessel[index_k_bessel], k_integrand);
                    printf("%d        %g\n", ptr->l[index_l], ptr->l[index_l]*(ptr->l[index_l]+1)*integq/(2*_PI_));
                    //printf("%g      %g\n", pgb2->k_bessel[index_k_bessel], Pk);
                    //printf("tau_minus = %g, tau_exact = %g, tau_plus = %g\n", pgb2->tau_sampling_cls[index_tau_firstq-1], tau_oneq, pgb2->tau_sampling_cls[index_tau_firstq]);
                    //printf("type1_minus = %g, type1_interpq = %g, type1_plus = %g\n", pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_firstq-1][index_k_bessel], type1_interpq, pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_firstq][index_k_bessel]);
                    //printf("source start =%g, source middle = %g, source end %g\n", pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_firstq-1][index_k_bessel], type1_interpq, type1q);
                    // WARNING k_integrand is changed
                    counter += 1.0;
                  }

                  if (index_l == 0 && index_alpha == middle_index_alpha && index_r == index_r_find) {
                    fprintf(k2,"%g       %g\n", pgb2->k_bessel[index_k_bessel],  k_integrand);
                  }

                  if (index_l == 5 && index_alpha == middle_index_alpha && index_r == index_r_find) {
                    fprintf(k125,"%g      %g\n", pgb2->k_bessel[index_k_bessel],  k_integrand);
                    counter += 1.0;
                  }

                  if (index_l == 8 && index_alpha == middle_index_alpha && index_r == index_r_find) {
                    fprintf(k275,"%g      %g\n", pgb2->k_bessel[index_k_bessel],  k_integrand);
                  }*/

                  //pgb2->k_integrand[index_l][bin1][bin2][index_alpha][index_r][index_k_bessel] = pow(pgb2->k_bessel[index_k_bessel],-1.0) * 4. * _PI_ *  Pk * type1_interpq * type2_interpq;

                }
                //printf("integq = %g\n", integq);

                pgb2->Dl3[index_type_first][index_type_second][index_l][bin1][bin2][index_alpha][index_r] = integq;

                /* Find indices in the selection grid */
                //printf("tau_oneq = %g\n",tau_oneq );
                //printf("tau_twoq = %g\n",tau_twoq );

                index_of_tau_sampling_selection(tau_oneq, bin1, &index_tau_first_polarq, pgb2);
                index_of_tau_sampling_selection(tau_twoq, bin2, &index_tau_second_polarq, pgb2);

                //printf("index_tau_first_polarq = %d\n",index_tau_first_polarq);
                //printf("index_tau_second_polarq = %d\n", index_tau_second_polarq );


                /* Interpolate between the two indices found on the tau_sampling_cls grid to find the exact tau_one and tau_two */


                window1_interpq = selection[bin1][index_tau_first_polarq-1]*(pgb2->tau_sampling_selection[bin1][index_tau_first_polarq]-tau_oneq)
                              + selection[bin1][index_tau_first_polarq]*(tau_oneq-pgb2->tau_sampling_selection[bin1][index_tau_first_polarq-1]);
                window1_interpq /= (pgb2->tau_sampling_selection[bin1][index_tau_first_polarq] - pgb2->tau_sampling_selection[bin1][index_tau_first_polarq-1]);

                window2_interpq = selection[bin2][index_tau_second_polarq-1]*(pgb2->tau_sampling_selection[bin2][index_tau_second_polarq]-tau_twoq )
                              + selection[bin2][index_tau_second_polarq]*(tau_twoq -pgb2->tau_sampling_selection[bin2][index_tau_second_polarq-1]);
                window2_interpq /= (pgb2->tau_sampling_selection[bin2][index_tau_second_polarq] - pgb2->tau_sampling_selection[bin2][index_tau_second_polarq-1]);

                //printf("window[%d][%d] = %g\n", bin1, index_tau_first_polarq-1,selection[bin1][index_tau_first_polarq-1]);
              //  printf("%g    %g\n", tau_oneq, window1_interpq);
              //  printf("window[%d][%d] = %g\n", bin1, index_tau_first_polarq, selection[bin1][index_tau_first_polarq]);
                //printf("%g      %g\n", r_bins[bin1][bin2][index_alpha][index_r], window1_interpq*window2_interpq);

                //inner_integrationq += tau0*tau0*r_bins[bin1][bin2][index_alpha][index_r]*window_pair[bin1][bin2][index_alpha][index_r]*w_trapz_r_bins[bin1][bin2][index_alpha][index_r]
                                      //  *pgb2->Dl3[index_type_first][index_type_second][index_l][bin1][bin2][index_alpha][index_r];

                inner_integrationq += tau0*tau0*r_bins[bin1][bin2][index_alpha][index_r]*renorm_factor[bin1][bin2]*window1_interpq*window2_interpq*w_trapz_r_bins[bin1][bin2][index_alpha][index_r]
                                      *pgb2->Dl3[index_type_first][index_type_second][index_l][bin1][bin2][index_alpha][index_r];

                //printf("WI: %g        %g        %g        %g        %g\n",alpha2[bin1][bin2][index_alpha], r_bins[bin1][bin2][index_alpha][index_r], tau_oneq, tau_twoq, renorm_factor[bin1][bin2]*window1_interpq*window2_interpq);
                //printf("WE: %g        %g        %g        %g        %g       %g\n",alpha2[bin1][bin2][index_alpha], r_bins[bin1][bin2][index_alpha][index_r], tau_oneq, tau_twoq, window_pair[bin1][bin2][index_alpha][index_r], (window_pair[bin1][bin2][index_alpha][index_r]-(renorm_factor[bin1][bin2]*window1_interpq*window2_interpq*100.))/window_pair[bin1][bin2][index_alpha][index_r]);

                integrand = tau0*tau0*r_bins[bin1][bin2][index_alpha][index_r]*renorm_factor[bin1][bin2]*window1_interpq*window2_interpq
                                        *pgb2->Dl3[index_type_first][index_type_second][index_l][bin1][bin2][index_alpha][index_r];

                integrand_no_window = tau0*tau0*r_bins[bin1][bin2][index_alpha][index_r]
                                                          *pgb2->Dl3[index_type_first][index_type_second][index_l][bin1][bin2][index_alpha][index_r];
                /*if (index_l == 0) {
                  fprintf(file3,"%g       %g       %g\n", tau_oneq, tau_twoq, integrand);
                }

                if (index_l == 5) {
                  fprintf(file4,"%g       %g       %g\n", tau_oneq, tau_twoq, integrand);
                }

                if (index_l == 8) {
                  fprintf(file5,"%g       %g       %g\n", tau_oneq, tau_twoq, integrand);
                }*/


                //printf("r_bins[%d][%d][%d][%d] = %g\n", bin1, bin2, index_alpha, index_r, r_bins[bin1][bin2][index_alpha][index_r] );

                /*Dl2_integrand[index_type_first][index_type_second][index_l][index_alpha][index_r] = tau0*tau0
                                      *r_bins[bin1][bin2][index_alpha][index_r]*window_pair[bin1][bin2][index_alpha][index_r]
                                        *pgb2->Dl3[index_type_first][index_type_second][index_l][bin1][bin2][index_alpha][index_r];*/
                 /*This is just to check that when integrates the window1*window2*jacobi factor over dr and dalpha we yield 1.                       */
                temp2q += tau0*tau0*r_bins[bin1][bin2][index_alpha][index_r]*window_pair[bin1][bin2][index_alpha][index_r]*w_trapz_r_bins[bin1][bin2][index_alpha][index_r];
                //temp2q += tau0*tau0*r_bins[bin1][bin2][index_alpha][index_r]*window1_interpq*window2_interpq*w_trapz_r_bins[bin1][bin2][index_alpha][index_r];

                /*if (index_l == 0 && index_r == index_r1) {
                  fprintf(d2r1,"%g       %g\n", alpha2[bin1][bin2][index_alpha],  integrand);
                }

                if (index_l == 0 && index_r == index_r_middle ) {
                  fprintf(d2r2,"%g       %g\n", alpha2[bin1][bin2][index_alpha],  integrand);
                }

                if (index_l == 0 && index_r == index_r3) {
                  fprintf(d2r3,"%g       %g\n", alpha2[bin1][bin2][index_alpha],  integrand);
                }

                if (index_l == 5 && index_r == index_r1) {
                  fprintf(d125a1,"%g      %g\n", alpha2[bin1][bin2][index_alpha],  integrand);
                }

                if (index_l == 5 && index_alpha == log_index_middle ) {
                  fprintf(d125a2,"%g      %g\n", r_bins[bin1][bin2][index_alpha][index_r],  integrand);
                }

                if (index_l == 5 && index_r == index_r3) {
                  fprintf(d125a3,"%g      %g\n", alpha2[bin1][bin2][index_alpha],  integrand);
                }

                if (index_l == 8 && index_r == index_r1) {
                  fprintf(d275r1,"%g      %g\n", alpha2[bin1][bin2][index_alpha],  integrand);
                }

                if (index_l == 8 && index_r == index_r_middle) {
                  fprintf(d275r2,"%g      %g\n", alpha2[bin1][bin2][index_alpha],  integrand);
                }

                if (index_l == 8 && index_r == index_r3) {
                  fprintf(d275r3,"%g      %g\n", alpha2[bin1][bin2][index_alpha],  integrand);
                }*/

              }





              outer_integrationq += inner_integrationq*w_trapz_alpha2[bin1][bin2][index_alpha];
              //outer_integrationq += inner_integrationq*w_trapz_r_bins[bin1][bin2][0][index_r];
              //printf("%g        %g        %g\n", alpha2[bin1][bin2][index_alpha], outer_integrationq, outer_integrationq-prev_outer_integration);
              prev_outer_integration = outer_integrationq;
              window_int_check2q += temp2q*w_trapz_alpha2[bin1][bin2][index_alpha];

            }


            printf("window_int_check2q = %g\n", window_int_check2q);
            printf("*outer_integrationq = %g\n", outer_integrationq );
            pgb2->Cl3[index_type_first][index_type_second][index_l][bin1][bin2] = outer_integrationq;
            printf("pgb2->Cl3[%d][%d][%d][%d][%d] = %g\n", index_type_first, index_type_second, index_l, bin1, bin2, ptr->l[index_l]*(ptr->l[index_l]+1)*pgb2->Cl3[index_type_first][index_type_second][index_l][bin1][bin2]/(2*_PI_));

          }
        }
      }
    }
  }


  printf("# AUTOMATIC BIN PARAMETERISATION (pgb2->Dl3, Cl3)\n" );
  if (k_sampling == 0) {
    printf("#k is linearly sampled with %d points\n", pgb2->k_size_bessel );
  }
  if (k_sampling == 1) {
    printf("#k is log-sampled with %d points\n", pgb2->k_size_bessel );
    //printf("#k_logbase = %g\n",k_logbase);
  }
  printf("#alpha size = %d\n", pgb2->alpha_size);
  printf("#k_max = %g\n",pgb2->k_bessel[pgb2->k_size_bessel-1]);
  if (alpha_log_sampling == 1) {
    printf("#alpha is log sampled: log spacing base = %g\n",alpha_log_base );
  }
  if (alpha_log_sampling != 1) {
    printf("#alpha is linearly sampled\n");
  }
  printf("#r size = %d\n", pgb2->r_size );
  printf("#log base: %g\n",log_b );
  printf("#epsilon: %g\n", epsilon);
  printf("#alpha_window_size = %d\n", pgb2->alpha_window_size);
  printf("#alpha_log_size = %d\n",pgb2->alpha_log_size);
  printf("#tau_size_selection = %d\n", pgb2->tau_size_selection);
  printf("#tau_size_cls= %d\n", pgb2->tau_size_cls);
  printf("#z1=%g, z2=%g (%g) Gaussian\n", ppt->selection_mean[0],  ppt->selection_mean[0],ppt->selection_width[0]);
  printf("#renorm_factor[%d][%d] =%g\n",0,0,renorm_factor[0][0] );
  printf("#l    l(l+1)Cl3/2pi  l(l+1)Dl3/2pi \n");


  //for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
  for (int index_l = 7; index_l < 9; index_l++) {


    printf("%d    %g\n",ptr->l[index_l],
                        ptr->l[index_l]*(ptr->l[index_l]+1)*(pgb2->Cl3[0][0][index_l][0][0]/*+pgb2->Cl3[1][0][index_l][0][0]+pgb2->Cl3[1][1][index_l][0][0]*/)/(2*_PI_));/*+pgb2->Cl3[0][0][index_l][1][0]*///));
  }
  printf("BULKED code\n");
  exit(0);

  double interpolated_Dl;
  //double interpolated_kl;
  for (int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++) {
    double tau_minus = tau0*(1.-r_bins[0][0][middle_index_alpha][index_r_find]*sin(alpha2[0][0][middle_index_alpha]));
    double tau_plus = tau0*(1.-r_bins[0][0][middle_index_alpha][index_r_find-1]*cos(alpha2[0][0][middle_index_alpha]));

    /*if (index_l == 0) {
      for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
        double kl_minus = pgb2->k_integrand[index_l][0][0][middle_index_alpha][index_r_find][index_k_bessel];
        double kl_plus = pgb2->k_integrand[index_l][0][0][middle_index_alpha][index_r_find][index_k_bessel-1];
        interpolated_kl = kl_minus+(selection_mean_tau_bin1-tau_minus)*(kl_plus-kl_minus)/(tau_plus-tau_minus);
        fprintf(k2, "%d     %g\n", pgb2->k_bessel[index_k], interpolated_kl);
      }
    }

    if (index_l == 5) {
      for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
        double kl_minus = pgb2->k_integrand[index_l][0][0][middle_index_alpha][index_r_find][index_k_bessel];
        double kl_plus = pgb2->k_integrand[index_l][0][0][middle_index_alpha][index_r_find][index_k_bessel-1];
        interpolated_kl = kl_minus+(selection_mean_tau_bin1-tau_minus)*(kl_plus-kl_minus)/(tau_plus-tau_minus);
        fprintf(k125, "%d     %g\n", pgb2->k_bessel[index_k], interpolated_kl);
      }
    }

    if (index_l == 8) {
      for (int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++) {
        double kl_minus = pgb2->k_integrand[index_l][0][0][middle_index_alpha][index_r_find][index_k_bessel];
        double kl_plus = pgb2->k_integrand[index_l][0][0][middle_index_alpha][index_r_find][index_k_bessel-1];
        interpolated_kl = kl_minus+(selection_mean_tau_bin1-tau_minus)*(kl_plus-kl_minus)/(tau_plus-tau_minus);
        fprintf(k275, "%d     %g\n", pgb2->k_bessel[index_k], interpolated_kl);
      }
    }*/


    double Dl_minus = pgb2->Dl3[0][0][index_l][0][0][middle_index_alpha][index_r_find];
    double Dl_plus = pgb2->Dl3[0][0][index_l][0][0][middle_index_alpha][index_r_find-1];
    interpolated_Dl = Dl_minus*(tau_plus-selection_mean_tau_bin1)
                  + Dl_plus*(selection_mean_tau_bin1-tau_minus);
    interpolated_Dl /= (tau_plus - tau_minus);
    //interpolated_Dl = Dl_minus+(selection_mean_tau_bin1-tau_minus)*(Dl_plus-Dl_minus)/(tau_plus-tau_minus);
    printf("Dl_minus = %g, Dl_interpolated = %g, Dl_plus = %g\n", Dl_minus, interpolated_Dl, Dl_plus);
    printf("tau_minus = %g, tau_exact = %g, tau_plus = %g, \n", tau_minus, selection_mean_tau_bin1, tau_plus);
    printf("index_r_find =%d\n", index_r_find);

    printf("%d    %g\n",ptr->l[index_l],
                        ptr->l[index_l]*(ptr->l[index_l]+1)*(interpolated_Dl/(2*_PI_)));
  }
  fclose(file3);
  exit(0);


  /* Fill the array pgb2->Dl[index_type_first][index_type_second][index_l][index_tau_first][index_tau_second] by calling the integral() function */

  /* NOTE THIS ONE IS ORIGINAL pgb2->Dl2 and has been checked */
/*  int index_tau_first, index_tau_second;
  double type1_interp, type2_interp;
  for(int index_type_first = 0; index_type_first < pgb2->type_size; index_type_first++){
    for (int k = 0; k < index_type_first+1; k++){
      printf("X");
    }
    printf("\n");
// alert last tau width is missing
  /* Allocate backgroundmake song vectors */
  /*  double * backvec_z1;
    double * backvec_z2;

    class_alloc(backvec_z1, pba->bg_size * sizeof(double), pba->error_message);
    class_alloc(backvec_z2, pba->bg_size * sizeof(double), pba->error_message);

    int index_tau_first_test = 1000;
    for(int index_type_second = 0; index_type_second < index_type_first+1; index_type_second++){
      //alert we don't run over last index_l
      // ALERT LOOPING ONL UP TO FIRST TIME INDEX
      // ALERT l running
      //ALERT l start
      //for(int index_l = test_index_l; index_l < test_index_l+1/*ptr->l_size[ppt->index_md_scalars]-1*///; index_l++){
    /*  for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]-1; index_l++){
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
      /*      double z1 = pba->a_today/backvec_z1[pba->index_bg_a]-1.;
            //for(int index_tau_second = 0; index_tau_second < index_tau_first+1/*pgb2->tau_size_cls*///; index_tau_second++){
            //for(int index_tau_second = 0; index_tau_second < pgb2->tau_size_cls; index_tau_second++){
            /* Zero the integration sum */
          /*  class_call(background_at_tau(pba,
                                         pgb2->tau_sampling_cls[index_tau_second],
                                         pba->long_info,
                                         pba->inter_normal,
                                         &last_index,
                                         backvec_z2),
                                         pba->error_message,
                                         ptr->error_message);

            /* infer redhsift */
            /*double z2 = pba->a_today/backvec_z2[pba->index_bg_a]-1.;
            integ = 0.0;

            for(int index_k_bessel = 0; index_k_bessel < pgb2->k_size_bessel; index_k_bessel++){

              /* Call primordial power spectrum (Fourier space)*/
        //      class_call(primordial_spectrum_at_k(ppm, ppt->index_md_scalars, linear, pgb2->k_bessel[index_k_bessel], &Pk), ppm->error_message, pgb2->error_message);

              //printf("tau_one = %g\n",tau_one );
              /*printf("%dx%dx%dx%dx%dx%d\n",
              index_type_first,
              index_type_second,
              index_l,
              index_alpha,
              index_r,
              index_k_bessel);
              printf("index_tau_first = %d, index_tau_second = %d \n", index_tau_first, index_tau_second);*/

          //    type1 = pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_first][index_k_bessel];
              //printf("type1 = %g\n", type1);
              //printf("type2 = %g\n", type2);
        //      type2 = pgb2->first_order_sources_integ[index_type_second][index_l][index_tau_second][index_k_bessel];

              /* Interpolate between the two indices found on the tau_sampling_cls grid to find the exact tau_one and tau_two */

          /*    type1_interp = pgb2->first_order_sources_integ[index_type_first][index_l][index_tau_first-1][index_k_bessel]*(pgb2->tau_sampling_cls[index_tau_first]-tau_one)
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
        /*      integ += pow(pgb2->k_bessel[index_k_bessel],-1.0) * 4. * _PI_ *  Pk * type1_interp * type2_interp  * pgb2->w_trapz_k[index_k_bessel];
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

    /*      }
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
    */






    // Interpolate to read exact Dl(z1,z2) given
  /*  double z_test = 1.0;
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


  printf("Done k-integration.\n");

*/

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


/*  printf("Starting time-integration.\n");

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

/*  if (ppt->selection == dirac) {
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

              /*    window1_interp = selection[0][index_tau_first_polar-1]*(pgb2->tau_sampling_selection[0][index_tau_first_polar]-tau_one_polar)
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
            /*    }
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
    printf("# HARDCODED BIN PARAMETERISATION (pgb2->Dl2, Cl)\n" );
    if (k_sampling == 0) {
      printf("#k is linearly sampled with %d points\n", pgb2->k_size_bessel );
    }
    if (k_sampling == 1) {
      printf("#k is log-sampled with %d points\n", pgb2->k_size_bessel );
      //printf("#k_logbase = %g\n",k_logbase);
    }
    printf("#alpha size = %d\n", pgb2->alpha_size);
    printf("#k_max = %g\n",pgb2->k_bessel[pgb2->k_size_bessel-1]);
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
    */




    /*double dummy_inner;
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
          /*    for(index_tau_second = 0; index_tau_second < pgb2->tau_size_selection; index_tau_second++){

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


              /*    double alpha2 = atan((1.-(tau1/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1]))/(1.-(tau2/pgb2->tau_sampling_cls[pgb2->tau_size_cls-1])));

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
/*printf("integrating between k =%g and %g\n",pgb2->k_bessel[index_k_bessel],pgb2->k_bessel[pgb2->k_size_bessel-1]);
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

  /*  printf("%d    %g\n",
           ptr->l[index_l],
           norm*pgb2->Cl[pgb2->index_type_rsd][pgb2->index_type_rsd][index_l][0][0]);
    }



  printf("Starting bispectrum computation.\n");


/***************************************************************************
=====================   Bispectrum Computation =============================
****************************************************************************/

/*  int index_l1, index_l2, index_l3;

//NOTE: Should each dimension in the asym_redgalbispectrum array be allocated equally?
/* Allocation for the array pgb2->asym_redgalbispectrum[index_type_SO1][index_type_SO2][index_type_FO1][index_type_FO2][index_l2][index_l3][bin1][bin2][bin3]. This quantity is then summed
  with different permutations to yield the full reduced galaxy bispectrum.*/
  /*class_alloc(pgb2->asym_redgalbispectrum,
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

/*  printf("Allocating size %ix%ix%ix%ix%ix%ix%ix%ix%ix%i bytes \n", pgb2->type_size, pgb2->type_size, pgb2->type_size, ptr->l_size[ppt->index_md_scalars], ptr->l_size[ppt->index_md_scalars], ptr->l_size[ppt->index_md_scalars], ppt->selection_num, ppt->selection_num, ppt->selection_num);

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
