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
