/* Module to compute galaxy bispectra */


#include "galbispectra2.h"
#include "perturbations2.h"







/*****************************************************************
##############       DEFINE K SAMPLING   #########################
******************************************************************


/**
 * Determine the Fourier grid in (k1,k2,k3) that will be used to sample the line of
 * sight sources.
 *
 * At second order, the perturbations depend on three Fourier wavemodes rather than one,
 * due to mode coupling. In principle, this amounts to 9 degrees of freedom; in practice,
 * we can enforce the statistical isotropy of the Universe to reduce them to 3, which we
 * choose to be the magnitudes of the three wavevectors, (k1,k2,k3). We shall then solve
 * and build the line-of-sight sources on a 3D grid in (k1,k2,k3).
 *
 * Here we define the (k1,k2,k3) grid using an optimised algorithm way; for an average
 * precision run, the final list will consist of about 10^5 triplets (compare to the
 * ~400 values of k needed at first order).
 *
 * The algorithm we use is described in Sec. 5.3.2 of http://arxiv.org/abs/1405.2280;
 * mode coupling is described in Sec 3.5.2; details on the geometry of the wavemodes
 * and statistical isotropy can be found in Appendix B.
 *
 * This function will write the following fields in the ppt2 structure:
 *  - ppt2->k
 *  - ppt2->k_size
 *  - ppt2->k3
 *  - ppt2->k3_size
 *  - ppt2->k_min
 *  - ppt2->k_max
 *
 * and the following ones in the ppt structure:
 *  - ppt->k
 *  - ppt->k_size
 *  - ppt->k_size_cl
 *  - ppt->k_size_cmb
 *  - ppt->k_min
 *  - ppt->k_max
 */

int galbispectra2_get_k_lists (
          struct precision * ppr,
          struct precision2 * ppr2,
          struct background * pba,
          struct thermo * pth,
          struct perturbs * ppt,
          struct perturbs2 * ppt2
          )
{


  // ====================================================================================
  // =                                 k1 and k2 grid                                   =
  // ====================================================================================

  /* The dummy wavemodes k1 and k2 share a common grid. We compute it here and store it
  in ppt2->k. */


  // -------------------------------------------------------------------------------
  // -                             Logarithmic sampling                            -
  // -------------------------------------------------------------------------------

  if (ppt2->k_sampling == log_k_sampling) {

    class_test (ppr2->k_min_custom >= ppr2->k_max_custom,
      ppt2->error_message,
      "k_min must be smaller than k_max");

    ppt2->k_size = ppr2->k_size_custom;

    class_alloc (ppt2->k, ppt2->k_size*sizeof(double), ppt2->error_message);

    class_call (log_space (ppt2->k, ppr2->k_min_custom, ppr2->k_max_custom, ppr2->k_size_custom),
      ppt2->error_message, ppt2->error_message);

  }


  // -------------------------------------------------------------------------------
  // -                               Linear sampling                               -
  // -------------------------------------------------------------------------------

  else if (ppt2->k_sampling == lin_k_sampling) {

    class_test (ppr2->k_min_custom >= ppr2->k_max_custom,
      ppt2->error_message,
      "k_min must be smaller than k_max");

    ppt2->k_size = ppr2->k_size_custom;

    class_alloc (ppt2->k, ppt2->k_size*sizeof(double), ppt2->error_message);

    class_call (lin_space (ppt2->k, ppr2->k_min_custom, ppr2->k_max_custom, ppr2->k_size_custom),
      ppt2->error_message, ppt2->error_message);

  }


  // -------------------------------------------------------------------------------
  // -                               Smart sampling                                -
  // -------------------------------------------------------------------------------

  /* Adopt the same k-sampling algorithm as the one adopted in vanilla CLASS, which
  consists of two linear regimes, with the addition of a logarithmic leg for very
  small k in order to capture the evolution on large scales. This method is described
  in Sec. 5.3.2.1 of http://arxiv.org/abs/1405.2280. */

  else if (ppt2->k_sampling == smart_sources_k_sampling) {

    class_test(ppr2->k_step_transition == 0,
      ppt2->error_message,
      "stop to avoid division by zero");

    class_test(pth->rs_rec == 0,
      ppt2->error_message,
      "stop to avoid division by zero");

    /* Since we do not know the size of ppt2->k yet, we allocate it with the largest
    possible number of k-modes */
    class_alloc (ppt2->k, K_SIZE_MAX*sizeof(double), ppt2->error_message);

    /* The smallest k-mode in the sampling is chosen to be inversely proportional to tau_0,
    with the constant of proportionality given by the user via ppr2->k_min_tau0. For the
    CMB, this parameter should be much smaller than 1 to sample the small l with high
    precision. */
    double k_min = ppr2->k_min_tau0 / pba->conformal_age;

    /* Setup the first point in the k-sampling*/
    int index_k = 0;
    ppt2->k[index_k++] = k_min;
    double k = k_min;

    /* The value of the maximum wavemode depends on what we are computing. For CMB
    observables such as the angular power spectrum C_l and the bispectrum B_l1_l2_l3,
    we require k_max to be big enough to sample the smallest needed angular scale, which
    is determined by the ppt->l_scalar_max parameter; for an average run with l_max=3000,
    usually k_max ~ 0.4/Mpc is enough. For LSS observables such as the power spectrum P(k)
    and the bispectrum B_k1_k2_k3, k_max is usually of order 10/Mpc or larger and is
    determined directly by the user with the P_k_max_h parameters. */
    double k_max_cmb = 0;
    double k_max_lss = 0;

    /* Comoving scale corresponding to the sound horizon at recombination, usually of
    order 2*pi/150 ~ 0.04 Mpc^-1. CMB quantities vary in k with roughly this step. */
    double k_rec = 2 * _PI_ / pth->rs_rec;


    /* - LSS sampling */

    /* The LSS observables (P_k, B_k1_k2_k3...) usually require a less dense sampling in k
    than the CMB ones (C_l, B_l1_l2_l3...), but usually have a larger k_max. Following CLASS
    algorithm, we sample k using the denser CMB sampling up to k_max_cmb, and then use the
    logarithmic LSS sampling up to k_max_lss. NOTE: This block  is not entered if
    k < k_max_cmb, which means that the LSS parameters such as k_per_decade_for_pk and
    k_per_decade_for_bao are ignored for k < k_max_cmb. */

    if (ppt2->has_lss == _TRUE_) {

      k_max_lss = ppt2->k_max_for_pk;

      while (k < k_max_lss) {

        /* Use a logarithmic sampling with step ppr2->k_per_decade_for_pk. When close to the
        BAO peak (k ~ 2*pi/150 ~ 0.04 Mpc^-1) make the sampling more dense according to the
        ppr2->k_per_decade_for_bao parameter. */
        k *= pow (10, 1/(ppr2->k_per_decade_for_pk
                         + (ppr2->k_per_decade_for_bao-ppr2->k_per_decade_for_pk)
                           * (1-tanh(pow((log(k)-log(ppr2->k_bao_center*k_rec))/log(ppr2->k_bao_width),4)))));

        ppt2->k[index_k++] = k; /* Update the k-sampling array with the new value */

        class_test ((index_k+1) > K_SIZE_MAX,
          ppt2->error_message,
          "k size for LSS is too large; check k_per_decade_for_pk and k_per_decade_for_bao");

      }

    } // if LSS sampling


    ppt2->k_size = index_k; /* Number of points in the k-sampling */
    ppt2->k[ppt2->k_size-1] = k_max_lss; /* Last sampling point = exactly k_max */

    /* Free the excess memory we have allocated in ppt2->k */
    class_realloc(ppt2->k,
                  ppt2->k,
                  ppt2->k_size*sizeof(double),
                  ppt2->error_message);


  } // end of if(smart_sources_k_sampling)


  // -------------------------------------------------------------------------------
  // -                             Add output points                               -
  // -------------------------------------------------------------------------------

  /* The user might have asked to output the perturbations at specific configurations
  of (k1,k2,k3) using the k1_out, k2_out and k3_out parameters. Here we add these
  k-values to the list of computed k in SONG, that is, to ppt2->k. */

  if (ppt2->k_out_size > 0) {

    /* If we are running in k_out_mode, SONG will compute only the k values specified
    here, and ignore the previously added k. Therefore, we erase the k grid completely */
    if (ppt2->k_out_mode == _TRUE_)
      ppt2->k_size = 0;

    /* Build a 1D array with the output values for k1 and k2 */
    double * k_out;
    class_alloc (k_out, (2*ppt2->k_out_size+2)*sizeof(double), ppt2->error_message);
    for (int index_k_out=0; index_k_out < ppt2->k_out_size; ++index_k_out) {
      k_out[2*index_k_out] = ppt2->k1_out[index_k_out];
      k_out[2*index_k_out+1] = ppt2->k2_out[index_k_out];
    }

    /* Also add the largest and smallest between the k3 output values. We won't
    be producing output for these values, we add them just to avoid SONG
    complaining about k3 being out of bounds. */
    double k3_out_min = _HUGE_;
    double k3_out_max = 0;
    for (int index_k_out=0; index_k_out < ppt2->k_out_size; ++index_k_out) {
      k3_out_min = MIN (MIN (MIN (k3_out_min, ppt2->k3_out[index_k_out]), ppt2->k2_out[index_k_out]), ppt2->k2_out[index_k_out]);
      k3_out_max = MAX (MAX (MAX (k3_out_max, ppt2->k3_out[index_k_out]), ppt2->k2_out[index_k_out]), ppt2->k2_out[index_k_out]);
    }
    k_out[2*ppt2->k_out_size] = k3_out_min;
    k_out[2*ppt2->k_out_size+1] = k3_out_max;

    /* Merge ppt2->k with the k1 and k2 output points, sort the resulting array and
    remove the duplicates in it */
    class_call (merge_arrays_double (
                  ppt2->k,
                  ppt2->k_size,
                  k_out,
                  2*ppt2->k_out_size+2,
                  &(ppt2->k),
                  &(ppt2->k_size),
                  compare_doubles,
                  ppt2->error_message
                  ),
      ppt2->error_message,
      ppt2->error_message);

    /* Assign to each output k the corresponding index in ppt2->k */
    for (int i=0; i < 2*ppt2->k_out_size; ++i) {

      /* Find index in ppt2->k corresponding to the current output k */
      int index_k = 0;
      while (ppt2->k[index_k] != k_out[i])
        index_k++;

      class_test (index_k >= ppt2->k_size,
        ppt2->error_message,
        "index_k out of bounds: something went wrong while adding k output values");

      /* Store the index in the index_k_out arrays */
      int index_k_out = (i - i%2)/2;

      if (i%2 == 0)
        ppt2->index_k1_out[index_k_out] = index_k;
      else if (i%2 == 1)
        ppt2->index_k2_out[index_k_out] = index_k;

      /* Debug - Print the k->k_out correspondence */
      // printf ("k_out=%g[%d] -> k=%g[%d]\n",
      //   k_out[i], i, ppt2->k[index_k], index_k);
    }

    /* Check that index_k2_out is always smaller or equal than index_k1_out */
    for (int index_k_out=0; index_k_out < ppt2->k_out_size; ++index_k_out)
      class_test (ppt2->index_k2_out[index_k_out] > ppt2->index_k1_out[index_k_out],
        ppt2->error_message,
        "found index_k2_out=%d larger than index_k1_out=%d for #%d output",
        ppt2->index_k2_out[index_k_out], ppt2->index_k1_out[index_k_out], index_k_out);

    free (k_out);

  } // end of if k_out


  /* The user might also ask to output the perturbations via the indices
  index_k1 and index_k2 rather than with the exact values k1 and k2. We
  implement this feature here. */

  if (ppt2->k_index_out_size > 0) {

    for (int k_index_out=0; k_index_out < ppt2->k_index_out_size; ++k_index_out) {

      /* Copy the index_k1 and index_k2 provided directly by the user at the end
      of the index_k_out arrays */
      int first_index_k_out = ppt2->k_out_size;
      ppt2->index_k1_out[first_index_k_out + k_index_out] = ppt2->k1_index_out[k_index_out];
      ppt2->index_k2_out[first_index_k_out + k_index_out] = ppt2->k2_index_out[k_index_out];

      /* If the user gave a value of index_k1 or index_k2 which is not included in the sampling,
      set them to the largest possible value */

      if (ppt2->index_k1_out[first_index_k_out + k_index_out] >= ppt2->k_size) {

        fprintf (ppt2->k_out_files[first_index_k_out + k_index_out],
          "NOTE: The requested k1_index_out=%d is too large; we set it to the highest possible value: ppt2->k_size-1=%d.\n",
          ppt2->index_k1_out[first_index_k_out + k_index_out], ppt2->k_size-1);

        ppt2->index_k1_out[first_index_k_out + k_index_out] = ppt2->k_size-1;
      }

      if (ppt2->index_k2_out[first_index_k_out + k_index_out] >= ppt2->k_size) {

        fprintf (ppt2->k_out_files[first_index_k_out + k_index_out],
          "NOTE: The requested k2_index_out=%d is too large; we set it to the highest possible value: ppt2->k_size-1=%d.\n",
           ppt2->index_k2_out[first_index_k_out + k_index_out], ppt2->k_size-1);

        ppt2->index_k2_out[first_index_k_out + k_index_out] = ppt2->k_size-1;
      }

      /* Add the k1 and k2 values to the ppt2->k1_out and ppt2->k2_out arrays */
      ppt2->k1_out[first_index_k_out + k_index_out] = ppt2->k[ppt2->index_k1_out[first_index_k_out + k_index_out]];
      ppt2->k2_out[first_index_k_out + k_index_out] = ppt2->k[ppt2->index_k2_out[first_index_k_out + k_index_out]];

    }

    /* Update the number of k-triplet to output */
    ppt2->k_out_size += ppt2->k_index_out_size;

  }


  /* Debug - Print out the k-list */
  // printf ("# ~~~ k-sampling for k1 and k2 (size=%d) ~~~\n", ppt2->k_size);
  // for (int index_k=0; index_k < ppt2->k_size; ++index_k) {
  //   printf ("%17d %17.7g", index_k, ppt2->k[index_k]);
  //   for (int index_k_out=0; index_k_out < ppt2->k_out_size; ++index_k_out) {
  //     if (index_k==ppt2->index_k1_out[index_k_out]) printf ("\t(triplet #%d, k1) ", index_k_out);
  //     if (index_k==ppt2->index_k2_out[index_k_out]) printf ("\t(triplet #%d, k2) ", index_k_out);
  //   }
  //   printf ("\n");
  // }

  /* Check that the minimum and maximum values of ppt2->k are different. This
  test might fire if the user set a custom time sampling with two equal k-values */
  class_test ((ppt2->k[0]>=ppt2->k[ppt2->k_size-1]) && (ppt2->k_out_mode==_FALSE_),
    ppt2->error_message,
    "first and last value of ppt2->k coincide; maybe you set k_min_custom>=k_max_custom?");

  /* Check that the k array is strictly ascending */
  for (int index_k=1; index_k < ppt2->k_size; ++index_k)
    class_test (ppt2->k[index_k]<=ppt2->k[index_k-1],
      ppt2->error_message,
      "the k-sampling should be strictly ascending");



  // ====================================================================================
  // =                                      k3 grid                                     =
  // ====================================================================================

  /* The third wavemode is constrained to be the sum of the first two. As a result, every
  (k1,k2) configuration will yield a different grid for k3. There are as many k3 grids as
  possible (k1,k2) pairs. We store these grids in the array ppt2->k3[index_k1][index_k2]. */

  /* Initialize counter of total k-configurations */
  ppt2->count_k_configurations = 0;

  /* Allocate k1 level */
  int k1_size = ppt2->k_size;

  class_alloc (ppt2->k3, k1_size*sizeof(double **), ppt2->error_message);
  class_alloc (ppt2->k3_size, k1_size*sizeof(int *), ppt2->error_message);

  for (int index_k1=0; index_k1 < ppt2->k_size; ++index_k1) {

    double k1 = ppt2->k[index_k1];

    /* Allocate k2 level */
    int k2_size = index_k1 + 1;

    class_alloc (ppt2->k3[index_k1], k2_size*sizeof(double *), ppt2->error_message);
    class_alloc (ppt2->k3_size[index_k1], k2_size*sizeof(int), ppt2->error_message);

    for (int index_k2=0; index_k2 <= index_k1; ++index_k2) {

      /* Remember that, given our loop choice, k2 is always smaller than k1 */
      double k2 = ppt2->k[index_k2];

      /* The maximum and minimum values of k3 are obtained when the values of the cosine
      between k1 and k2 are respectively +1 and -1. In an numerical code, this is not exactly
      achievable. Relying on this would ultimately give rise to nan, for example when
      computing the Legendre polynomials or square roots. Hence, we shall never sample k3 too close
      to its exact minimum or maximum. */
      double k3_min = fabs(k1 - k2) + fabs(_MIN_K3_DISTANCE_);
      double k3_max = k1 + k2 - fabs(_MIN_K3_DISTANCE_);

      /* The differential system dies when k1=k2 and k3 is very small. These configurations
      are irrelevant, so we set a minimum ratio between k1=k2 and k3. TODO: remove? */
      k3_min = MAX (k3_min, (k1+k2)/_MIN_K3_RATIO_);

      /* We take k3 in the same range as k1 and k2. Comment it out if you prefer a range that
      goes all the way to the limits of the triangular condition. If you do so, remember to
      double pbs2->xx_max in input2.c */
      k3_min = MAX (k3_min, ppt2->k[0]);
      k3_max = MIN (k3_max, ppt2->k[ppt2->k_size-1]);

      /* Check that the chosen k3_min and k3_max make sense */
      class_test ((k3_min >= k3_max) && (ppt2->k_out_mode == _FALSE_),
        ppt2->error_message,
        "found k3_min=%g>k3_max=%g for k1(%d)=%g and k2(%d)=%g",
        k3_min, k3_max, index_k1, k1, index_k2, k2);

      /* Shortcuts for ppt2->k3[index_k1][index_k2] and ppt2->k3_size[index_k1][index_k2] */
      double * k3_grid;
      int k3_size;


      // -------------------------------------------------------------------------------
      // -                            Lin/log sampling for k3                          -
      // -------------------------------------------------------------------------------

      /* Adopt a simple sampling with a fixed number of points for each (k1,k2) configuration.
      This is not efficient, as low values of k1 and k2 do not need a sampling as good as the
      one needed for high values */

      if ((ppt2->k3_sampling == lin_k3_sampling) || (ppt2->k3_sampling == log_k3_sampling)) {

        /* The size of the k3 array is the same for every (k1,k2) configuration, and is read
        from the precision structure */
        k3_size = ppr2->k3_size;
        class_alloc (k3_grid, ppr2->k3_size*sizeof(double), ppt2->error_message);

        /* Build the grids */
        if (ppt2->k3_sampling == log_k3_sampling) {
          log_space (k3_grid, k3_min, k3_max, ppr2->k3_size);
        }
        else if (ppt2->k3_sampling == lin_k3_sampling) {
          lin_space (k3_grid, k3_min, k3_max, ppr2->k3_size);
        }

      } // end of lin/log sampling


      // -------------------------------------------------------------------------------
      // -                       Linear sampling for the angle                         -
      // -------------------------------------------------------------------------------

      /* Sample linearly the angle between k1 and k3 */

      else if (ppt2->k3_sampling == theta12_k3_sampling) {

        /* The size of the k3 array is the same for every (k1,k2) configuration, and is read
        from the precision structure */
        k3_size = ppr2->k3_size;
        class_alloc (k3_grid, ppr2->k3_size*sizeof(double), ppt2->error_message);

        double cosk1k2_min = (k3_min*k3_min - k1*k1 - k2*k2)/(2*k1*k2);
        double cosk1k2_max = (k3_max*k3_max - k1*k1 - k2*k2)/(2*k1*k2);

        class_test (fabs(cosk1k2_min)>1, ppt2->error_message, "stop to prevent nans");
        class_test (fabs(cosk1k2_max)>1, ppt2->error_message, "stop to prevent nans");

        double theta12_min = acos(cosk1k2_min);
        double theta12_max = acos(cosk1k2_max);
        double theta_step = (theta12_max - theta12_min)/(k3_size-1);

        if (index_k1==index_k2) {
          double angle_factor = 5;
          theta12_min += theta_step/angle_factor;
          theta_step = (theta12_max - theta12_min)/(k3_size-1);
        }

        for (int index_k3=0; index_k3 < k3_size; ++index_k3) {
          double cosk1k2 = cos (theta12_min + index_k3*theta_step);
          k3_grid[index_k3] = sqrt (k1*k1 + k2*k2 + 2*cosk1k2*k1*k2);
        }

      } // end of theta13 sampling


      // -------------------------------------------------------------------------------
      // -                              Symmetric sampling                             -
      // -------------------------------------------------------------------------------

      /* Use for k3 exactly the same sampling as for k1 an k2. See documentation
      for sym_k3_sampling in perturbations2.h for more details. Note that the
      output produced in this way might be unprecise if the requested k-values
      are isolated from the remaining k-sampling, because the symmetric sampling
      relies on interpolation of the quadratic sources over k. */

      else if (ppt2->k3_sampling == sym_k3_sampling) {

        /* The size of the k3 array is the same for every (k1,k2) configuration*/
        k3_size = ppt2->k_size;

        class_alloc (k3_grid, k3_size*sizeof(double), ppt2->error_message);

        /* Build the grids */
        for (int index_k3=0; index_k3 < ppt2->k_size; ++index_k3)
          k3_grid[index_k3] = ppt2->k[index_k3];

        /* There is no triangular condition in the symmetric sampling */
        k3_min = ppt2->k[0];
        k3_max = ppt2->k[ppt2->k_size-1];

      } // end of sym sampling


      // -------------------------------------------------------------------------------
      // -                             Smart sampling for k3                           -
      // -------------------------------------------------------------------------------

      /* Use for k3 the same grid we used for k1 and k2, making sure that at least
      ppr2->k3_size_min points are included for each (k1,k2) configuration. */

      else if (ppt2->k3_sampling == smart_k3_sampling) {

        /* Test that k3_min and k3_max are within bounds */
        class_test ((k3_min < ppt2->k[0]) || (k3_min > ppt2->k[ppt2->k_size-1]),
          ppt2->error_message,
          "k3_min=%g out of bounds", k3_min);

        class_test ((k3_max < ppt2->k[0]) || (k3_max > ppt2->k[ppt2->k_size-1]),
          ppt2->error_message,
          "k3_max=%g out of bounds", k3_max);

        /* Find the minimum allowed index of k3 inside ppt2->k. If k3_min is smaller than
        the smallest element in ppt2->k, take the latter. */
        int index_k3_min = 0;
        while (ppt2->k[index_k3_min] < k3_min)
          ++index_k3_min;

        /* Find the maximum allowed index of k. If k_max is larger than the largest element
        in ppt2->k, take the latter. */
        int index_k3_max = ppt2->k_size-1;
        while (ppt2->k[index_k3_max] > k3_max)
          --index_k3_max;

        /* Number of points in ppt2->k between k3_min and k3_max */
        k3_size = index_k3_max - index_k3_min + 1;

        /* If k3_min and k3_max are too close, the number of points between them will be zero.
        However, since both of them are within the bounds of ppt2->k, k3_size can be zero but
        not negative */
        class_test (k3_size < 0,
          ppt2->error_message,
          "something went wrong in bracketing k3_min and k3_max in ppt2->k");

        /* Copy the bracketed points from ppt2->k to k3 */
        class_alloc (k3_grid, k3_size*sizeof(double), ppt2->error_message);
        for (int index_k3=0; index_k3 < k3_size; ++index_k3)
          k3_grid[index_k3] = ppt2->k[index_k3_min+index_k3];

        /* Add by hand the points corresponding to k3_min and k3_max, because they are
        most probably not part of ppt2->k. The function merge_arrays_double() will
        take care of potential duplicates. */

        double k3_bounds[2] = {k3_min, k3_max};

        class_call (merge_arrays_double (
                      k3_grid,
                      k3_size,
                      &k3_bounds[0],
                      2,
                      &(k3_grid),
                      &(k3_size),
                      compare_doubles,
                      ppt2->error_message
                      ),
          ppt2->error_message,
          ppt2->error_message);

        /* We choose the grid to have at least ppr2->k3_size_min values for every (k1,k2).
        If this is not possible using the standard ppt2->k grid, include as many linearly
        sampled points between k3_min and k3_max. */

        if (k3_size < ppr2->k3_size_min) {

          /* We want a total of k3_size_min points, therefore the number of points to add
          should be k3_size_min-k3_size. We add two more points because k3_min and k3_max
          are already included in k3_grid and therefore will not be considered. */
          int n_extra_points = ppr2->k3_size_min - k3_size + 2;

          /* Fill the k3 array with linearly-spaced points between k3_min and k3_max */
          double * extra_points;
          class_alloc (extra_points, n_extra_points*sizeof(double), ppt2->error_message);
          lin_space (extra_points, k3_min, k3_max, n_extra_points);

          /* Include the linearly spaced points in the k3 grid */
          class_call (merge_arrays_double (
                        k3_grid,
                        k3_size,
                        extra_points,
                        n_extra_points,
                        &(k3_grid),
                        &(k3_size),
                        compare_doubles,
                        ppt2->error_message
                        ),
            ppt2->error_message,
            ppt2->error_message);

          free (extra_points);

        }

      } // end of smart sampling


      // -------------------------------------------------------------------------------
      // -                             Add k3 output points                            -
      // -------------------------------------------------------------------------------

      if (ppt2->k_out_size > 0) {

        /* If SONG is running in k_out_mode, we ignore the k3 grid computed above and
        start over. If an output time or redshift is requested, we do keep the k3 grid,
        lest the output files have only one entry. */

        if ((ppt2->k_out_mode == _TRUE_) && ((ppt2->tau_out_size+ppt2->z_out_size) <= 0))
          k3_size = 0;

        /* Add the output values to the k3 sampling. These values are contained in
        ppt2->k3_out and satisfy the triangular condition.  */

        for (int index_k_out=0; index_k_out < ppt2->k_out_size; ++index_k_out) {

          /* Add the k3 output points to the k3 grid, but only for those k1 and k2 that are
          inside k1_out and k2_out. Given a (k1,k2) pair, we will be adding the k3_out values
          one at a time. If two triplets have the same (k1,k2) pair, this will mess up the
          assignation of the index_k3_out indices, because each addition will shift the grid
          and with it the indices assigned previously. We prevent this from happening by
          enforcing in input2.c that the k3 values are given in ascending order. */

          if ((index_k1==ppt2->index_k1_out[index_k_out]) && (index_k2==ppt2->index_k2_out[index_k_out])) {

            /* If this k3 output value was added by specifying an exact value, then we do need to
            add this value to the k3 grid */

            if (index_k_out < (ppt2->k_out_size - ppt2->k_index_out_size)) {

              /* k3 value to be added to the grid */
              double k3 = ppt2->k3_out[index_k_out];

              /* Alert the user if k3 is either larger than k3_max or smaller than k3_min, which
              are the limits that SONG would would have imposed on k3 if no output were requested */
              if (k3 < k3_min)
                fprintf (ppt2->k_out_files[index_k_out],
                  "# NOTE: the output value chosen for this file, k3=%.17f, was smaller than k3_min=%.17f (|k1-k2|=%.17f)\n",
                  k3, k3_min, fabs(k1-k2));

              if (k3 > k3_max)
                fprintf (ppt2->k_out_files[index_k_out],
                  "# NOTE: the output value chosen for this file, k3=%.17f, was larger than k3_max=%.17f (|k1-k2|=%.17f)\n",
                  k3, k3_max, k1+k2);

              /* Add the considered k3 point to the k3 grid */
              class_call (merge_arrays_double (
                            k3_grid,
                            k3_size,
                            &k3,
                            1,
                            &(k3_grid),
                            &(k3_size),
                            compare_doubles,
                            ppt2->error_message
                            ),
                ppt2->error_message,
                ppt2->error_message);

              /* Find the index corresponding to k3 in the new k3 grid */
              int index_k3 = k3_size - 1;
              while ((index_k3 >= 0) && (k3_grid[index_k3] != k3))
                index_k3--;

              class_test (index_k3<0,
                ppt2->error_message,
                "k3_out=%g not found for index_k1=%d, index_k2=%d; bug in merge_arrays_double?",
                index_k1,index_k2);

              ppt2->index_k3_out[index_k_out] = index_k3;

            } // end of if


            /* If this k3 output value was added by specifying an index rather than an exact
            value, then the value is already in the k3 grid and we do not need to add it */

            else {

              /* Copy the index_k3 provided directly by the user */
              int first_index_k_out = ppt2->k_out_size - ppt2->k_index_out_size;
              ppt2->index_k3_out[index_k_out] = ppt2->k3_index_out[index_k_out-first_index_k_out];

              /* If the user gave a value of index_k3 which is not included in the sampling,
              set index_k3 to the largest possible value */

              if (ppt2->index_k3_out[index_k_out] >= k3_size) {

                fprintf (ppt2->k_out_files[index_k_out],
                  "NOTE: The requested k3_index_out=%d is larger than the size of the k3 grid for the requested\
 index_k1=%d and index_k2=%d; we have set index_k3 to the highest possible value, k3_size-1=%d.\n",
                  ppt2->index_k3_out[index_k_out], index_k1, index_k2, k3_size-1);

                ppt2->index_k3_out[index_k_out] = k3_size-1;
              }

              /* Add the k3 value to the ppt2->k3_out array */
              ppt2->k3_out[index_k_out] = k3_grid[ppt2->index_k3_out[index_k_out]];

            } // end of if(index or value)

            /* Debug - Print the k3 grid for all (k1,k2) configurations requested as output */
            // fprintf (stderr, "k1[%d]=%g, k2[%d]=%g, k3_size=%d, k3_min=%g, k3_max=%g\n",
            //   index_k1, k1, index_k2, k2, k3_size, k3_min, k3_max);
            //
            // for (int index_k3=0; index_k3 < k3_size; ++index_k3) {
            //   double k3 = k3_grid[index_k3];
            //   fprintf(stderr, "%3d %12g ", index_k3, k3);
            //   if (index_k3 == ppt2->index_k3_out[index_k_out])
            //     fprintf (stderr, "\t(triplet #%d) ", index_k_out);
            //   fprintf (stderr, "\n");
            // }
            // fprintf (stderr, "\n\n");

          } // end of if (k1==k1_out && k2==k2_out)

        } // end of for k_out

      } // end of if k_out


      // -------------------------------------------------------------------------------
      // -                                 Update grid                                 -
      // -------------------------------------------------------------------------------

      /* Convert the shortcuts to the real stuff */
      ppt2->k3[index_k1][index_k2] = k3_grid;
      ppt2->k3_size[index_k1][index_k2] = k3_size;

      /* Update counter of k-configurations */
      ppt2->count_k_configurations += k3_size;

      /* Debug - Print out the k3 list for a special configuration */
      // if ((index_k1==9) && (index_k2==5)) {
      //
      //   fprintf (stderr, "k1[%d]=%.17f, k2[%d]=%.17f, k3_size=%d, k3_min=%.17f, k3_max=%.17f\n",
      //     index_k1, k1, index_k2, k2, k3_size, k3_min, k3_max);
      //
      //   for (int index_k3=0; index_k3 < k3_size; ++index_k3)
      //     fprintf(stderr, "%d %.17f\n", index_k3, k3_grid[index_k3]);
      //
      //   fprintf (stderr, "\n\n");
      // }

    } // end of for (index_k2)

  } // end of for (index_k1)


  /* Set the minimum and maximum k-values that will be fed to the differential system.
  This corresponds to the minimum and maximum k-values ever used in SONG. */

  ppt2->k_min = _HUGE_;
  ppt2->k_max = 0;

  double kt_min = _HUGE_;
  double kt_max = 0;

  for (int index_k1=0; index_k1 < ppt2->k_size; ++index_k1) {

    double k1 = ppt2->k[index_k1];

    for (int index_k2=0; index_k2 <= index_k1; ++index_k2) {

      double k2 = ppt2->k[index_k2];

      for (int index_k3=0; index_k3 < ppt2->k3_size[index_k1][index_k2]; ++index_k3) {

        double k3 = ppt2->k3[index_k1][index_k2][index_k3];

        ppt2->k_min = MIN (MIN (MIN (ppt2->k_min, k3), k2), k1);
        ppt2->k_max = MAX (MAX (MAX (ppt2->k_max, k3), k2), k1);

        /* In case of symmetric sampling, the differential system evolves a different
        set of wavemodes. */

        if (ppt2->k3_sampling == sym_k3_sampling) {

          double K[4] = {0, k1, k2, k3};
          double KT[4];

          class_call (symmetric_sampling (K, KT, ppt2->error_message),
            ppt2->error_message, ppt2->error_message);

          kt_min = MIN (MIN (MIN (kt_min, KT[3]), KT[2]), KT[1]);
          kt_max = MAX (MAX (MAX (kt_max, KT[3]), KT[2]), KT[1]);

        }
      }
    }
  }

  /* The symmetric sampling should have the same maximum and minimum k, because of
  the structure of the transformation. Here we check that this is indeed the case. */

  if (ppt2->k3_sampling == sym_k3_sampling) {

    class_test (fabs(1-kt_min/ppt2->k_min) > _SMALL_,
      ppt2->error_message,
      "inconsistency in kt transformation");

    class_test (fabs(1-kt_max/ppt2->k_max) > _SMALL_,
      ppt2->error_message,
      "inconsistency in kt transformation");
  }

  class_test ((ppt2->k_min==0) || (ppt2->k_max==0),
    ppt2->error_message,
    "found vanishing value in either k_min (%g) or k_max (%g)",
    ppt2->k_min, ppt2->k_max);

  /* Debug - Print the k-output indices */
  // for (int index_k_out=0; index_k_out < ppt2->k_out_size; ++index_k_out) {
  //   printf ("index_k_out=%3d:  (%3d,%3d,%3d)\n",
  //     index_k_out,
  //     ppt2->index_k1_out[index_k_out],
  //     ppt2->index_k2_out[index_k_out],
  //     ppt2->index_k3_out[index_k_out]);
  // }



  // ====================================================================================
  // =                                    First-order k-grid                            =
  // ====================================================================================

  /* The 1st-order system must be solved for the k-values needed at second-order. Here we
  set all the relevant k-arrays in the first-order perturbation module to match the
  k-sampling contained in ppt2->k. */

  int md_size = 1; /* only scalars at first supported so far */
  int index_md_scalars = 0; /* only scalars at first supported so far */

  class_alloc (ppt->k_size, md_size*sizeof(int *), ppt2->error_message);
  class_alloc (ppt->k_size_cl, md_size*sizeof(int *), ppt2->error_message);
  class_alloc (ppt->k_size_cmb, md_size*sizeof(int *), ppt2->error_message);
  class_alloc (ppt->k, md_size*sizeof(double *), ppt2->error_message);

  for (int index_md=0; index_md < md_size; ++index_md) {

    ppt->k_size[index_md] = ppt2->k_size;
    ppt->k_size_cl[index_md] = ppt2->k_size;
    ppt->k_size_cmb[index_md] = ppt2->k_size;

    class_alloc (ppt->k[index_md], ppt->k_size[index_md]*sizeof(double), ppt2->error_message);
    for (int index_k=0; index_k < ppt2->k_size; ++index_k)
      ppt->k[index_md][index_k] = ppt2->k[index_k];
  }

  /* Determine ppt->k_min and ppt->k_max. The following block is copied from perturbations.c */
  ppt->k_min = _HUGE_;
  ppt->k_max = 0;
  if (ppt->has_scalars == _TRUE_) {
    ppt->k_min = MIN(ppt->k_min,ppt->k[index_md_scalars][0]);
    ppt->k_max = MAX(ppt->k_max,ppt->k[index_md_scalars][ppt->k_size[index_md_scalars]-1]);
  }

  /* Debug - Print first-order k-sampling */
  // for (int index_k=0; index_k < ppt->k_size[index_md_scalars]; ++index_k) {
  //   printf ("%5d %10g\n", index_k, ppt->k[index_md_scalars][index_k]);
  // }


  /* Make CLASS output the first-order perturbations in the same k where SONG outputs
  the second-order ones */

  if ((ppt2->k_out_size > 0) && (ppt2->output_class_perturbations == _TRUE_)) {

    /* Allocate and fill the array with the output values for k1, k2 and k3 */
    double * k_out_class;
    class_alloc (k_out_class, 2*ppt2->k_out_size*sizeof(double), ppt2->error_message);
    for (int index_k_out=0; index_k_out < ppt2->k_out_size; ++index_k_out) {
      k_out_class[2*index_k_out] = ppt2->k1_out[index_k_out];
      k_out_class[2*index_k_out + 1] = ppt2->k2_out[index_k_out];
    }

    /* Sort k_out_class and remove duplicates from it */
    class_call (merge_arrays_double (
                  k_out_class,
                  2*ppt2->k_out_size,
                  NULL,
                  0,
                  &(k_out_class),
                  &(ppt->k_output_values_num),
                  compare_doubles,
                  ppt2->error_message
                  ),
      ppt2->error_message,
      ppt2->error_message);

    class_test (ppt->k_output_values_num > _MAX_NUMBER_OF_K_FILES_,
     ppt2->error_message,
     "reached limit of k output files - increase _MAX_NUMBER_OF_K_FILES_ in perturbations.h");

    /* Copy k_out_class into ppt->k_output_values */
    for (int index_k_out=0; index_k_out < ppt->k_output_values_num; ++index_k_out)
      ppt->k_output_values[index_k_out] = k_out_class[index_k_out];

    /* Allocate and fill array with indices of output values */
    class_alloc(ppt->index_k_output_values, sizeof(double)*ppt->k_output_values_num, ppt2->error_message);

    for (int index_k_out=0; index_k_out < ppt->k_output_values_num; ++index_k_out) {

      int index_md_scalars = 0;

      /* Find index in ppt->k corresponding to the current output k */
      int index_k = 0;
      while (ppt->k[index_md_scalars][index_k] != ppt->k_output_values[index_k_out])
        index_k++;

      class_test (index_k >= ppt->k_size[index_md_scalars],
        ppt2->error_message,
        "index_k=%d out of bounds: something went wrong while adding k output values", index_k);

      ppt->index_k_output_values[index_k_out] = index_k;

    }

    /* Debug - Print the first-order output values and indices */
    // for (int index_k_out=0; index_k_out < ppt->k_output_values_num; ++index_k_out)
    //   printf ("%4d %12d %12g\n", index_k_out, ppt->index_k_output_values[index_k_out], ppt->k_output_values[index_k_out]);

    free (k_out_class);

  } // end of if k_out

  /* Debug - Print out the first-order k-list */
  // printf ("# ~~~ first-order k-sampling (size=%d) ~~~\n", ppt->k_size[index_md_scalars]);
  // for (int index_k=0; index_k < ppt->k_size[index_md_scalars]; ++index_k) {
  //   printf ("%17d %17.7g", index_k, ppt->k[index_md_scalars][index_k]);
  //   for (int index_k_out=0; index_k_out < ppt->k_output_values_num; ++index_k_out)
  //     if (index_k==ppt->index_k_output_values[index_k_out])
  //       printf ("\t(output #%d) ", index_k_out);
  //   printf ("\n");
  // }


  return _SUCCESS_;

}
// End of galbispectra2_get_k_lists















/*==================================================================
-------------------- Define Time Sampling --------------------------
====================================================================*/

/**
 *
 * Define the sampling in conformal time tau for the first order line of sight sources.
 *
 * In detail, this function does:
 *
 *  -# Define the time sampling for the line of sight sources in ppt2->tau_sampling.
 *     For the CMB, these are the sources that will be integrated over time in the
 *     transfer2.c module to obtain today's value of the transfer functions.
*/

int galbispectra2_timesampling_for_sources (
             struct precision * ppr,
             struct precision2 * ppr2,
             struct background * pba,
             struct thermo * pth,
             struct perturbs * ppt,
             struct perturbs2 * ppt2
             )
{

  /* Temporary arrays used to store background and thermodynamics quantities, defined on the
      stack. Only usable inside of galbispectra2_timesampling_for_sources()*/
  double *pvecback, *pvecthermo;
  class_alloc (pvecback, pba->bg_size*sizeof(double), ppt2->error_message);
  class_alloc (pvecthermo, pth->th_size*sizeof(double), ppt2->error_message);
  int dump;


  // =====================================================================================
  // =                              Custom sources sampling                              =
  // =====================================================================================

  /* The user can specify a time sampling for the sources via the parameter file,
  by providing a start time, and end time and a preference for the sampling method
  (either linear or logarithmic). */

  if (ppt2->has_custom_timesampling == _TRUE_) {

    ppt2->tau_size = ppt2->custom_tau_size;
    class_alloc (ppt2->tau_sampling, ppt2->tau_size*sizeof(double), ppt2->error_message);

    /* If the user set the custom end-time to 0, we assume that they want to compute
    the sources all the way to today */
    if (ppt2->custom_tau_end == 0)
      ppt2->custom_tau_end = pba->conformal_age;

    /* Linear sampling */
    if (ppt2->custom_tau_mode == lin_tau_sampling) {
      lin_space(ppt2->tau_sampling, ppt2->custom_tau_ini, ppt2->custom_tau_end, ppt2->tau_size);
    }
    /* Logarithmic sampling */
    else if (ppt2->custom_tau_mode == log_tau_sampling) {
      log_space(ppt2->tau_sampling, ppt2->custom_tau_ini, ppt2->custom_tau_end, ppt2->tau_size);
    }
  }



  // =====================================================================================
  // =                            Standard sources sampling                              =
  // =====================================================================================

  /* Determine the time sampling of the line of sight sources using the same algorithm
  as in standard CLASS, whereby the sampling is denser where the perturbations evolve
  faster, based on the visibility function and on the Hubble rate. For details, please
  refer to CLASS function perturb_timesampling_for_sources(). With respect to that
  function, we use the same parameter to determine the starting time
  (ppr->start_sources_at_tau_c_over_tau_h) but a different parameter to determine the
  sampling frequency (ppr2->perturb_sampling_stepsize_song rather than
  ppr2->perturb_sampling_stepsize_song). */

  else {

    // -----------------------------------------------------------------------------
    // -                            Find first point                               -
    // -----------------------------------------------------------------------------

    double tau_ini;

    /* Using bisection, search the time such that the ratio between the Hubble
    time-scale tau_h = 1/aH and the Compton time-scale 1/kappa_dot is equal to
    ppr->start_sources_at_tau_c_over_tau_h. Usually, this parameter is about
    0.01, which means that we start sampling the sources in the tight coupling
    regime (tau_c<<tau_h), where the visibility function is still small.  */

    if (ppt2->has_cmb == _TRUE_) {

      double tau_lower = pth->tau_ini;

      class_call (background_at_tau(pba,
                   tau_lower,
                   pba->short_info,
                   pba->inter_normal,
                   &dump,
                   pvecback),
        pba->error_message,
        ppt2->error_message);

      double a = pvecback[pba->index_bg_a];
      double Hc = a*pvecback[pba->index_bg_H];

      class_call (thermodynamics_at_z(pba,
                    pth,
                    1/a-1,
                    pth->inter_normal,
                    &dump,
                    pvecback,
                    pvecthermo),
        pth->error_message,
        ppt2->error_message);

      double kappa_dot = pvecthermo[pth->index_th_dkappa];

      class_test (
        Hc/kappa_dot > ppr->start_sources_at_tau_c_over_tau_h,
        ppt2->error_message,
        "your choice of initial time for computing sources is inappropriate: it corresponds\
 to an earlier time than the one at which the integration of thermodynamical variables\
 started (tau=%g). You should increase either 'start_sources_at_tau_c_over_tau_h' or\
 'recfast_z_initial'\n",
        tau_lower);

      /* The upper limit is when the visibility function peaks */
      double tau_upper = pth->tau_rec;

      class_call (background_at_tau(pba,
                    tau_upper,
                    pba->short_info,
                    pba->inter_normal,
                    &dump,
                    pvecback),
        pba->error_message,
        ppt2->error_message);

      a = pvecback[pba->index_bg_a];
      Hc = a*pvecback[pba->index_bg_H];

      class_call (thermodynamics_at_z(pba,
                    pth,
                    1/a-1,
                    pth->inter_normal,
                    &dump,
                    pvecback,
                    pvecthermo),
        pth->error_message,
        ppt2->error_message);

      kappa_dot = pvecthermo[pth->index_th_dkappa];

      class_test (Hc/kappa_dot < ppr->start_sources_at_tau_c_over_tau_h,
        ppt2->error_message,
        "your choice of initial time for computing sources is inappropriate: it corresponds\
 to a time after recombination. You should decrease 'start_sources_at_tau_c_over_tau_h'\n");

      double tau_mid = 0.5*(tau_lower + tau_upper);

      while (tau_upper - tau_lower > ppr->tol_tau_approx) {

        class_call (background_at_tau(pba,
                      tau_mid,
                      pba->short_info,
                      pba->inter_normal,
                      &dump,
                      pvecback),
          pba->error_message,
          ppt2->error_message);

        a = pvecback[pba->index_bg_a];
        Hc = a*pvecback[pba->index_bg_H];

        class_call (thermodynamics_at_z(pba,
                      pth,
                      1/a-1,
                      pth->inter_normal,
                      &dump,
                      pvecback,
                      pvecthermo),
          pth->error_message,
          ppt2->error_message);

        kappa_dot = pvecthermo[pth->index_th_dkappa];

        if (Hc/kappa_dot > ppr->start_sources_at_tau_c_over_tau_h)
          tau_upper = tau_mid;
        else
          tau_lower = tau_mid;

        tau_mid = 0.5*(tau_lower + tau_upper);
      }

      tau_ini = tau_mid;
    }

    /* If the CMB is not requested, just start sampling the perturbations at
    recombination time (copied from CLASS) */

    else {

      tau_ini = pth->tau_rec;

    } // if(has_cmb)


    // -----------------------------------------------------------------------------
    // -                           Determine time-grid                             -
    // -----------------------------------------------------------------------------

    /* Since we do not know yet how many points to include in the time sampling,
    we first allocate ppt2->tau_sampling with a very large value */
    class_alloc (ppt2->tau_sampling, TAU_SIZE_MAX*sizeof(double), ppt2->error_message);

    /* Set the value of the first point in the time sampling */
    int index_tau = 0;
    ppt2->tau_sampling[index_tau] = tau_ini;
    double tau = tau_ini;
    index_tau++;

    /* Add points to the time-sampling until we reach today */
    while (tau < pba->conformal_age) {

      /* The next sampling point is determined by the lowest of two timescales: the time
      variation of the visibility function and the acceleration parameter. Schematically:

        next = previous + ppr2->perturb_sampling_stepsize_song * timescale_source

      where:

        timescale_source = 1 / (1/timescale_source1 + 1/timescale_source2)
        timescale_source1 = g/g_dot
        timescale_source2 = 1/sqrt |2*a_dot_dot/a - Hc^2| */

      class_call (background_at_tau(pba,
                    tau,
                    pba->short_info,
                    pba->inter_normal,
                    &dump,
                    pvecback),
        pba->error_message,
        ppt2->error_message);

      double a = pvecback[pba->index_bg_a];
      double Hc = a*pvecback[pba->index_bg_H];
      double H_prime = pvecback[pba->index_bg_H_prime];

      class_call (thermodynamics_at_z(pba,
                    pth,
                    1./pvecback[pba->index_bg_a]-1.,  /* redshift z=1/a-1 */
                    pth->inter_normal,
                    &dump,
                    pvecback,
                    pvecthermo),
        pth->error_message,
        ppt2->error_message);

      double kappa_dot = pvecthermo[pth->index_th_dkappa];

      /* If the CMB is requested, the time sampling needs to be denser at recombination
      and when the ISW effect is important */

      double timescale_source;

      if (ppt2->has_cmb == _TRUE_) {

        /* Variation rate of thermodynamics variables */
        double rate_thermo = pvecthermo[pth->index_th_rate];

        /* Variation rate of metric due to late ISW effect (important at late times) */
        double a_primeprime_over_a = a*H_prime + 2*Hc*Hc;
        double rate_isw_squared = fabs (2*a_primeprime_over_a - Hc*Hc);

        /* Add points to the late-time part of the time sampling, in view of the
        time integration in transfer2.c */
        rate_isw_squared *= pow (ppr2->perturb_sampling_late_time_boost, 2);

        /* Compute rate */
        timescale_source = 1/sqrt(rate_thermo*rate_thermo + rate_isw_squared);

        /* Debug - Print ratio between timescales */
        // printf ("%12g %12g %12g %12g%12g\n", tau, rate_thermo/sqrt(rate_isw_squared),
        //   rate_thermo, sqrt(rate_isw_squared), timescale_source);

      }

      /* If the CMB is not requested, we use as sampling frequency the conformal Hubble
      rate aH, which is the natural time scale of the differential system. Since 1/aH is
      proportional to the conformal time tau, this sampling roughly corresponds to a
      logarithmic sampling in tau. */

      else {

        /* Variation rate given by Hubble time */
        timescale_source = 1/Hc;

      }

      /* Update the time-sampling array with the new value */
      double step = ppr2->perturb_sampling_stepsize_song * timescale_source;

      class_test(
        fabs(step/tau) < ppr->smallest_allowed_variation,
        ppt2->error_message,
        "integration step =%e < machine precision: leads to infinite loop",
        ppr2->perturb_sampling_stepsize_song*timescale_source);

      tau = tau + step;
      ppt2->tau_sampling[index_tau] = tau;
      index_tau++;

      class_test ((index_tau+1) > TAU_SIZE_MAX,
        ppt2->error_message,
        "ppt2->tau_sampling size is too large; check the perturb_sampling_stepsize_song parameter");
    }

    /* Total number of time steps */
    ppt2->tau_size = index_tau;

    /* Last sampling point = exactly today */
    ppt2->tau_sampling[ppt2->tau_size-1] = pba->conformal_age;

    /* Free the excess memory we have allocated in ppt2->tau_sampling */
    class_realloc(ppt2->tau_sampling,
                  ppt2->tau_sampling,
                  ppt2->tau_size*sizeof(double),
                  ppt2->error_message);

  } // if has_custom_timesampling


  // ====================================================================================
  // =                                Add output points                                 =
  // ====================================================================================

  /* The user might have asked to output the perturbations at specific tau times
  using the tau_out parameter. Here we add these tau-values to the list of computed
  tau in SONG, that is, to ppt2->tau_sampling. */

  /* Convert the requested redshifts in z_out to conformal times, and add them to ppt2->tau_out */

  for (int index_z=0; index_z < ppt2->z_out_size; ++index_z) {

    if (ppt2->z_out[index_z] < 0) {

      for (int index_k_out=0; index_k_out < ppt2->k_out_size; ++index_k_out)
        fprintf (ppt2->tau_out_files[index_k_out][ppt2->tau_out_size+index_z],
          "# NOTE: z_out was negative for this file, so we set it to z=0\n");

      ppt2->z_out[index_z] = 0;

    }

    class_call (background_tau_of_z (
                  pba,
                  ppt2->z_out[index_z],
                  &ppt2->tau_out[ppt2->tau_out_size+index_z]),
      ppt2->error_message,
      ppt2->error_message);

  }

  ppt2->tau_out_size += ppt2->z_out_size;


  /* Merge ppt2->tau_out with ppt2->tau_sampling */

  if (ppt2->tau_out_size > 0) {

    /* Check that the requested times can be computed by SONG */
    for (int index_tau_out=0; index_tau_out < ppt2->tau_out_size; ++index_tau_out) {

      class_test (ppt2->tau_out[index_tau_out] < pth->tau_ini,
        ppt2->error_message,
        "choose tau_out to be larger than %g (you choose %g)",
          pth->tau_ini, ppt2->tau_out[index_tau_out]);


      /* If the requested time is too large, set it to the largest available time. Comment the if
      block to add the requested time anyway. */

      ppt2->tau_out_was_reduced[index_tau_out] = _FALSE_;

      sprintf (ppt2->tau_out_reduction_message,
        "NOTE: the requested time for this file was too large, or redshift was too low; we set tau_out to the largest available value in our time sampling (%g).",
        ppt2->tau_sampling[ppt2->tau_size-1]);

      if (ppt2->tau_out[index_tau_out] > ppt2->tau_sampling[ppt2->tau_size-1]) {

        ppt2->tau_out_was_reduced[index_tau_out] = _TRUE_;

        ppt2->tau_out[index_tau_out] = ppt2->tau_sampling[ppt2->tau_size-1];

      }
    }

    /* Merge ppt2->tau_sampling with the tau output points, sort the resulting array
    and remove the duplicates in it */
    class_call (merge_arrays_double (
                  ppt2->tau_sampling,
                  ppt2->tau_size,
                  ppt2->tau_out,
                  ppt2->tau_out_size,
                  &(ppt2->tau_sampling),
                  &(ppt2->tau_size),
                  compare_doubles,
                  ppt2->error_message
                  ),
      ppt2->error_message,
      ppt2->error_message);

    /* Assign to each output tau the corresponding index in ppt2->tau_sampling */
    for (int index_tau_out=0; index_tau_out < ppt2->tau_out_size; ++index_tau_out) {

      int index_tau = 0;
      while (ppt2->tau_sampling[index_tau] != ppt2->tau_out[index_tau_out])
        index_tau++;

      class_test (index_tau >= ppt2->tau_size,
        ppt2->error_message,
        "index_tau out of bounds: something went wrong while adding tau output values");

      ppt2->index_tau_out[index_tau_out] = index_tau;

      /* Debug - Print the tau->tau_out correspondence */
      // printf ("tau_out=%g[%d] -> tau=%g[%d]\n",
      //   ppt2->tau_out[index_tau_out], index_tau_out,
      //   ppt2->tau_sampling[index_tau], index_tau);
    }

  } // end of if tau_out


  /* Debug - Print time sampling */
  // fprintf (stderr, "# ~~~ tau-sampling for the source function ~~~\n");
  // for (int index_tau=0; index_tau < ppt2->tau_size; ++index_tau) {
  //   fprintf (stderr, "%12d %16g", index_tau, ppt2->tau_sampling[index_tau]);
  //   for (int index_tau_out=0; index_tau_out < ppt2->tau_out_size; ++index_tau_out) {
  //     if (index_tau == ppt2->index_tau_out[index_tau_out])
  //       if (index_tau_out < (ppt2->tau_out_size-ppt2->z_out_size))
  //         fprintf (stderr, "\toutput #%d", index_tau_out);
  //       else
  //         fprintf (stderr, "\toutput #%d (z=%g)",
  //           index_tau_out, ppt2->z_out[index_tau_out-(ppt2->tau_out_size-ppt2->z_out_size)]);
  //   }
  //   fprintf (stderr, "\n");
  // }

  free (pvecback);
  free (pvecthermo);

  return _SUCCESS_;

}// end of galbispectra2_timesampling_for_sources()



/* Define integrand we wish to integrate over (integrand of angular power spectrum) */
double integrand(
  struct background * pba, /* input pointer to background */
  struct bessels * pbs, /* input, pointer to bessels structure */
  struct galbispectra2 * pgb2, /* input pointer to galaxy bispectra structure */
  struct perturbs * ppt, /* input pointer to galaxy bispectra structure */
  int index_type, /* index corresponding to type (e.g. delta_m) */
  int index_tau_first,
  int index_tau_second, /* index corresponding to ppt->tau_sampling[index_tau] */
  int index_k,
  int index_l) /**< order of the Bessel functions from pbs->l[index_l] */
  {
    double j1;
    double j2;
    double x1 = ppt->k[ppt->index_md_scalars][index_k]*(pba->conformal_age - ppt->tau_sampling[index_tau_first]);
    double x2 = ppt->k[ppt->index_md_scalars][index_k]*(pba->conformal_age - ppt->tau_sampling[index_tau_second]);
    class_call(bessel_at_x(pbs,	x1, index_l, &j1),pgb2->error_message, pbs->error_message);
    class_call(bessel_at_x(pbs, x2 , index_l, &j2), pgb2->error_message, pbs->error_message);

    double result = pow(ppt->k[ppt->index_md_scalars][index_k],-4.) * 4. * _PI_ * pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau_first][index_k]
      * pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau_second][index_k]
        * j1 * j2 ;

    return result;
  };

/* Define integral over k to give angular power spectrum without window functions */
double integral(
  struct background * pba, /* input pointer to background */
  struct bessels * pbs, /* input, pointer to bessels structure */
  struct galbispectra2 * pgb2, /* input pointer to galaxy bispectra structure */
  struct perturbs * ppt,
  int index_type, /* index corresponding to type (e.g. delta_m) */
  int index_tau_first,
  int index_tau_second, /* index corresponding to ppt->tau[index_tau] */
  int index_l) /**< order of the Bessel functions from pbs->l[index_l] */
  {
  int index_md;
  double tmp = 0.;
  int k_size = ppt->k_size[ppt->index_md_scalars];


  class_alloc(pgb2->w_trapz_k,
              k_size * sizeof(double),
              ppt->error_message);

  class_call(array_trapezoidal_weights(ppt->k[ppt->index_md_scalars],
                                       k_size,
                                       pgb2->w_trapz_k,
                                       pgb2->error_message),
                                       pgb2->error_message,
                                       pgb2->error_message);

                                       /*
  for (int index_k; index_k < k_size; index_k++){
    printf("pgb2->w_trapz_k[index_k] = %g\n", pgb2->w_trapz_k[index_k]);
  }
  */
  for(int index_k = 0; index_k < k_size; index_k++){

/* Using area of trapezoid, we can sum a number of areas of trapezoids to approximate the integral */

    tmp+=integrand(pba, pbs, pgb2, ppt, index_type, index_tau_first, index_tau_second, index_k, index_l)
      * pgb2->w_trapz_k[index_k];

    }
  return tmp;
} //End of integral()

/*
  double angbispectrum(
    struct background * pba, /* input pointer to background */
    //struct bessels * pbs, /* input, pointer to bessels structure */
    //struct galbispectra2 * pgb2, /* input pointer to galaxy bispectra structure */
    //struct perturbations * ppt, /* input pointer to galaxy bispectra structure */
    //int index_type, /* index corresponding to type (e.g. delta_m) */
    //int index_tau1, /* index corresponding to integral(struct background * pba, /* input pointer to background */
    //struct bessels * pbs, /* input, pointer to bessels structure */
    //struct galbispectra2 * pgb2, /* input pointer to galaxy bispectra structure */
    //struct perturbations * ppt, /* input pointer to galaxy bispectra structure */
    //int index_type, /* index corresponding to type (e.g. delta_m) */
    //int index_tau_first,
    //int index_tau_second, /* index corresponding to ppt->tau_sampling[index_tau] */
  ///  int index_k, /* index corresponding to ppt->k[index_md][index_k] */
    //int index_l) /**< order of the Bessel functions from pbs->l[index_l] */
  //  {

    //int index_tau2, /* index corresponding to ppt->tau_sampling[index_tau] */
    //int index_tau3,
    //int index_k;
    //int index_k1, /* index corresponding to ppt->k[index_md][index_k] */
    //int index_k2, /* index corresponding to ppt->k[index_md][index_k] */
  ///  int index_k3, /* index corresponding to ppt->k[index_md][index_k] */
    //int index_l1, /**< order of the Bessel functions from pbs->l[index_l] */
  //  int index_l2, /**< order of the Bessel functions from pbs->l[index_l] */
//    int index_l3)


  /* Different combinations of spherical bessel arguments */

    /* NOTE: It is probably more efficient to loop over index_type, index_k and index_tau */
    /* NOTE: If transfer functions are distinct, there needs to be 12 pairs of C_l's */

  /*  double tmp = 2.; /* integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau2, index_k1, index_l2) * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau3, index_k2, index_l3)
      * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau2, index_k1, index_l1) * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau3, index_k2, index_l3)
        * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau2, index_k1, index_l2) * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau3, index_k2, index_l1)
          * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau2, index_k1, index_l3) * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau3, index_k2, index_l2)
            * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau2, index_k1, index_l3) * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau3, index_k2, index_l1)
              * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau2, index_k1, index_l1) * integral(pba, pbs, pgb2, ppt, index_type, index_tau1, index_tau3, index_k2, index_l2);
*/

  //  return index_k;

//  };






/* Initialise the module, this function is called by the song.c main file. */

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
      double tau0 = pba->conformal_age;
      double * w_trapz;
      double * tau0_minus_tau;
  printf("We are here 1!\n");
  /*
  if (ppt2->has_perturbations2 == _FALSE_) {

    if (ppt2->perturbations2_verbose > 0)
      printf("No second-order sources requested. Second-order perturbations module skipped.\n");

    return _SUCCESS_;
  }
  */
  printf("We are here 2!\n");
  class_alloc(pvecback,pba->bg_size*sizeof(double), pgb2->error_message);
  if (ppt2->perturbations2_verbose > 0)
    printf("Computing second-order perturbations\n");




  // ====================================================================================
  // =                              Indices and samplings                               =
  // ====================================================================================

  /* Determine which sources need to be computed and their k-sampling */
  printf("We are here 3!\n");

/*  class_call( perturb2_indices_of_perturbs(
                ppr,
                ppr2,
                pba,
                pth,
                ppt,
                ppt2),
  ppt2->error_message,
  ppt2->error_message); */

  printf("We are here 4!\n");

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

  printf("We are here 5!\n");

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

  printf("We are here 6!\n");

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

  printf("We are here 7!\n");
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

  printf("We are here 8!\n");

  /* Do not evaluate the subsequent modules if ppt2->stop_at_perturbations2 == _TRUE_ */
  if (ppt2->stop_at_perturbations2 == _TRUE_) {
    ppt->has_perturbations = _FALSE_;
    ppt->has_cls = _FALSE_;
    ppt->has_cmb_bispectra = _FALSE_;
    ppt2->has_cmb_spectra = _FALSE_;
    ppt2->has_cmb_bispectra = _FALSE_;
  }

  printf("We are here 11!\n");




  /* Define an array of values of first order transfer functions:
              pgb2->first_order_sources[index_type][index_tau][index_k] */

  class_alloc(pgb2->first_order_sources, ppt->qs_size[ppt->index_md_scalars] * sizeof(double **), ppt->error_message);
    for (int index_type = 0; index_type < ppt->qs_size[ppt->index_md_scalars]; index_type++) {
      /* Allocate memory for pgb2->first_order_sources[index_type][index_tau] */
      class_alloc(pgb2->first_order_sources[index_type],
                  ppt->tau_size * sizeof(double *),
                  ppt->error_message);
      /* Allocate memory for pgb2->first_order_sources[index_type] */
      for (int index_tau = 0; index_tau < ppt->tau_size; index_tau++) {
          /* Loop over types, tau and k. For each of them, allocate memory
           for pgb2->first_order_sources[index_type][index_tau][index_k]  */
           class_alloc(pgb2->first_order_sources[index_type][index_tau],
                      ppt->k_size[ppt->index_md_scalars] * sizeof(double),
                      ppt->error_message);
    }
  }

  printf("We are here 12!\n");
  printf("here %p %p \n",ptr,ptr->l_size);

  /* Allocate array for Cl[index_l][index_tau_first][index_tau_second] */
  class_alloc(pgb2->Cl, ptr->l_size[ppt->index_md_scalars] * sizeof(double **), ppt->error_message);
  for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
    class_alloc(pgb2->Cl[index_l], ppt->tau_size * sizeof(double *), ppt->error_message);
    for (int index_tau_first = 0; index_tau_first < ppt->tau_size; index_tau_first++) {
      class_alloc(pgb2->Cl[index_l][index_tau_first], ppt->tau_size * sizeof(double), ppt->error_message);
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

  printf("We are here 13!\n");
/*
  class_alloc(pgb2->index_tau_first,
              ppt->tau_size * sizeof(double),
              ppt->error_message);

  class_alloc(pgb2->index_tau_second,
              ppt->tau_size * sizeof(double),
              ppt->error_message);
*/
  class_alloc(tau0_minus_tau,
              ppt->tau_size * sizeof(double),
              ppt->error_message);

  class_alloc(w_trapz,
              ppt->tau_size * sizeof(double),
              ppt->error_message);

  class_call(array_trapezoidal_mweights(tau0_minus_tau,ppt->tau_size,w_trapz,pgb2->error_message),
                                        ppt2->error_message,
                                        ppt2->error_message);
  printf("We are here 14!\n");

/* Fill the array ptw->tau0_minus_tau[index_tau] */

  for(int index_tau = 0; index_tau < ppt->tau_size; index_tau++){
    tau0_minus_tau[index_tau] = pba->conformal_age-ppt->tau_sampling[index_tau];
  }
  /* Declaration of temporary pointer */
  double ** selection;

  /* Allocation of first dimension selection[bin] */
  class_alloc(selection,
              ppt->selection_num * sizeof(double*),
              ppt->error_message);

  printf("ppt->selection_num = %d\n",ppt->selection_num);

  /* Allocation of second dimension selection[bin][index_tau] */
  for(int bin = 0; bin < ppt->selection_num; bin++){
    printf("bin = %d\n",bin);
    class_alloc(selection[bin],
                ppt->tau_size * sizeof(double),
                ppt->error_message);
    /* transfer_selection_compute prints in to selection[bin] */
    class_call(transfer_selection_compute(ppr, pba, ppt, ptr, selection[bin], tau0_minus_tau, w_trapz, ppt->tau_size, pvecback, tau0, bin),
               pgb2->error_message,
               pgb2->error_message);
    printf("selection[%d] = %g\n", bin, selection[bin]);
  }

  printf("We are here 15!\n");
  printf("ppt->tp_size[index_md] = %d\n", ppt->tp_size[index_md]);

  printf("ppt->tp_size[ppt->index_md_scalars] = %d\n", ppt->tp_size[ppt->index_md_scalars]);
  printf("ppt->index_ic_ad = %d\n", ppt->index_ic_ad);
  printf("ppt->index_md_scalars = %d\n",ppt->index_md_scalars);
  printf("ppt->index_tp_delta_m = %d\n", ppt->index_qs_delta_cdm);
  printf("tau_size = %d\n", ppt->tau_size);
  printf("k_size = %d\n", ppt->k_size[index_md]);


  printf("We are here 16!\n");
  printf("%g\n",ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][100 * ppt->k_size[ppt->index_md_scalars] + 20]);
  //for (int index_type = 0; index_type < ppt->tp_size[index_md]; index_type++) {
    for (int index_tau = 0; index_tau < ppt->tau_size; index_tau++) {
      for (int index_k = 0; index_k < ppt->k_size[ppt->index_md_scalars]; index_k++) {
        pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau][index_k] =
          ppt->quadsources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->qs_size[ppt->index_md_scalars]+ppt->index_qs_delta_cdm][index_tau * ppt->k_size[ppt->index_md_scalars] + index_k];


      }
       /*printf("first_order_sources[%d][%d][%d] = %g , %g\n",ppt->index_qs_delta_cdm, index_tau, index_k,
        pgb2->first_order_sources[ppt->index_qs_delta_cdm][index_tau][10]);*/

    }
  //}


/*  for (index_md = 0; index_md < ppt->md_size; index_md++) {
    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
      for (index_type = 0; index_type < ppt->tp_size[index_md]; index_type++) {

        class_alloc(ppt->sources[index_md][index_ic*ppt->tp_size[index_md]+index_type],
                    ppt->k_size[index_md] * ppt->tau_size * sizeof(double),
                    ppt->error_message);

      }
    }
  }*/

  printf("We are here 17!\n");

  /* Now fill this new pointer-array pgb2->first_order_sources[index_type][index_tau][index_k]
    with information from the pre computed ppt->sources[index_md][index_ic*ppt->tp_size[index_md]+index_type]
    [index_tau * ppt->k_size[index_md] + index_k]. */


  /* Set the pointer pgb2->first_order_sources equal to ppt->sources array. */

  for (int index_type = 0; index_type < ppt->tp_size[index_md]; index_type++) {
    for (int index_tau = 0; index_tau < ppt->tau_size; index_tau++) {
      for (int index_k = 0; index_k < ppt->k_size[index_md]; index_k++) {
        pgb2->first_order_sources[index_type][index_tau][index_k] =
          ppt->sources[ppt->index_md_scalars][ppt->index_ic_ad*ppt->tp_size[ppt->index_md_scalars]+index_type][index_tau * ppt->k_size[ppt->index_md_scalars] + index_k];


      }
    }
  }
  printf("We are here 18!\n");

  ppt->selection_mean[0] = 1.0;
  ppt->selection_mean[1] = 1.5;
  ppt->selection_mean[2] = 2.0;



  /* Integrate over our integrand w.r.t. k*/

  printf("ptr->l_size[ppt->index_md_scalar] = %d\n", ptr->l_size[ppt->index_md_scalars] );
  printf("ppt->selection_num = %d\n", ppt->selection_num);
  printf("ppt->selection_mean[0] = %g\n", ppt->selection_mean[0]);
  printf("ppt->selection_width[0] = %g\n", ppt->selection_width[0]);



  /* Write the result of the angular power spectrum integral into an array */
  for(int index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
    for (int index_tau_first = 0; index_tau_first < ppt->tau_size; index_tau_first++){
      for(int index_tau_second = 0; index_tau_second < ppt->tau_size; index_tau_second++){
        pgb2->Cl[index_l][index_tau_first][index_tau_second] = integral(pba, pbs, pgb2, ppt, ppt->index_qs_delta_cdm, index_tau_first, index_tau_second, index_l);
        //printf("Cl[index_l = %d][index_tau_first = %d][index_tau_second = %d] = %g\n",index_l, index_tau_first, index_tau_second, pgb2->Cl[index_l][index_tau_first][index_tau_second]);
        //printf("%g\n", integrand(pba,pbs,pgb2, ppt, index_type, index_tau_first, index_tau_second, index_k, index_l));
        //printf("integral(index_l = %d, index_tau_first = %d, index_tau_second = %d) = %g\n",index_l,index_tau_first, index_tau_second, integral(pba, pbs, pgb2, ppt, ppt->index_qs_delta_cdm, index_tau_first, index_tau_second, index_l));

      }
    }
  }
    printf("We are here 19!\n");
// NOTE: Do these bins need to be linked to the .ini file that lists the bins? selection_num needs to be defined correctly (check it's not null).

  int  bin1 = 0;
  int  bin2 = 0;


   /*
  for(int index_l = 0; index_l < 120; index_l++){
  pgb2->Cl_final[index_l][ppt->selection_num][ppt->selection_num] += pgb2->Cl[index_l][50][50] * w_trapz[50] * w_trapz[50]
      * selection[ppt->selection_num][50] * selection[ppt->selection_num][50];
  }
  */

  /* Integrate over the window function */
  for(int index_l =0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
    pgb2->Cl_final[index_l][bin1][bin2] = 0;
    for (index_tau_first = 0; index_tau_first < ppt->tau_size; index_tau_first++){
      for(index_tau_second = 0; index_tau_second < ppt->tau_size; index_tau_second++){
        pgb2->Cl_final[index_l][bin1][bin2] += pgb2->Cl[index_l][index_tau_first][index_tau_second] * w_trapz[index_tau_first] * w_trapz[index_tau_second]
            * selection[ppt->selection_num][index_tau_first] * selection[ppt->selection_num][index_tau_second];
        printf("pgb2->Cl_final[index_l = %d][bin1 = %d][bin2 = %d] = %g\n", ppt->selection_num, ppt->selection_num, pgb2->Cl_final[index_l][bin1][bin2]);
      }
    }
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


  for(index_l = 0; index_l < ptr->l_size[ppt->index_md_scalars]; index_l++){
    for (index_tau_first = 0; index_tau_first < ppt->tau_size; index_tau_first++){
      for(index_tau_second = 0; index_tau_second < ppt->tau_size; index_tau_second++){
        printf("Test dens-dens angular power spectrum C_%d(%g)(%g) \n", ptr->l[index_l] , ppt->tau_sampling[index_tau_first], ppt->tau_sampling[index_tau_second],
          pgb2->Cl_final[index_l][index_tau_first][index_tau_second]);
      }
    }
  }

printf("End of galbispectra2!\n");
  return _SUCCESS_;

} // end of galbispectra2_init()
