/** @file thermodynamics.c Documented thermodynamics module
 *
 * Julien Lesgourgues, 6.09.2010    
 *
 * Deals with the thermodynamical evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the thermodynamics, i.e. to
 *   integrate the thermodynamical equations, and store all
 *   thermodynamical quantities as a function of redshift inside an
 *   interpolation table. The current version of recombination is
 *   based on RECFAST v1.5. The current version of reionization is
 *   based on exactly the same reionization function as in CAMB, in
 *   order to make allow for comparison. It should be easy to
 *   generalize the module to more complicated reionization histories.
 *
 * - to provide a routine which allow other modules to evaluate any
 *   thermodynamical quantitites at a given redshft value (by
 *   interpolating within the interpolation table).
 *
 * 
 * The logic is the following:
 *
 * - in a first step, the code assumes that there is no reionization,
 *   and computes the ionization fraction, Thomson scattering rate,
 *   baryon temperature, etc., using RECFAST. The result is stored in
 *   a temporary table 'recombination_table' (within a temporary
 *   structure of type 'recombination') for each redshfit in a range 0
 *   < z < z_initial.  The sampling in z space is done with a simple
 *   linear step size.
 * - in a second step, the code adds the reionization history,
 *   starting from a redshift z_reio_start. The ionization fraction at
 *   this redshift is read in the previous recombination table in
 *   order to ensure a perfect matching. The code computes the
 *   ionization fraction, Thomson scattering rate, baryon temperature,
 *   etc., using a given parametrization of the reionization
 *   history. The result is stored in a temporary table
 *   'reionization_table' (within a temporary structure of type
 *   'reionization') for each redshift in the range 0 < z <
 *   z_reio_start. The sampling in z space is found automatically,
 *   given the precision parameter 'reionization_sampling'.
 *
 * - in a third step, the code merges the two tables
 *   'recombination_table' and 'reionization_table' inside the table
 *   'thermodynamics_table', and the temporary structures
 *   'recombination' and 'reionization' are freed. In
 *   'thermodynamics_table', the sampling in z space is the one
 *   defined in the recombination algorithm for z_reio_start < z <
 *   z_initial, and the one defined in the reionization algorithm for
 *   0 < z < z_reio_start.
 * 
 * - at this stage, only a few columns in the table
 *   'thermodynamics_table' have been filled. In a fourth step, the
 *   remaining columns are filled, using some numerical
 *   integration/derivation routines from the 'array.c' tools module.
 * 
 * - small detail: one of the columns contains the maximum variation
 *   rate of a few relevant thermodynamical quantitites. This rate
 *   will be used for defining automatically the sampling step size in
 *   the perturbation module. Hence, the exact value of this rate is
 *   unimportant, but its order of magnitude at a given z defines the
 *   sampling precision of the perturbation module. Hence, it is
 *   harmless to use a smoothing routine in order to make this rate
 *   look nicer, although this will not affect the final result
 *   significantly. The last step in the thermodynamics_init module is
 *   to perform this smoothing. 
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# thermodynamics_init() at the beginning (but after background_init()) 
 * -# thermodynamics_at_z() at any later time
 * -# thermodynamics_free() at the end, when no more calls to thermodynamics_at_z() are needed
 */

#include "thermodynamics.h"

#ifdef HYREC
#include "hyrec.h"
#endif

/** 
 * Thermodynamics quantities at given redshift z. 
 *
 * Evaluates all thermodynamics quantities at a given value of
 * the redshift by reading the pre-computed table and interpolating.
 *
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure (containing pre-computed table)
 * @param z          Input: redshift
 * @param intermode  Input: interpolation mode (normal or growing_closeby)
 * @param last_index Input/Ouput: index of the previous/current point in the interpolation array (input only for closeby mode, output for both) 
 * @param pvecback   Input: vector of background quantitites (used only in case z>z_initial for getting ddkappa and dddkappa; in that case, should be already allocated and filled, with format short_info or larger; in other cases, will be ignored)
 * @param pvecthermo Output: vector of thermodynamics quantities (assumed to be already allocated)
 * @return the error status
 */

int thermodynamics_at_z(
      struct background * pba,
      struct thermo * pth,
      double z,
      short inter_mode,
      int * last_index,
      double * pvecback,
      double * pvecthermo
      ) {

  /** Summary: */

  /** - define local variables */

  double x0;

  /* - the fact that z is in the pre-computed range 0 <= z <= z_initial
     will be checked in the interpolation routines below. Before
     trying to interpolate, allow the routine to deal with the case z
     > z_intial: then, all relevant quantitites can be extrapolated
     using simple analytic approximations */

  if (z >= pth->z_table[pth->tt_size-1]) {

    /* ionization fraction assmued to reamin constant at large z. it is larger than one if we consider Helium. */
    x0 = pth->thermodynamics_table[(pth->tt_size-1)*pth->th_size+pth->index_th_xe];
    pvecthermo[pth->index_th_xe] = x0;

    /* Calculate dkappa/dtau (dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T in units of 1/Mpc) */
    pvecthermo[pth->index_th_dkappa] = (1.+z) * (1.+z) * pth->n_e * x0 * _sigma_ * _Mpc_over_m_;

    /* Calculate d2kappa/dtau2 = dz/dtau d/dz[dkappa/dtau] given that [dkappa/dtau] proportional to (1+z)^2 and dz/dtau = -H */
    pvecthermo[pth->index_th_ddkappa] = -pvecback[pba->index_bg_H] * 2. / (1.+z) * pvecthermo[pth->index_th_dkappa];

    /* Calculate d3kappa/dtau3 given that [dkappa/dtau] proportional to (1+z)^2 */
    pvecthermo[pth->index_th_dddkappa] = (pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]/ (1.+z) - pvecback[pba->index_bg_H_prime]) * 2. / (1.+z) * pvecthermo[pth->index_th_dkappa];

    /* \f$ exp^{-\kappa}, g, g', g'' \f$ can be set to zero: they are
       used only for computing the source functions in the
       perturbation module; but source functions only need to be
       sampled below z_initial (the condition that
       z_start_sources<z_initial is checked in the perturbation
       module) */
    pvecthermo[pth->index_th_exp_m_kappa] = 0.;
    pvecthermo[pth->index_th_g]=0.;
    pvecthermo[pth->index_th_dg]=0.;
    pvecthermo[pth->index_th_ddg]=0.;

    /* Calculate Tb */
    pvecthermo[pth->index_th_Tb] = pba->T_cmb*(1.+z);

    /* Calculate cb2 (cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz)) */
    /* note that m_H / mu = 1 + (m_H/m_He-1) Y_p + x_e (1-Y_p) */
    pvecthermo[pth->index_th_cb2] = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + x0 * (1.-pth->YHe)) * pba->T_cmb * (1.+z) * 4. / 3.;

    /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
    if (pth->compute_cb2_derivatives == _TRUE_) {

      /* since cb2 proportional to (1+z) or 1/a, its derivative wrt conformal time is given by dcb2 = - a H cb2 */ 
      pvecthermo[pth->index_th_dcb2] = - pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] * pvecthermo[pth->index_th_cb2];
      
      /* then its second derivative is given by ddcb2 = - a H' cb2 */ 
      pvecthermo[pth->index_th_ddcb2] = - pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a] * pvecthermo[pth->index_th_cb2];
    }
    
    /* in this regime, variation rate = dkappa/dtau */
    pvecthermo[pth->index_th_rate] = pvecthermo[pth->index_th_dkappa];
    

    // *** MY MODIFICATIONS ***
    
    // *** Perturbed recombination ***
    if (pth->compute_xe_derivatives == _TRUE_) {
      
      /* The ionization fraction is assumed to be constant at very large z (>10000) */
      pvecthermo[pth->index_th_dxe] = 0;
      pvecthermo[pth->index_th_ddxe] = 0;      
    }
    
    if (pth->has_perturbed_recombination == _TRUE_) {

      /* It makes no sense to compute perturbed recombination early on, as the Q function features a discontinuity
      just before recombination, ad x ~ 34 (x is the inverse temperature normalized to the ionization energy of
      the hydrogen atom). The system, however, is stable enough that we can provide Q only for x > 34. */
      pvecthermo[pth->index_th_dQ_dx] = 0;
      pvecthermo[pth->index_th_dQ_dX] = 0;
      pvecthermo[pth->index_th_dQ_dn] = 0;
      pvecthermo[pth->index_th_dQ_dH] = 0;
      pvecthermo[pth->index_th_Q] = 0;

      /* Uncomment to compute the Q function and its derivatives at early times. Expect a lot of oscillations
      and numerical noise */
      // class_call (thermodynamics_compute_Q_derivatives_at_z (
      //               // ppr,
      //               pba,
      //               pth,
      //               x0,
      //               z),
      //   pth->error_message,
      //   pth->error_message);
      // 
      // pvecthermo[pth->index_th_dQ_dx] = pth->dQ_dx;
      // pvecthermo[pth->index_th_dQ_dX] = pth->dQ_dX;
      // pvecthermo[pth->index_th_dQ_dn] = pth->dQ_dn;
      // pvecthermo[pth->index_th_dQ_dH] = pth->dQ_dH;
      // pvecthermo[pth->index_th_Q] = pth->Q;

    }
    
    // *** Rayleigh scattering ***
    if (pth->has_rayleigh_scattering == _TRUE_) {

      /* The neutral hydrogen and helium fractions are negligible because everything is ionised in the early Universe,
      hence there is no Rayleigh scattering. */
      pvecthermo[pth->index_th_rayleigh_dkappa] = 0;
      pvecthermo[pth->index_th_rayleigh_ddkappa] = 0;
      pvecthermo[pth->index_th_rayleigh_dddkappa] = 0;
      pvecthermo[pth->index_th_rayleigh_exp_m_kappa] = 0; /* TODO: can we safely assume that kappa(tau_0)>>1? */
      pvecthermo[pth->index_th_rayleigh_g] = 0;
      pvecthermo[pth->index_th_rayleigh_dg] = 0;
      pvecthermo[pth->index_th_rayleigh_ddg] = 0;
      
      /* The Thomson scattering rate is equal to the total scattering rate */
      pvecthermo[pth->index_th_thomson_dkappa] = pvecthermo[pth->index_th_dkappa];
      pvecthermo[pth->index_th_thomson_ddkappa] = pvecthermo[pth->index_th_ddkappa];
      pvecthermo[pth->index_th_thomson_dddkappa] = pvecthermo[pth->index_th_dddkappa];
      pvecthermo[pth->index_th_thomson_exp_m_kappa] = pvecthermo[pth->index_th_exp_m_kappa];
      pvecthermo[pth->index_th_thomson_g] = pvecthermo[pth->index_th_g];
      pvecthermo[pth->index_th_thomson_dg] = pvecthermo[pth->index_th_dg];
      pvecthermo[pth->index_th_thomson_ddg] = pvecthermo[pth->index_th_ddg];

    }

    
    // *** END OF MY MODIFICATIONS ***    
    
  }

  /** - interpolate in table with array_interpolate_spline() (normal
      mode) or array_interpolate_spline_growing_closeby() (closeby
      mode) */

  else {

    if (inter_mode == pth->inter_normal) {
      
      class_call(array_interpolate_spline(
            pth->z_table,
            pth->tt_size,
            pth->thermodynamics_table,
            pth->d2thermodynamics_dz2_table,
            pth->th_size,
            z,
            last_index,
            pvecthermo,
            pth->th_size,
            pth->error_message),
     pth->error_message,
     pth->error_message);
      
    }
    
    if (inter_mode == pth->inter_closeby) {
      
      class_call(array_interpolate_spline_growing_closeby(
                pth->z_table,
                pth->tt_size,
                pth->thermodynamics_table,
                pth->d2thermodynamics_dz2_table,
                pth->th_size,
                z,
                last_index,
                pvecthermo,
                pth->th_size,
                pth->error_message),
     pth->error_message,
     pth->error_message);
    }
  }
  return _SUCCESS_;
}

/** 
 * Initialize the thermo structure, and in particular the
 * thermodynamics interpolation table.
 * 
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background structure
 * @param pth Input/Output : pointer to initialized thermo structure
 * @return the error status
 */
int thermodynamics_init(
      struct precision * ppr,
      struct background * pba,
      struct thermo * pth
      ) {

  /** Summary: */

  /** - define local variables */

  /* index running over time*/
  int index_tau;
  /* temporary variables related to visibility function */
  double g;
  /* vector of background values for calling background_at_tau() */
  double * pvecback;
  /* index for calling background_at_tau() */
  int last_index_back;
  /* temporary table of values of tau associated with z values in pth->z_table */
  double * tau_table;
  /* conformal time of reionization */
  double tau_reio;
  /* structures for storing temporarily information on recombination
     and reionization */
  struct recombination reco;
  struct reionization reio;
  struct recombination * preco;
  struct reionization * preio;

  double tau;

  /** - initialize pointers, allocate background vector */

  preco=&reco;
  preio=&reio;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  if (pth->thermodynamics_verbose > 0)
    printf("Computing thermodynamics");

  /** - compute and check primordial Helium fraction  */

  /* Y_He */
  if (pth->YHe == _BBN_) {
    class_call(thermodynamics_helium_from_bbn(ppr,pba,pth),
         pth->error_message,
         pth->error_message);
    if (pth->thermodynamics_verbose > 0)
      printf(" with Y_He=%.4f\n",pth->YHe);
  }
  else {
    if (pth->thermodynamics_verbose > 0)
      printf("\n");
  }

  class_test((pth->YHe < _YHE_SMALL_)||(pth->YHe > _YHE_BIG_),
        pth->error_message,
        "Y_He=%g out of bounds (%g<Y_He<%g)",pth->YHe,_YHE_SMALL_,_YHE_BIG_);

  /* tests in order to prevent segmentation fault in the following */
  class_test(_not4_ == 0.,
       pth->error_message,
       "stop to avoid division by zero");
  class_test(pth->YHe == 1.,
       pth->error_message,
       "stop to avoid division by zero");

  /** - assign values to all indices in the structures with thermodynamics_indices()*/

  class_call(thermodynamics_indices(pth,preco,preio),
       pth->error_message,
       pth->error_message);

  /** - solve recombination and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f $ with thermodynamics_recombination() */

  class_call(thermodynamics_recombination(ppr,pba,pth,preco,pvecback),
       pth->error_message,
       pth->error_message);

  /** - if there is reionization, solve reionization and store values of \f$ z, x_e, d \kappa / d \tau, T_b, c_b^2 \f $ with thermodynamics_reionization()*/

  if (pth->reio_parametrization != reio_none) {
    class_call(thermodynamics_reionization(ppr,pba,pth,preco,preio,pvecback),
         pth->error_message,
         pth->error_message);
  }
  else {
    preio->rt_size=0;
    preio->index_reco_when_reio_start=-1;
  }

  /** - merge tables in recombination and reionization structures into
        a single table in thermo structure */

  class_call(thermodynamics_merge_reco_and_reio(ppr,pth,preco,preio),
       pth->error_message,
       pth->error_message);

  /** - compute table of corresponding conformal times */

  class_alloc(tau_table,pth->tt_size*sizeof(double),pth->error_message);

  for (index_tau=0; index_tau < pth->tt_size; index_tau++) {
    class_call(background_tau_of_z(pba,pth->z_table[index_tau],tau_table+index_tau),
         pba->error_message,
         pth->error_message);
  }

  /** - store initial value of conformal time in the structure */

  pth->tau_ini = tau_table[pth->tt_size-1];

  /** - fill missing columns (quantities not computed previously but related) */

  /** -> second derivative with respect to tau of dkappa (in view of spline interpolation) */
  class_call(array_spline_table_line_to_line(tau_table,
               pth->tt_size,
               pth->thermodynamics_table,
               pth->th_size,
               pth->index_th_dkappa,
               pth->index_th_dddkappa,
               _SPLINE_EST_DERIV_,
               pth->error_message),
       pth->error_message,
       pth->error_message);

  /** -> first derivative with respect to tau of dkappa (using spline interpolation) */
  class_call(array_derive_spline_table_line_to_line(tau_table,
                pth->tt_size,
                pth->thermodynamics_table,
                pth->th_size,
                pth->index_th_dkappa,
                pth->index_th_dddkappa,
                pth->index_th_ddkappa,
                pth->error_message),
       pth->error_message,
       pth->error_message);

  /** -> compute -kappa = [int_{tau_today}^{tau} dtau dkappa/dtau], store temporarily in column "g" */ 
  class_call(array_integrate_spline_table_line_to_line(tau_table,
                   pth->tt_size,
                   pth->thermodynamics_table,
                   pth->th_size,
                   pth->index_th_dkappa,
                   pth->index_th_dddkappa,
                   pth->index_th_g,
                   pth->error_message),
       pth->error_message,
       pth->error_message);
  
  /** -> derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives == _TRUE_) {

    /** -> second derivative with respect to tau of cb2 */
    class_call(array_spline_table_line_to_line(tau_table,
                 pth->tt_size,
                 pth->thermodynamics_table,
                 pth->th_size,
                 pth->index_th_cb2,
                 pth->index_th_ddcb2,
                 _SPLINE_EST_DERIV_,
                 pth->error_message),
         pth->error_message,
         pth->error_message);
    
    
    /** -> first derivative with respect to tau of cb2 (using spline interpolation) */
    class_call(array_derive_spline_table_line_to_line(tau_table,
                  pth->tt_size,
                  pth->thermodynamics_table,
                  pth->th_size,
                  pth->index_th_cb2,
                  pth->index_th_ddcb2,
                  pth->index_th_dcb2,
                  pth->error_message),
         pth->error_message,
         pth->error_message);
  }


  /** -> derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_xe_derivatives == _TRUE_) {

    /** -> compute_xe_derivatives respect to tau of cb2 */
    class_call(array_spline_table_line_to_line(tau_table,
                 pth->tt_size,
                 pth->thermodynamics_table,
                 pth->th_size,
                 pth->index_th_xe,
                 pth->index_th_ddxe,
                 _SPLINE_EST_DERIV_,
                 pth->error_message),
         pth->error_message,
         pth->error_message);
    
    
    /** -> first derivative with respect to tau of cb2 (using spline interpolation) */
    class_call(array_derive_spline_table_line_to_line(tau_table,
                  pth->tt_size,
                  pth->thermodynamics_table,
                  pth->th_size,
                  pth->index_th_xe,
                  pth->index_th_ddxe,
                  pth->index_th_dxe,
                  pth->error_message),
         pth->error_message,
         pth->error_message);
  }

  /** -> compute visibility : \f$ g= (d \kappa/d \tau) e^{- \kappa} */

  /* loop on z (decreasing z, increasing time) */
  for (index_tau=pth->tt_size-1; index_tau>=0; index_tau--) {

    /** -> compute g */  
    g = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);
    
    /** -> compute exp(-kappa) */    
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_exp_m_kappa] = 
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);

    /** -> compute g' (the plus sign of the second term is correct, see def of -kappa in thermodynamics module!) */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dg] = 
      (pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] +
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa]) *
      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_exp_m_kappa];
    
    /** -> compute g''  */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddg] = 
      (pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dddkappa] +
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa] * 3. +
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] *
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa]) *
      exp(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]);

    /** -> store g */
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g] = g;

    /** -> compute variation rate */
    class_test(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] == 0.,
         pth->error_message,
         "variation rate diverges");

    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rate] =
      sqrt(pow(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa],2)
     +pow(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_ddkappa]/ 
    pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa],2)
     +fabs(pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dddkappa]/ 
     pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa]));

  }


  // *** MY MODIFICATIONS ***

  /* Compute g and e^-k for the Rayleigh scattering. This is basically the same computation done above
  for the Thomson scattering with minimal modifications. */
  if (pth->has_rayleigh_scattering == _TRUE_) {
  
    /** -> second derivative with respect to tau of dkappa (in view of spline interpolation). Make sure
    to use the _SPLINE_NATURAL_ condition because for Rayleigh scattering kappa_dot=0 at early times. */
    class_call(array_spline_table_line_to_line(
                 tau_table,
                 pth->tt_size,
                 pth->thermodynamics_table,
                 pth->th_size,
                 pth->index_th_rayleigh_dkappa,
                 pth->index_th_rayleigh_dddkappa,
                 _SPLINE_NATURAL_,
                 pth->error_message),
         pth->error_message,
         pth->error_message);
      
    /** -> first derivative with respect to tau of dkappa (using spline interpolation) */
    class_call(array_derive_spline_table_line_to_line(
                  tau_table,
                  pth->tt_size,
                  pth->thermodynamics_table,
                  pth->th_size,
                  pth->index_th_rayleigh_dkappa,
                  pth->index_th_rayleigh_dddkappa,
                  pth->index_th_rayleigh_ddkappa,
                  pth->error_message),
         pth->error_message,
         pth->error_message);
      
    /** -> compute -kappa = [int_{tau_today}^{tau} dtau dkappa/dtau], store temporarily in column "g" */ 
    class_call(array_integrate_spline_table_line_to_line(
                     tau_table,
                     pth->tt_size,
                     pth->thermodynamics_table,
                     pth->th_size,
                     pth->index_th_rayleigh_dkappa,
                     pth->index_th_rayleigh_dddkappa,
                     pth->index_th_rayleigh_g,
                     pth->error_message),
         pth->error_message,
         pth->error_message);
         

    /* loop on z (decreasing z, increasing time) */
    for (index_tau=pth->tt_size-1; index_tau>=0; index_tau--) {
      
      /* Compute e^-k factors */
      double kappa_rayleigh = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_g];
      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_exp_m_kappa] = exp(kappa_rayleigh);
      double exp_m_kappa_total = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_exp_m_kappa];
      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_exp_m_kappa] = 
        exp_m_kappa_total / pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_exp_m_kappa];
      
      /* Compute visibility functions */  
      double dkappa_rayleigh = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_dkappa];
      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_g] = 
        pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_exp_m_kappa] * dkappa_rayleigh;
      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_dkappa] = 
        pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_dkappa] - dkappa_rayleigh;
      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_g] = 
        pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_exp_m_kappa] * 
        pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_dkappa];

      /* Compute the derivative of the Rayleigh visibility function */
      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_dg] = 
        (pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_ddkappa] +
         pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_dkappa] *
         pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_dkappa]) *
        pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_rayleigh_exp_m_kappa];

    } // end of loop on z
  
  
    /* Compute the derivative of the Thomson visibility function */
    class_call(array_spline_table_line_to_line(
                 tau_table,
                 pth->tt_size,
                 pth->thermodynamics_table,
                 pth->th_size,
                 pth->index_th_thomson_dkappa,
                 pth->index_th_thomson_dddkappa,
                 _SPLINE_EST_DERIV_,
                 pth->error_message),
         pth->error_message,
         pth->error_message);
      
    class_call(array_derive_spline_table_line_to_line(
                  tau_table,
                  pth->tt_size,
                  pth->thermodynamics_table,
                  pth->th_size,
                  pth->index_th_thomson_dkappa,
                  pth->index_th_thomson_dddkappa,
                  pth->index_th_thomson_ddkappa,
                  pth->error_message),
         pth->error_message,
         pth->error_message);
      
    for (index_tau=pth->tt_size-1; index_tau>=0; index_tau--) {
      pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_dg] = 
        (pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_ddkappa] +
         pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_dkappa] *
         pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_dkappa]) *
        pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_exp_m_kappa];
    
    }
    
  } // end of Rayleigh scattering computation

  // *** END OF MY MODIFICATIONS ***



  /** -> smooth the rate (dtauils of smoothing unimportant: only the
         order of magnitude of the rate matters) */ 
  class_call(array_smooth(pth->thermodynamics_table,
        pth->th_size,
        pth->tt_size,
        pth->index_th_rate,
        ppr->thermo_rate_smoothing_radius,
        pth->error_message),
       pth->error_message,
       pth->error_message);

  /** - fill tables of second derivatives with respect to z (in view of spline interpolation) */

  class_call(array_spline_table_lines(pth->z_table,
              pth->tt_size,
              pth->thermodynamics_table,
              pth->th_size,
              pth->d2thermodynamics_dz2_table,
              _SPLINE_EST_DERIV_,
              pth->error_message),
       pth->error_message,
       pth->error_message);

  /** - find maximum of g */

  index_tau=pth->tt_size-1;
  while (pth->z_table[index_tau]>_Z_REC_MAX_) {
    index_tau--;
  }

  class_test(pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g] >
       pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g],
       pth->error_message,
       "found a recombination redshift greater or equal to the maximum value imposed in thermodynamics.h, z_rec_max=%g",_Z_REC_MAX_);

  while (pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g] <
   pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_g]) {
    index_tau--;
  }
  
  /* approximation for maximum of g, using cubic interpolation, assuming equally spaced z's */
  pth->z_rec=pth->z_table[index_tau+1]+0.5*(pth->z_table[index_tau+1]-pth->z_table[index_tau])*(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g]-pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_g])/(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_g]-2.*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_g]+pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_g]);

  class_test(pth->z_rec+ppr->smallest_allowed_variation >= _Z_REC_MAX_,
       pth->error_message,
       "found a recombination redshift greater or equal to the maximum value imposed in thermodynamics.h, z_rec_max=%g",_Z_REC_MAX_);

  class_test(pth->z_rec-ppr->smallest_allowed_variation <= _Z_REC_MIN_,
       pth->error_message,
       "found a recombination redshift smaller or equal to the maximum value imposed in thermodynamics.h, z_rec_min=%g",_Z_REC_MIN_);


   // *** MY MODIFICATIONS ***

   /** - find what would have been the maximum of g without Rayleigh scattering */
   if (pth->has_rayleigh_scattering == _TRUE_) {

     int index_tau;

     index_tau=pth->tt_size-1;
     while (pth->z_table[index_tau]>_Z_REC_MAX_) {
       index_tau--;
     }

     class_test(pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_thomson_g] >
   	     pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_g],
   	     pth->error_message,
   	     "found a recombination redshift greater or equal to the maximum value imposed in thermodynamics.h, z_rec_max=%g",_Z_REC_MAX_);

     while (pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_thomson_g] <
   	 pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_g]) {
       index_tau--;
     }

     double g_max = pth->thermodynamics_table[index_tau*pth->th_size+pth->index_th_thomson_g];
     int index_tau_max = index_tau;

     /* approximation for maximum of g, using cubic interpolation, assuming equally spaced z's */
     pth->z_rec_thomson=pth->z_table[index_tau+1]+0.5*(pth->z_table[index_tau+1]-pth->z_table[index_tau])*(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_thomson_g]-pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_thomson_g])/(pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_thomson_g]-2.*pth->thermodynamics_table[(index_tau+1)*pth->th_size+pth->index_th_thomson_g]+pth->thermodynamics_table[(index_tau+2)*pth->th_size+pth->index_th_thomson_g]);

   }
   // *** END OF MY MODIFICATIONS ***

  /** - find conformal recombination time using background_tau_of_z() **/

  class_call(background_tau_of_z(pba,pth->z_rec,&(pth->tau_rec)),
       pba->error_message,
       pth->error_message);

  class_call(background_at_tau(pba,pth->tau_rec, pba->long_info, pba->inter_normal, &last_index_back, pvecback),
       pba->error_message,
       pth->error_message);

  pth->rs_rec=pvecback[pba->index_bg_rs];
  pth->ds_rec=pth->rs_rec*pba->a_today/(1.+pth->z_rec);
  pth->da_rec=pvecback[pba->index_bg_ang_distance];

  /** - find time (always after recombination) at which tau_c/tau
        falls below some threshold, defining tau_free_streaming */

  class_call(background_tau_of_z(pba,pth->z_table[index_tau],&tau),
       pba->error_message,
       pth->error_message);
  
  while (1./pth->thermodynamics_table[(index_tau)*pth->th_size+pth->index_th_dkappa]/tau 
   < ppr->radiation_streaming_trigger_tau_c_over_tau) {

    index_tau--;
    
    class_call(background_tau_of_z(pba,pth->z_table[index_tau],&tau),
         pba->error_message,
         pth->error_message);
    
  }

  pth->tau_free_streaming = tau;


  // *** MY MODIFICATIONS ***
  
  /* Compute functions needed to obtain the first-order fraction of free electrons */
  if (pth->has_perturbed_recombination == _TRUE_) {

    class_call (thermodynamics_compute_Q_derivatives (
                  ppr,
                  pba,
                  pth),
      pth->error_message,
      pth->error_message);
    
  }

  /* Debug - print thermodynamical quantities to file or screen */
  // if (pth->has_rayleigh_scattering == _FALSE_) {
  // 
  //   fprintf (stderr, "%16s %16s %16s %16s %16s\n",
  //     "tau","z","dkappa_t","g_t","exp_m_kappa_t\n");
  //   int i;
  //   for (i=0; i < pth->tt_size; ++i) {
  // 
  //     if (i%4!=0) continue;
  // 
  //     double z = pth->z_table[i];
  //     double tau;
  // 
  //     class_call(background_tau_of_z(pba,z,&tau),
  //            pba->error_message,
  //            pth->error_message);
  // 
  //     double dkappa_t = pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa];
  //     double g_t = pth->thermodynamics_table[i*pth->th_size+pth->index_th_g];
  //     double exp_m_kappa_t = pth->thermodynamics_table[i*pth->th_size+pth->index_th_exp_m_kappa];
  // 
  //     fprintf (stderr, "%16g %16g %16g %16g %16g\n",
  //       tau, z, dkappa_t, g_t, exp_m_kappa_t);
  //   }
  // 
  // }
  // else {
  // 
  //   fprintf (stderr, "%16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s \n",
  //     "tau","z","dkappa_t","dkappa_r","dkappa_tot","g_t","g_r","g_tot","dg_t","dg_r","dg_tot","exp_m_kappa_t","exp_m_kappa_r","exp_m_kappa_tot","ddkappa_t","ddkappa_r","ddkappa_tot","dddkappa_t","dddkappa_r","dddkappa_tot");
  //   int i;
  //   for (i=0; i < pth->tt_size; ++i) {
  // 
  //     if (i%4!=0) continue;
  // 
  //     double z = pth->z_table[i];
  //     double tau;
  // 
  //     class_call(background_tau_of_z(pba,z,&tau),
  //            pba->error_message,
  //            pth->error_message);
  // 
  //     double dkappa_t = pth->thermodynamics_table[i*pth->th_size+pth->index_th_thomson_dkappa];
  //     double dkappa_r = pth->thermodynamics_table[i*pth->th_size+pth->index_th_rayleigh_dkappa];
  //     double dkappa_tot = pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa];
  // 
  //     double ddkappa_t = pth->thermodynamics_table[i*pth->th_size+pth->index_th_thomson_ddkappa];
  //     double ddkappa_r = pth->thermodynamics_table[i*pth->th_size+pth->index_th_rayleigh_ddkappa];
  //     double ddkappa_tot = pth->thermodynamics_table[i*pth->th_size+pth->index_th_ddkappa];
  // 
  //     double dddkappa_t = pth->thermodynamics_table[i*pth->th_size+pth->index_th_thomson_dddkappa];
  //     double dddkappa_r = pth->thermodynamics_table[i*pth->th_size+pth->index_th_rayleigh_dddkappa];
  //     double dddkappa_tot = pth->thermodynamics_table[i*pth->th_size+pth->index_th_dddkappa];
  // 
  //     double g_t = pth->thermodynamics_table[i*pth->th_size+pth->index_th_thomson_g];
  //     double g_r = pth->thermodynamics_table[i*pth->th_size+pth->index_th_rayleigh_g];
  //     double g_tot = pth->thermodynamics_table[i*pth->th_size+pth->index_th_g];
  // 
  //     double dg_t = pth->thermodynamics_table[i*pth->th_size+pth->index_th_thomson_dg];
  //     double dg_r = pth->thermodynamics_table[i*pth->th_size+pth->index_th_rayleigh_dg];
  //     double dg_tot = pth->thermodynamics_table[i*pth->th_size+pth->index_th_dg];
  // 
  //     double exp_m_kappa_t = pth->thermodynamics_table[i*pth->th_size+pth->index_th_thomson_exp_m_kappa];
  //     double exp_m_kappa_r = pth->thermodynamics_table[i*pth->th_size+pth->index_th_rayleigh_exp_m_kappa];
  //     double exp_m_kappa_tot = pth->thermodynamics_table[i*pth->th_size+pth->index_th_exp_m_kappa];
  //   
  //     fprintf (stderr, "%16g %16g %16g %16g %16g %16g %16g %16g %16g %16g %16g %16g %16g %16g %16g %16g %16g %16g %16g %16g \n",
  //       tau, z, dkappa_t, dkappa_r, dkappa_tot, g_t, g_r, g_tot, dg_t, dg_r, dg_tot, exp_m_kappa_t, exp_m_kappa_r, exp_m_kappa_tot, ddkappa_t, ddkappa_r, ddkappa_tot, dddkappa_t, dddkappa_r, dddkappa_tot);
  //   }
  // }

  // *** END OF MY MODIFICATIONS ***


  /** - if verbose flag set to next-to-minimum value, print the main results */

  if (pth->thermodynamics_verbose > 0) {
    printf(" -> recombination at z = %f\n",pth->z_rec);
    printf("    corresponding to conformal time = %f Mpc\n",pth->tau_rec);
    printf("    with sound horizon = %f Mpc\n",pth->ds_rec);
    printf("    with comoving sound horizon = %f Mpc\n",pth->rs_rec);
    printf("    and angular diameter distance = %f Mpc\n",pth->da_rec);
    // *** MY MODIFICATIONS ***
    if (pth->has_rayleigh_scattering == _TRUE_) {
      printf(" -> Rayleigh scattering with nu = %f Ghz\n",pth->rayleigh_frequency);
      printf("    moves recombination from z=%f to z=%f\n",pth->z_rec,pth->z_rec_thomson);
    }
    // *** END OF MY MODIFICATIONS ***
    if (pth->reio_parametrization == reio_camb) {
      if (pth->reio_z_or_tau==reio_tau) 
  printf(" -> reionization  at z = %f\n",pth->z_reio);
      if (pth->reio_z_or_tau==reio_z)
  printf(" -> reionization with optical depth = %f\n",pth->tau_reio);
      class_call(background_tau_of_z(pba,pth->z_reio,&tau_reio),
     pba->error_message,
     pth->error_message);
      printf("    corresponding to conformal time = %f Mpc\n",tau_reio);
    }
    if (pth->reio_parametrization == reio_bins_tanh) {
      printf(" -> binned reionization gives optical depth = %f\n",pth->tau_reio);
    }
    if (pth->thermodynamics_verbose > 1) {
      printf(" -> free-streaming approximation can be turned on as soon as tau=%g Mpc\n",
       pth->tau_free_streaming);
    }
  }

  free(pvecback);
  free(tau_table);

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by thermodynamics_init().
 * 
 *
 * @param pth Input/Output : pointer to thermo structure (to be freed)
 * @return the error status
 */

int thermodynamics_free(
      struct thermo * pth
      ) {

  free(pth->z_table);
  free(pth->thermodynamics_table);
  free(pth->d2thermodynamics_dz2_table);

  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of thermodynamical quantities,
 * as well as in vector containing reionization parameters.
 *
 *
 * @param pth   Input/Output: pointer to thermo structure
 * @param preco Input/Output: pointer to recombination structure
 * @param preio Input/Output: pointer to reionization structure 
 * @return the error status
 */

int thermodynamics_indices(
         struct thermo * pth,
         struct recombination * preco,
         struct reionization * preio
         ) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of thermodynamics quantities */
  int index;

  /** - intialization of all indices and flags in thermo structure */
  index = 0;
  
  pth->index_th_xe = index;
  index++;
  
  // *** MY MODIFICATIONS ***

  /* Derivatives of ionization fraction (needed for the approximated perturbed recombination) */  
  if (pth->compute_xe_derivatives == _TRUE_) {

    pth->index_th_dxe=index++;
    pth->index_th_ddxe = index++;
  }

  /* Derivatives of the Q function, needed for the full perturbed recombination */
  if (pth->has_perturbed_recombination == _TRUE_) {
  
    pth->index_th_dQ_dx = index++;
    pth->index_th_dQ_dX = index++;
    pth->index_th_dQ_dn = index++;
    pth->index_th_dQ_dH = index++;
    pth->index_th_Q = index++;
  }
  
  // *** END OF MY MODIFICATIONS ***
  
  
  pth->index_th_dkappa = index;
  index++;
  pth->index_th_ddkappa = index;
  index++;
  pth->index_th_dddkappa = index;
  index++;
  pth->index_th_exp_m_kappa = index;
  index++;
  pth->index_th_g = index;
  index++;
  pth->index_th_dg = index;
  index++;
  pth->index_th_ddg = index;
  index++;
  pth->index_th_Tb = index;
  index++;
  pth->index_th_cb2 = index;
  index++;
  // *** MY MODIFICATIONS ***
  if (pth->has_rayleigh_scattering == _TRUE_) {
    pth->index_th_rayleigh_dkappa = index++;
    pth->index_th_rayleigh_ddkappa = index++;
    pth->index_th_rayleigh_dddkappa = index++;
    pth->index_th_rayleigh_exp_m_kappa = index++;
    pth->index_th_rayleigh_g = index++;
    pth->index_th_rayleigh_dg = index++;
    pth->index_th_rayleigh_ddg = index++;

    pth->index_th_thomson_dkappa = index++;
    pth->index_th_thomson_ddkappa = index++;
    pth->index_th_thomson_dddkappa = index++;
    pth->index_th_thomson_exp_m_kappa = index++;
    pth->index_th_thomson_g = index++;
    pth->index_th_thomson_dg = index++;
    pth->index_th_thomson_ddg = index++;
  }
  // *** END OF MY MODIFICATIONS ***

  /* derivatives of baryon sound speed (only computed if some non-minimal tight-coupling schemes is requested) */
  if (pth->compute_cb2_derivatives == _TRUE_) {
    pth->index_th_dcb2 = index;
    index++;
    pth->index_th_ddcb2 = index;
    index++;
  }

  pth->index_th_rate = index;
  index++;
  

  /* end of indices */
  pth->th_size = index;

  /** - intialization of all indices and flags in recombination structure */
  index = 0;

  preco->index_re_z = index;
  index++;
  preco->index_re_xe = index;
  index++;
  preco->index_re_dkappadtau = index;
  index++;
  preco->index_re_Tb = index;
  index++;
  preco->index_re_cb2 = index;
  index++;
  // *** MY MODIFICATIONS ***
  if (pth->has_rayleigh_scattering == _TRUE_)
    preco->index_re_rayleigh_dkappadtau = index++;
  // *** END OF MY MODIFICATIONS ***

  /* end of indices */
  preco->re_size = index;

  /** - intialization of all indices and flags in reionization structure */
  index = 0;

  preio->index_re_z = index;
  index++;
  preio->index_re_xe = index;
  index++;
  preio->index_re_Tb = index;
  index++;
  preio->index_re_cb2 = index;
  index++;
  preio->index_re_dkappadtau = index;
  index++;
  preio->index_re_dkappadz = index;
  index++;
  preio->index_re_d3kappadz3 = index;
  index++;
  // *** MY MODIFICATIONS ***
  if (pth->has_rayleigh_scattering == _TRUE_)
    preio->index_re_rayleigh_dkappadtau = index++;
  // *** END OF MY MODIFICATIONS ***


  /* end of indices */
  preio->re_size = index;

  /* same with parameters of the function x_e(z) */
  
  index=0;

  preio->index_reio_start = index;
  index++;

  /* case where x_e(z) taken like in CAMB (other cases can be added) */ 
  if (pth->reio_parametrization == reio_camb) {

    preio->index_reio_redshift = index;
    index++;
    preio->index_reio_exponent = index;
    index++;
    preio->index_reio_width = index;
    index++;
    preio->index_reio_xe_before = index;
    index++;
    preio->index_reio_xe_after = index;
    index++;
    preio->index_helium_fullreio_fraction = index;
    index++;
    preio->index_helium_fullreio_redshift = index;
    index++;
    preio->index_helium_fullreio_width = index;
    index++;

  }

  /* case where x_e(z) is binned */ 
  if (pth->reio_parametrization == reio_bins_tanh) {

    /* the code will not only copy here the "bin centers" passed in
       input. It will add an initial and final value for (z,xe). So
       this array has a dimension bigger than the bin center array */

    preio->reio_num_z=pth->binned_reio_num+2; /** add two values: beginning and end of reio */

    preio->index_reio_first_z = index;
    index+= preio->reio_num_z;
    preio->index_reio_first_xe = index;
    index+= preio->reio_num_z; 
    preio->index_reio_step_sharpness = index;
    index++; 

    preio->reio_num_params = index;
  }

  preio->reio_num_params = index;

  /* flags for calling the interpolation routine */

  pth->inter_normal=0;
  pth->inter_closeby=1;

  return _SUCCESS_;
}











// *** MY MODIFICATIONS ***
  
/**
 * Compute functions needed to obtain the perturbed fraction of free electrons, delta_Xe.
 * These quantities will be used in the perturbation module to construct the evolution
 * equation for delta_Xe.
 *
 * All computations are made in eV to better relate with what is done in CMBquick.
 */
int thermodynamics_compute_Q_derivatives_at_z (
          struct background * pba,
          struct thermo * pth,
          double X,
          double z
          )
{
  
  double recfast_fudge_H = 1.14;
      
  /* Obtain all time variables */
  double a = 1/(z+1);
  double T_kelvin = pba->T_cmb/a;
  double T_ev = T_kelvin * _kBoltz_;
  double x = _EI_/T_ev;

  /* Some debug - print all time variables */
  // double tau;
  // class_call(background_tau_of_z(pba, z, &tau),
  //   pba->error_message,
  //   pth->error_message);
  // fprintf(stderr, "%3d %12.7g %12.7g %12.7g %12.7g %12.7g %12.7g\n", index_tau, tau, z, a, T_kelvin, T_ev, x);

  /* Allocate vector of background quantities */
  double * pvecback;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
    
  /* Compute background quantities at the given redshift */
  class_call(background_functions(pba, a, pba->long_info, pvecback),
    pba->error_message,
    pth->error_message);

  /* Factor that converts the Hubble factor in pvecback (defined as [H/c] in inverse Mpc) to
    the Hubble factor in electron volts (eV). */ 
  double inverse_seconds_to_eV = (_c_ / _Mpc_over_m_) / _s_in_inverse_eV_;
    
  /* Hubble factor in eV (the one in pvecback is [H/c] in inverse Mpc) */
  double H = pvecback[pba->index_bg_H] * inverse_seconds_to_eV;
  double H0 = pba->H0 * inverse_seconds_to_eV;
  double Hc = a*H;

  /* Baryon density in eV^-3 */
  double H0_MKS = pba->H0 * (_c_ / _Mpc_over_m_);
  double rho_crit_MKS = 3*H0_MKS*H0_MKS/(8.*_PI_*_G_);
  double n_b = pba->Omega0_b/pow(a/pba->a_today, 3) * (rho_crit_MKS/_m_H_/pow(_m_in_eV_,3));

  /* Helium fraction */
  double Y = pth->YHe;

  /* Some debug */
  // fprintf (stderr, "%g %g %g %g %g %g\n", z, a, x, X, n_b, H);


  // *** Derived thermodynamics quantities ***

  /* G_alp, depending only on H */
  double G_alp_constant = pow(3*_EI_,3)/pow(8*_PI_,2);
  double G_alp = H * G_alp_constant;

  /* alpha_2p, depending only on x */
  double alpha_2p_constant = 0.448*64*sqrt(_PI_/27)*(_FINE_/_m_e_eV_)*(_FINE_/_m_e_eV_);
  double alpha_2p = alpha_2p_constant * sqrt(x) * log(x);
    
  /* beta_p, depending only on x */
  double beta_o_constant = pow(_m_e_eV_*_EI_/(2*_PI_), 1.5);
  double beta_p = beta_o_constant * alpha_2p * exp(-0.25*x) / pow(x, 1.5);
    
  /* As done in CMBquick, we divide the Q function in a front coefficient and a source term,
    that is we write Q = F * S. */
        
  /* Transition rate of 2S->1S in electronvolts */
  double lambda = _L2s1s_/_s_in_inverse_eV_;
      
  /* We further divide the front coefficient as F = F1*F2 in order to better deal with the
    derivatives. Note that F2 features a singularity around x~30. */
  double F1 = 1 + (1-X)*(1-Y)*lambda/G_alp * n_b;
  double lambda_beta_factor = lambda/recfast_fudge_H + beta_p;
  double F2 = 1./(1/recfast_fudge_H + (1-X)*(1-Y)*n_b/G_alp * lambda_beta_factor);
  double F = F1*F2;
    
  /* To compute the source term, we need the function beta_o. */
  double beta_o = beta_o_constant * alpha_2p * exp(-x) / pow(x, 1.5);
  double S = 1/X * ((1-X)*beta_o - X*X*(1-Y)*n_b*alpha_2p);
    
  /* Finally, we compute and store the background Q function */
  double Q = F*S;
  pth->Q = Q/inverse_seconds_to_eV;
    
  /* Some debug */
  // fprintf (stderr, "%g %g %g %g %g %g\n", x, F, S, Q, beta_o, (1-X)*(1-Y)*n_b/G_alp * lambda_beta_factor);
    


  // *** Derivatives of thermodynamics quantities ***
    
  /* Derivative with respect to X. We add a pre-factor to account for the fact that X^(1) = X * delta_X,
  where X = X_e and delta_X = delta_X_e.  */
  double dF_dX  = - (1-Y)/G_alp * n_b * F2 * (lambda - F1*F2 * lambda_beta_factor);
  double dS_dX  = - 1/X * (S + beta_o + 2*X*(1-Y)*n_b*alpha_2p);
  double dQ_dX  = F*dS_dX + S*dF_dX;
  pth->dQ_dX = X * dQ_dX/inverse_seconds_to_eV;

  /* Derivative with respect to H. We add a pre-factor to account for the fact that H^(1) = H * delta_H,
  where delta_H is given by eq. A. 67 of Pitrou et al. 2010.  */
  double dF_dH  = -(1-X)*(1-Y)/(G_alp*G_alp) * G_alp_constant * n_b * F2 * (lambda - F1*F2 * lambda_beta_factor);
  double dS_dH  = 0;
  double dQ_dH  = F*dS_dH + S*dF_dH;
  pth->dQ_dH = H * dQ_dH/inverse_seconds_to_eV;

  /* Derivative with respect to n_b. We add a pre-factor to account for the fact that n_b^(1) = n_b * delta_b */
  double dF_dn  = (1-X)*(1-Y)/G_alp * F2 * (lambda - F1*F2 * lambda_beta_factor);
  double dS_dn  = - X * (1-Y) * alpha_2p;
  double dQ_dn  = F*dS_dn + S*dF_dn;
  pth->dQ_dn = n_b * dQ_dn/inverse_seconds_to_eV;

  /* Derivative with respect to x. We add a pre-factor to account for the fact that x^(1) = -x/4 * delta_g */
  double d_alpha2p_dx = alpha_2p_constant * (sqrt(x)/x + log(x)/(2*sqrt(x)));
  double d_betao_dx = - beta_o * (1 + 1.5/x - d_alpha2p_dx/alpha_2p);
  double d_betap_dx = - beta_p * (0.25 + 1.5/x - d_alpha2p_dx/alpha_2p);
  double dF_dx  = -F1 * F2 * F2 * (1-X) * (1-Y) * n_b/G_alp * d_betap_dx;
  double dS_dx  = 1/X * ((1-X)*d_betao_dx - X*X*(1-Y)*n_b*d_alpha2p_dx);
  double dQ_dx  = F*dS_dx + S*dF_dx;
  pth->dQ_dx = -x/4 * dQ_dx/inverse_seconds_to_eV;

  /* Some debug */
  // fprintf (stderr, "%g %g %g %g %g %g \n", x, dQ_dX, dQ_dH, dQ_dn, dQ_dx, dF_dx);

  free (pvecback);
  
  return _SUCCESS_;
  
}


int thermodynamics_compute_Q_derivatives(
          struct precision * ppr,
          struct background * pba,
          struct thermo * pth
          )
{
  
  /* Temporary table with the conformal time values corresponding to the redshift steps in pth->z_table */
  double * tau_table;
  class_alloc (tau_table, pth->tt_size*sizeof(double), pth->error_message);

  /* Assign the times to each redshift, sorting the times in ascending order */
  int index_tau;
  for (index_tau=0; index_tau < pth->tt_size; ++index_tau)
    class_call(background_tau_of_z(pba, pth->z_table[index_tau], tau_table+index_tau),
      pba->error_message,
      pth->error_message);
  
  /* Allocate vector of background quantities */
  double * pvecback;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  
  /* Loop on time (decreasing z, increasing time) */
  for (index_tau=pth->tt_size-1; index_tau>=0; --index_tau) {
    
    /* Obtain all time variables */
    double tau = tau_table[index_tau];
    double z = pth->z_table[index_tau];
    double a = 1/(z+1);
    double T_kelvin = pba->T_cmb/a;
    double T_ev = T_kelvin * _kBoltz_;
    double x = _EI_/T_ev;

    /* Some debug - print all time variables */
    // fprintf(stderr, "%3d %12.7g %12.7g %12.7g %12.7g %12.7g %12.7g\n", index_tau, tau, z, a, T_kelvin, T_ev, x);
    
    /* Compute background quantities at the given time */
    class_call(background_functions(pba, a, pba->long_info, pvecback),
      pba->error_message,
      pth->error_message);


    /* The Q function features a discontinuity around x ~ 34, hence we integrate only after that time. */
    if (x < pth->perturbed_recombination_turnx) {
      pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_Q] = 0;
      pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_dQ_dX] = 0;
      pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_dQ_dH] = 0;
      pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_dQ_dn] = 0;
      pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_dQ_dx] = 0;
      continue;
    }

    // *** Simple thermodynamics quantities ***

    /* Fraction of free electrons */
    double X = pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_xe];

    /* Factor that converts the Hubble factor in pvecback (defined as [H/c] in inverse Mpc) to
      the Hubble factor in electron volts (eV). */ 
    double inverse_seconds_to_eV = (_c_ / _Mpc_over_m_) / _s_in_inverse_eV_;
    
    /* Hubble factor in eV (the one in pvecback is [H/c] in inverse Mpc) */
    double H = pvecback[pba->index_bg_H] * inverse_seconds_to_eV;
    double H0 = pba->H0 * inverse_seconds_to_eV;
    double Hc = a*H;

    /* Baryon density in eV^-3 */
    double H0_MKS = pba->H0 * (_c_ / _Mpc_over_m_);
    double rho_crit_MKS = 3*H0_MKS*H0_MKS/(8.*_PI_*_G_);
    double n_b = pba->Omega0_b/pow(a/pba->a_today, 3) * (rho_crit_MKS/_m_H_/pow(_m_in_eV_,3));

    /* Helium fraction */
    double Y = pth->YHe;

    /* Some debug */
    // fprintf (stderr, "%g %g %g %g %g %g\n", z, a, x, X, n_b, H);


    // *** Derived thermodynamics quantities ***

    /* G_alp, depending only on H */
    double G_alp_constant = pow(3*_EI_,3)/pow(8*_PI_,2);
    double G_alp = H * G_alp_constant;

    /* alpha_2p, depending only on x */
    double alpha_2p_constant = 0.448*64*sqrt(_PI_/27)*(_FINE_/_m_e_eV_)*(_FINE_/_m_e_eV_);
    double alpha_2p = alpha_2p_constant * sqrt(x) * log(x);
    
    /* beta_p, depending only on x */
    double beta_o_constant = pow(_m_e_eV_*_EI_/(2*_PI_), 1.5);
    double beta_p = beta_o_constant * alpha_2p * exp(-0.25*x) / pow(x, 1.5);
    
    /* As done in CMBquick, we divide the Q function in a front coefficient and a source term,
      that is we write Q = F * S. */
        
    /* Transition rate of 2S->1S in electronvolts */
    double lambda = _L2s1s_/_s_in_inverse_eV_;
      
    /* We further divide the front coefficient as F = F1*F2 in order to better deal with the
      derivatives. Note that F2 features a singularity around x~30. */
    double F1 = 1 + (1-X)*(1-Y)*lambda/G_alp * n_b;
    double lambda_beta_factor = lambda/ppr->recfast_fudge_H + beta_p;
    double F2 = 1./(1/ppr->recfast_fudge_H + (1-X)*(1-Y)*n_b/G_alp * lambda_beta_factor);
    double F = F1*F2;
    
    /* To compute the source term, we need the function beta_o. */
    double beta_o = beta_o_constant * alpha_2p * exp(-x) / pow(x, 1.5);
    double S = 1/X * ((1-X)*beta_o - X*X*(1-Y)*n_b*alpha_2p);
    
    /* Finally, we compute and store the background Q function */
    double Q = F*S;
    pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_Q] = Q/inverse_seconds_to_eV;
    
    /* Some debug */
    // fprintf (stderr, "%g %g %g %g %g %g\n", x, F, S, Q, beta_o, (1-X)*(1-Y)*n_b/G_alp * lambda_beta_factor);
    


    // *** Derivatives of thermodynamics quantities ***
    
    /* Derivative with respect to X. We add a pre-factor to account for the fact that X^(1) = X * delta_Xe,
      where X = Xe and delta_X = delta_Xe  */
    double dF_dX  = - (1-Y)/G_alp * n_b * F2 * (lambda - F1*F2 * lambda_beta_factor);
    double dS_dX  = - 1/X * (S + beta_o + 2*X*(1-Y)*n_b*alpha_2p);
    double dQ_dX  = F*dS_dX + S*dF_dX;
    pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_dQ_dX] = X * dQ_dX/inverse_seconds_to_eV;

    /* Derivative with respect to H. We add a pre-factor to account for the fact that H^(1) = H * delta_H,
      where delta_H is given by eq. A. 67 of Pitrou et al. 2010.  */
    double dF_dH  = -(1-X)*(1-Y)/(G_alp*G_alp) * G_alp_constant * n_b * F2 * (lambda - F1*F2 * lambda_beta_factor);
    double dS_dH  = 0;
    double dQ_dH  = F*dS_dH + S*dF_dH;
    pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_dQ_dH] = H * dQ_dH/inverse_seconds_to_eV;

    /* Derivative with respect to n_b. We add a pre-factor to account for the fact that n_b^(1) = n_b * delta_b */
    double dF_dn  = (1-X)*(1-Y)/G_alp * F2 * (lambda - F1*F2 * lambda_beta_factor);
    double dS_dn  = - X * (1-Y) * alpha_2p;
    double dQ_dn  = F*dS_dn + S*dF_dn;
    pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_dQ_dn] = n_b * dQ_dn/inverse_seconds_to_eV;

    /* Derivative with respect to x. We add a pre-factor to account for the fact that x^(1) = -x/4 * delta_g */
    double d_alpha2p_dx = alpha_2p_constant * (sqrt(x)/x + log(x)/(2*sqrt(x)));
    double d_betao_dx = - beta_o * (1 + 1.5/x - d_alpha2p_dx/alpha_2p);
    double d_betap_dx = - beta_p * (0.25 + 1.5/x - d_alpha2p_dx/alpha_2p);
    double dF_dx  = -F1 * F2 * F2 * (1-X) * (1-Y) * n_b/G_alp * d_betap_dx;
    double dS_dx  = 1/X * ((1-X)*d_betao_dx - X*X*(1-Y)*n_b*d_alpha2p_dx);
    double dQ_dx  = F*dS_dx + S*dF_dx;
    pth->thermodynamics_table[index_tau*pth->th_size + pth->index_th_dQ_dx] = -x/4 * dQ_dx/inverse_seconds_to_eV;

    /* Some debug */
    // fprintf (stderr, "%g %g %g %g %g %g \n", x, dQ_dX, dQ_dH, dQ_dn, dQ_dx, dF_dx);

    
  } // end of for(index_tau)

  /* Free vectors */
  free (pvecback);
  free (tau_table);
  
  return _SUCCESS_;
  
}

// *** MY MODIFICATIONS ***
















/** 
 * Infer the primordial helium fraction from standard BBN, as a
 * function of the baryon density and expansion rate during BBN.
 *
 * This module is simpler then the one used in arXiv:0712.2826 because
 * it neglects the impact of a possible significant chemical
 * potentials for electron neutrinos. The full code with xi_nu_e could
 * be introduced here later.
 * 
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background structure
 * @param pth Input/Output : pointer to initialized thermo structure
 * @return the error status
 */
int thermodynamics_helium_from_bbn(
             struct precision * ppr,
             struct background * pba,
             struct thermo * pth
           ) {

  FILE * fA;
  char line[_LINE_LENGTH_MAX_];
  char * left;

  int num_omegab=0;
  int num_deltaN=0;

  double * omegab;
  double * deltaN;
  double * YHe;
  double * ddYHe;
  double * YHe_at_deltaN;
  double * ddYHe_at_deltaN;

  int array_line,i;
  double DeltaNeff;
  int last_index;

  /* compute Delta N_eff as defined in bbn file, i.e. Delta N_eff=0 means N_eff=3.046 */
  DeltaNeff = pba->Neff - 3.046;

  /* the following file is assumed to contain (apart from comments and blank lines):
     - the two numbers (num_omegab, num_deltaN) = number of values of BBN free paramters
     - three columns (omegab, deltaN, YHe) where omegab = Omega0_b h^2 and deltaN = Neff-3.046 by definition
     - omegab and deltaN are assumed to be arranged as:
          omegab1 deltaN1 YHe
    omegab2 deltaN1 YHe
                 .....
    omegab1 delatN2 YHe
    omegab2 deltaN2 YHe
                 .....
  */

  class_open(fA,ppr->sBBN_file, "r",pth->error_message);

  /* go through each line */
  while (fgets(line,_LINE_LENGTH_MAX_-1,fA) != NULL) {
    
    /* eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
      left++;
    }

    /* check that the line is neiyher blank neither a comment. In
       ASCII, left[0]>39 means that first non-blank charachter might
       be the beginning of some data (it is not a newline, a #, a %,
       etc.) */  
    if (left[0] > 39) {

      /* if the line contains data, we must interprete it. If
   (num_omegab, num_deltaN)=(0,0), the current line must contain
   their values. Otherwise, it must contain (omegab, delatN,
   YHe). */
      if ((num_omegab==0) && (num_deltaN==0)) {
  
  /* read (num_omegab, num_deltaN), infer size of arrays and allocate them */
  class_test(sscanf(line,"%d %d",&num_omegab,&num_deltaN) != 2,
       pth->error_message,
       "could not read value of parameters (num_omegab,num_deltaN) in file %s\n",ppr->sBBN_file);

  class_alloc(omegab,num_omegab*sizeof(double),pth->error_message);
  class_alloc(deltaN,num_deltaN*sizeof(double),pth->error_message);
  class_alloc(YHe,num_omegab*num_deltaN*sizeof(double),pth->error_message);
  class_alloc(ddYHe,num_omegab*num_deltaN*sizeof(double),pth->error_message);
  class_alloc(YHe_at_deltaN,num_omegab*sizeof(double),pth->error_message);
  class_alloc(ddYHe_at_deltaN,num_omegab*sizeof(double),pth->error_message);
        array_line=0;

      }
      else {

  /* read (omegab, deltaN, YHe) */
  class_test(sscanf(line,"%lg %lg %lg",
        &(omegab[array_line%num_omegab]),
        &(deltaN[array_line/num_omegab]),
        &(YHe[array_line])
        ) != 3,
       pth->error_message,
       "could not read value of parameters (omegab,deltaN,YHe) in file %s\n",ppr->sBBN_file);
  array_line ++;
      }
    }
  }

  fclose(fA);

  /* spline in one dimension (along deltaN) */
  class_call(array_spline_table_lines(deltaN,
              num_deltaN,
              YHe,
              num_omegab,    
              ddYHe,
              _SPLINE_NATURAL_,
              pth->error_message),
       pth->error_message,
       pth->error_message);

  /* interpolate in one dimension (along deltaN) */
  class_call(array_interpolate_spline(deltaN,
              num_deltaN,
              YHe,
              ddYHe,
              num_omegab,
              DeltaNeff,
              &last_index,
              YHe_at_deltaN,
              num_omegab,
              pth->error_message),
       pth->error_message,
       pth->error_message);

  /* spline in remaining dimension (along omegab) */
  class_call(array_spline_table_lines(omegab,
              num_omegab,
              YHe_at_deltaN,
              1,    
              ddYHe_at_deltaN,
              _SPLINE_NATURAL_,
              pth->error_message),
       pth->error_message,
       pth->error_message);

  /* interpolate in remaining dimension (along omegab) */
  class_call(array_interpolate_spline(omegab,
              num_omegab,
              YHe_at_deltaN,
              ddYHe_at_deltaN,
              1,
              pba->Omega0_b*pba->h*pba->h,
              &last_index,
              &(pth->YHe),
              1,
              pth->error_message),
       pth->error_message,
       pth->error_message);

  /* deallocate arrays */
  free(omegab);
  free(deltaN);
  free(YHe);
  free(ddYHe);
  free(YHe_at_deltaN);
  free(ddYHe_at_deltaN);

  return _SUCCESS_;

}


/**
 * This subroutine contains the reionization function \f$ X_e(z) \f$
 * (one for each scheme; so far, only the function corresponding to
 * the reio_camb scheme is coded)
 *
 * @param z     Input : redshift
 * @param pth   Input : pointer to thermo structure, to know which scheme is used
 * @param preio Input : pointer to reionization structure, containing the parameters of the function \f$ X_e(z) \f$
 * @param xe    Output: \f$ X_e(z) \f$
 */

int thermodynamics_reionization_function(
           double z,
           struct thermo * pth,
           struct reionization * preio,
           double * xe
           ) {

  /** Summary: */

  /** - define local variables */
  double argument;
  int i;

  /** - implementation of ionization function similar to the one in CAMB */

  if (pth->reio_parametrization == reio_camb) {

    /** -> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_start]) {
      *xe = preio->reionization_parameters[preio->index_reio_xe_before];
    }

    else {
      
      /** -> case z < z_reio_start: hydrogen contribution (tanh of complicated argument) */

      argument = (pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
          preio->reionization_parameters[preio->index_reio_exponent]) 
      - pow((1.+z),preio->reionization_parameters[preio->index_reio_exponent]))
  /(preio->reionization_parameters[preio->index_reio_exponent] 
    /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
    *pow((1.+preio->reionization_parameters[preio->index_reio_redshift]),
         (preio->reionization_parameters[preio->index_reio_exponent]-1.)))
  /preio->reionization_parameters[preio->index_reio_width]; 
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
      
      *xe = (preio->reionization_parameters[preio->index_reio_xe_after]
       -preio->reionization_parameters[preio->index_reio_xe_before])
  *(tanh(argument)+1.)/2.
  +preio->reionization_parameters[preio->index_reio_xe_before];
      
      /** -> case z < z_reio_start: helium contribution (tanh of simpler argument) */

      argument = (preio->reionization_parameters[preio->index_helium_fullreio_redshift] - z)
  /preio->reionization_parameters[preio->index_helium_fullreio_width];
      /* no possible segmentation fault: checked to be non-zero in thermodynamics_reionization() */
      *xe += preio->reionization_parameters[preio->index_helium_fullreio_fraction] 
  * (tanh(argument)+1.)/2.;

    }

    return _SUCCESS_;

  }

 /** - implementation of binned ionization function similar to astro-ph/0606552 */

  if (pth->reio_parametrization == reio_bins_tanh) {

    /** -> case z > z_reio_start */

    if (z > preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1];
    }

    else if (z < preio->reionization_parameters[preio->index_reio_first_z]) {
      *xe = preio->reionization_parameters[preio->index_reio_first_xe];
    }

    else {

      i = 0;
      while (preio->reionization_parameters[preio->index_reio_first_z+i+1]<z) i++;

      *xe = preio->reionization_parameters[preio->index_reio_first_xe+i]
  +0.5*(tanh((2.*(z-preio->reionization_parameters[preio->index_reio_first_z+i])
  /(preio->reionization_parameters[preio->index_reio_first_z+i+1]
    -preio->reionization_parameters[preio->index_reio_first_z+i])-1.)
       /preio->reionization_parameters[preio->index_reio_step_sharpness])
        /tanh(1./preio->reionization_parameters[preio->index_reio_step_sharpness])+1.)
  *(preio->reionization_parameters[preio->index_reio_first_xe+i+1]
    -preio->reionization_parameters[preio->index_reio_first_xe+i]);

      
    }

    return _SUCCESS_;

  }

  class_test(0 == 0,
       pth->error_message,
       "value of reio_parametrization=%d unclear",pth->reio_parametrization);
}

/**
 * This subroutine reads \f$ X_e(z) \f$ in the recombination table at
 * the time at which reionization starts. Hence it provides correct
 * initial conditions for the reionization function.
 *
 * @param ppr   Input : pointer to precision structure
 * @param pth   Input : pointer to thermo structure
 * @param preco Input : pointer to recombination structure
 * @param z     Input : redshift z_reio_start
 * @param xe    Output: \f$ X_e(z) \f$ at z
 */

int thermodynamics_get_xe_before_reionization(
                struct precision * ppr,
                struct thermo * pth,
                struct recombination * preco,
                double z,
                double * xe
                ) {

  int i;

  i=0;
  while (preco->recombination_table[i*preco->re_size+preco->index_re_z] < z) {
    i++;
    class_test(i == ppr->recfast_Nz0,
         pth->error_message,
         "z = %e > largest redshift in thermodynamics table \n",ppr->reionization_z_start_max);
  }
  
  *xe = preco->recombination_table[i*preco->re_size+preco->index_re_xe];

  return _SUCCESS_;

}
  

/** 
 * This routine computes the reionization history. In the reio_camb
 * scheme, this is straightforward if the input parameter is the
 * reionization redshift. If the input is the optical depth, need to
 * find z_reio by dichotomy (trying several z_reio until the correct
 * tau_reio is approached).
 *
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background structure
 * @param pth Input : pointer to thermo structure
 * @param preco Input : pointer to filled recombination structure
 * @param preio Input/Output: pointer to reionization structure (to be filled)
 * @param pvecback   Input: vector of background quantitites (used as workspace: must be already allocated, with format short_info or larger, but does not need to be filled)
 * @return the error status
 */

int thermodynamics_reionization(
        struct precision * ppr,
        struct background * pba,
        struct thermo * pth,
        struct recombination * preco,
        struct reionization * preio,
        double * pvecback
        ) {

  /** Summary: */

  /** - define local variables */

  int counter;
  double z_sup,z_mid,z_inf;
  double tau_sup,tau_mid,tau_inf;
  int bin;

  /** - allocate the vector of parameters defining the function \f$ X_e(z) \f$ */

  class_alloc(preio->reionization_parameters,preio->reio_num_params*sizeof(double),pth->error_message); 

  /** (a) if reionization implemented like in CAMB */

  if (pth->reio_parametrization == reio_camb) {
    
    /** - set values of these parameters, excepted those depending on the reionization redshift */

    preio->reionization_parameters[preio->index_reio_xe_after] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_reio_exponent] = pth->reionization_exponent; /* reio_exponent */
    preio->reionization_parameters[preio->index_reio_width] = pth->reionization_width;    /* reio_width */
    preio->reionization_parameters[preio->index_helium_fullreio_fraction] = pth->YHe/(_not4_*(1.-pth->YHe)); /* helium_fullreio_fraction (note: segmentation fault impossible, checked before that denominator is non-zero) */
    preio->reionization_parameters[preio->index_helium_fullreio_redshift] = pth->helium_fullreio_redshift; /* helium_fullreio_redshift */
    preio->reionization_parameters[preio->index_helium_fullreio_width] = pth->helium_fullreio_width;    /* helium_fullreio_width */

    class_test(preio->reionization_parameters[preio->index_reio_exponent]==0,
         pth->error_message,
         "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_reio_width]==0,
         pth->error_message,
         "stop to avoid division by zero");

    class_test(preio->reionization_parameters[preio->index_helium_fullreio_width]==0,
         pth->error_message,
         "stop to avoid division by zero");

    /** - if reionization redshift given as an input, initialize the remaining values and fill reionization table*/

    if (pth->reio_z_or_tau == reio_z) {
      
      /* reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = pth->z_reio; 
      /* infer starting redshift for hydrogen */
      preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_redshift]+ppr->reionization_start_factor*pth->reionization_width;
      /* if starting redshift for helium is larger, take that one
   (does not happen in realistic models) */
      if (preio->reionization_parameters[preio->index_reio_start] < 
    pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width)

  preio->reionization_parameters[preio->index_reio_start] =
    pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width;
  
      class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
     pth->error_message,
     "starting redshift for reionization > reionization_z_start_max = %e\n",ppr->reionization_z_start_max);

      /* infer xe_before_reio */
      class_call(thermodynamics_get_xe_before_reionization(ppr,
                 pth,
                 preco,
                 preio->reionization_parameters[preio->index_reio_start],
                 &(preio->reionization_parameters[preio->index_reio_xe_before])),
     pth->error_message,
     pth->error_message);

      /* fill reionization table */
      class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
     pth->error_message,
     pth->error_message);

      pth->tau_reio=preio->reionization_optical_depth;

    }

    /** - if reionization optical depth given as an input, find reionization redshift by dichotomy and initialize the remaining values */

    if (pth->reio_z_or_tau == reio_tau) {

      /* upper value */

      z_sup = ppr->reionization_z_start_max-ppr->reionization_start_factor*pth->reionization_width;
      class_test(z_sup < 0.,
     pth->error_message,
     "parameters are such that reionization cannot take place before today while starting after z_start_max; need to increase z_start_max");

      /* maximum possible reionization redshift */
      preio->reionization_parameters[preio->index_reio_redshift] = z_sup; 
      /* maximum possible starting redshift */
      preio->reionization_parameters[preio->index_reio_start] = ppr->reionization_z_start_max;
      /* infer xe_before_reio */
      class_call(thermodynamics_get_xe_before_reionization(ppr,
                 pth,
                 preco,
                 preio->reionization_parameters[preio->index_reio_start],
                 &(preio->reionization_parameters[preio->index_reio_xe_before])),
     pth->error_message,
     pth->error_message);

      /* fill reionization table */
      class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
     pth->error_message,
     pth->error_message);

      tau_sup=preio->reionization_optical_depth;

      class_test(tau_sup < pth->tau_reio,
     pth->error_message,
     "parameters are such that reionization cannot start after z_start_max");

      /* lower value */

      z_inf = 0.;
      tau_inf = 0.;      

      /* try intermediate values */
    
      counter=0;
      while ((tau_sup-tau_inf) > pth->tau_reio * ppr->reionization_optical_depth_tol) {
  z_mid=0.5*(z_sup+z_inf);
      
  /* reionization redshift */
  preio->reionization_parameters[preio->index_reio_redshift] = z_mid;
  /* infer starting redshift for hygrogen */
  preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_redshift]+ppr->reionization_start_factor*pth->reionization_width;
      /* if starting redshift for helium is larger, take that one
   (does not happen in realistic models) */
  if (preio->reionization_parameters[preio->index_reio_start] < 
      pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width)
    
    preio->reionization_parameters[preio->index_reio_start] =
      pth->helium_fullreio_redshift+ppr->reionization_start_factor*pth->helium_fullreio_width;
  
  class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
       pth->error_message,
       "starting redshift for reionization > reionization_z_start_max = %e",ppr->reionization_z_start_max);
  
  /* infer xe_before_reio */
  class_call(thermodynamics_get_xe_before_reionization(ppr,
                   pth,
                   preco,
                   preio->reionization_parameters[preio->index_reio_start],
                   &(preio->reionization_parameters[preio->index_reio_xe_before])),
       pth->error_message,
       pth->error_message);

  /* clean and fill reionization table */
  free(preio->reionization_table);
  class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
       pth->error_message,
       pth->error_message);

  tau_mid=preio->reionization_optical_depth;
  
  /* trial */

  if (tau_mid > pth->tau_reio) {
    z_sup=z_mid;
    tau_sup=tau_mid;
  }
  else {
    z_inf=z_mid;
    tau_inf=tau_mid;
  }

  counter++;
  class_test(counter > _MAX_IT_,
       pth->error_message,
       "while searching for reionization_optical_depth, maximum number of iterations exceeded");
      }

      /* store z_reionization in thermodynamics structure */
      pth->z_reio=preio->reionization_parameters[preio->index_reio_redshift];

    }

    free(preio->reionization_parameters);

    return _SUCCESS_;

  }

  if (pth->reio_parametrization == reio_bins_tanh) {

    /* this algorithm requires at least two bin centers (i.e. at least
       4 values in the (z,xe) array, counting the edges). */
    class_test(pth->binned_reio_num<2,
         pth->error_message,
         "current implementation of binned reio requires at least two bin centers");

    /* check that this input can be interpreted by the code */
    for (bin=1; bin<pth->binned_reio_num; bin++) {
      class_test(pth->binned_reio_z[bin]<pth->binned_reio_z[bin],
     pth->error_message,
     "value of reionization bin centers z_i expected to be passed in growing order");
    }

    /* the code will not only copy here the "bin centers" passed in
       input. It will add an initial and final value for (z,xe).   
       First, fill all entries except the first and the last */

    for (bin=1; bin<preio->reio_num_z-1; bin++) {
      preio->reionization_parameters[preio->index_reio_first_z+bin] = pth->binned_reio_z[bin-1]; 
      preio->reionization_parameters[preio->index_reio_first_xe+bin] = pth->binned_reio_xe[bin-1]; 
    }


    /* find largest value of z in the array. We choose to define it as
       z_(i_max) + (the distance between z_(i_max) and z_(i_max-1)). E.g. if
       the bins are in 10,12,14, the largest z will be 16. */
    preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1] = 
      2.*preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-2]
      -preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-3];

    /* copy this value in reio_start */
    preio->reionization_parameters[preio->index_reio_start] = preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1];

    /* check it's not too big */
    class_test(preio->reionization_parameters[preio->index_reio_start] > ppr->reionization_z_start_max,
         pth->error_message,
         "starting redshift for reionization = %e, reionization_z_start_max = %e, you must change the binning or increase reionization_z_start_max",
         preio->reionization_parameters[preio->index_reio_start],
         ppr->reionization_z_start_max);
    
    /* find smallest value of z in the array. We choose
       to define it as z_0 - (the distance between z_1 and z_0). E.g. if
       the bins are in 10,12,14, the stop redshift will be 8. */

    preio->reionization_parameters[preio->index_reio_first_z] = 
      2.*preio->reionization_parameters[preio->index_reio_first_z+1]
      -preio->reionization_parameters[preio->index_reio_first_z+2];

    /* check it's not too small */
    class_test(preio->reionization_parameters[preio->index_reio_first_z] < 0,
         pth->error_message,
         "final redshift for reionization = %e, you must change the binning or redefine the way in which the code extrapolates below the first value of z_i",preio->reionization_parameters[preio->index_reio_first_z]);

    /* infer xe before reio */
    class_call(thermodynamics_get_xe_before_reionization(ppr,
               pth,
               preco,
               preio->reionization_parameters[preio->index_reio_first_z+preio->reio_num_z-1],
               &(preio->reionization_parameters[preio->index_reio_first_xe+preio->reio_num_z-1])),
         pth->error_message,
         pth->error_message);

    /* infer xe after reio */
    preio->reionization_parameters[preio->index_reio_first_xe] = 1. + pth->YHe/(_not4_*(1.-pth->YHe));    /* xe_after_reio: H + singly ionized He (note: segmentation fault impossible, checked before that denominator is non-zero) */

    /* pass step sharpness parameter */
    preio->reionization_parameters[preio->index_reio_step_sharpness] = pth->binned_reio_step_sharpness;

    /* fill reionization table */
    class_call(thermodynamics_reionization_sample(ppr,pba,pth,preco,preio,pvecback),
         pth->error_message,
         pth->error_message);
    
    pth->tau_reio=preio->reionization_optical_depth;
    
    return _SUCCESS_;
    
  }

  class_test(0 == 0,
       pth->error_message,
       "value of reio_z_or_tau=%d unclear",pth->reio_z_or_tau);

}

/** 
 * For fixed input reionization parameters, this routine computes the
 * reionization history and fills the reionization table.
 *
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background structure
 * @param pth Input : pointer to thermo structure
 * @param preco Input : pointer to filled recombination structure
 * @param preio Input/Output: pointer to reionization structure (to be filled)
 * @param pvecback   Input: vector of background quantitites (used as workspace: must be already allocated, with format short_info or larger, but does not need to be filled)
 * @return the error status
 */

int thermodynamics_reionization_sample(
             struct precision * ppr,
             struct background * pba,
             struct thermo * pth,
             struct recombination * preco,
             struct reionization * preio,
             double * pvecback
             ) {

  /** Summary: */

  /** - define local variables */

  /* a growing table (since the number of redshift steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* current vector of values related to reionization */
  double * reio_vector;
  /* running index inside thermodynamics table */
  int i;
  int number_of_redshifts;
  /* values of z, dz, X_e */
  double dz,dz_max;
  double z,z_next;
  double xe,xe_next;
  double dkappadz,dkappadz_next;
  double Tb,Tba2,Yp,dTdz,opacity,mu;
  double dkappadtau,dkappadtau_next;

  Yp = pth->YHe;

  /** (a) allocate vector of values related to reionization */
  class_alloc(reio_vector,preio->re_size*sizeof(double),pth->error_message);

  /** (b) create a growTable with gt_init() */
  class_call(gt_init(&gTable),
       gTable.error_message,
       pth->error_message);

  /** (c) first line is taken from thermodynamics table, just before reionization starts */

  /** - look where to start in current thermodynamics table */
  i=0;
  while (preco->recombination_table[i*preco->re_size+preco->index_re_z] < preio->reionization_parameters[preio->index_reio_start]) {
    i++;
    class_test(i == ppr->recfast_Nz0,
         pth->error_message,
         "reionization_z_start_max = %e > largest redshift in thermodynamics table",ppr->reionization_z_start_max);
  }

  /** - get redshift */
  z=preco->recombination_table[i*preco->re_size+preco->index_re_z];
  reio_vector[preio->index_re_z]=z;
  preio->index_reco_when_reio_start=i;

  /** - get \f$ X_e \f$ */
  class_call(thermodynamics_reionization_function(z,pth,preio,&xe),
       pth->error_message,
       pth->error_message);

  reio_vector[preio->index_re_xe] = xe;

  /** - get \f$ d kappa / d z = (d kappa / d tau) * (d tau / d z) = - (d kappa / d tau) / H \f$ */
  class_call(background_functions(pba,1./(1.+z),pba->short_info,pvecback),
       pba->error_message,
       pth->error_message);

  reio_vector[preio->index_re_dkappadtau] = (1.+z) * (1.+z) * pth->n_e * xe * _sigma_ * _Mpc_over_m_;

  class_test(pvecback[pba->index_bg_H] == 0.,
       pth->error_message,
       "stop to avoid division by zero");

  reio_vector[preio->index_re_dkappadz] = reio_vector[preio->index_re_dkappadtau] / pvecback[pba->index_bg_H];

  dkappadz = reio_vector[preio->index_re_dkappadz];
  dkappadtau = reio_vector[preio->index_re_dkappadtau];

  /** - get baryon temperature **/
  Tb = preco->recombination_table[i*preco->re_size+preco->index_re_Tb];
  reio_vector[preio->index_re_Tb] = Tb;

  /** - after recombination, Tb scales like (1+z)**2. Compute constant factor Tb/(1+z)**2. */
  Tba2 = Tb/(1+z)/(1+z);

  /** - get baryon sound speed */
  reio_vector[preio->index_re_cb2] = 5./3. * _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * Yp + xe * (1.-Yp)) * Tb;

  /** - store these values in growing table */
  class_call(gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)),
       gTable.error_message,
       pth->error_message);

  number_of_redshifts=1;

  /** (d) set the maximum step value (equal to the step in thermodynamics table) */
  dz_max=preco->recombination_table[i*preco->re_size+preco->index_re_z]
    -preco->recombination_table[(i-1)*preco->re_size+preco->index_re_z];

  /** (e) loop over redshift values in order to find values of z, x_e, kappa' (Tb and cb2 found later by integration). The sampling in z space is found here. */

  while (z > 0.) {

    /** - try default step */
    dz = dz_max;
    z_next=z-dz;
    if (z_next < 0.) z_next=0.;

    class_call(thermodynamics_reionization_function(z_next,pth,preio,&xe_next),
         pth->error_message,
         pth->error_message);

    class_call(background_functions(pba,1./(1.+z_next),pba->short_info,pvecback),
         pba->error_message,
         pth->error_message);

    class_test(pvecback[pba->index_bg_H] == 0.,
         pth->error_message,
         "stop to avoid division by zero");

    dkappadz_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_ / pvecback[pba->index_bg_H];

    dkappadtau_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_;

    class_test((dkappadz == 0.) || (dkappadtau == 0.),
         pth->error_message,
         "stop to avoid division by zero");

    /** - reduce step if necessary */
    while (((fabs(dkappadz_next-dkappadz)/dkappadz) > ppr->reionization_sampling) || 
     ((fabs(dkappadtau_next-dkappadtau)/dkappadtau) > ppr->reionization_sampling)) {

      dz*=0.9;

      class_test(dz < ppr->smallest_allowed_variation,
     pth->error_message,
     "integration step =%e < machine precision : leads either to numerical error or infinite loop",dz);

      z_next=z-dz;
      if (z_next < 0.) z_next=0.;

      class_call(thermodynamics_reionization_function(z_next,pth,preio,&xe_next),
     pth->error_message,
     pth->error_message);

      class_call(background_functions(pba,1./(1.+z_next),pba->short_info,pvecback),
     pba->error_message,
     pth->error_message);

      class_test(pvecback[pba->index_bg_H] == 0.,
     pth->error_message,
     "stop to avoid division by zero");

      dkappadz_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_ / pvecback[pba->index_bg_H];

      dkappadtau_next= (1.+z_next) * (1.+z_next) * pth->n_e * xe_next * _sigma_ * _Mpc_over_m_;
    }

    /** - get \f$ z, X_e, d kappa / d z \f$ and store in growing table */
    z=z_next;
    xe=xe_next;
    dkappadz=dkappadz_next;
    dkappadtau= dkappadtau_next;

    class_test((dkappadz == 0.) || (dkappadtau == 0.),
         pth->error_message,
         "dkappadz=%e, dkappadtau=%e, stop to avoid division by zero",dkappadz,dkappadtau);

    reio_vector[preio->index_re_z] = z;   
    reio_vector[preio->index_re_xe] = xe;
    reio_vector[preio->index_re_dkappadz] = dkappadz;
    reio_vector[preio->index_re_dkappadtau] = dkappadz * pvecback[pba->index_bg_H];

    class_call(gt_add(&gTable,_GT_END_,(void *) reio_vector,sizeof(double)*(preio->re_size)),
         gTable.error_message,
         pth->error_message);

    number_of_redshifts++;
  }
  
  /** (f) allocate reionization_table with correct size */
  class_alloc(preio->reionization_table,preio->re_size*number_of_redshifts*sizeof(double),pth->error_message);

  preio->rt_size=number_of_redshifts;

  /** (g) retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
       gTable.error_message,
       pth->error_message);

  /** (h) copy growTable to reionization_temporary_table (invert order of lines, so that redshift is growing, like in recombination table) */
  for (i=0; i < preio->rt_size; i++) {
    memcopy_result = memcpy(preio->reionization_table+i*preio->re_size,pData+(preio->rt_size-i-1)*preio->re_size,preio->re_size*sizeof(double));
    class_test(memcopy_result != preio->reionization_table+i*preio->re_size,
         pth->error_message,
         "cannot copy data back to reionization_temporary_table");

  }

  /** (i) free the growTable with gt_free() , free vector of reionization variables */
  class_call(gt_free(&gTable),
       gTable.error_message,
       pth->error_message);
  
  free(reio_vector);

  /** (j) another loop on z, to integrate equation for Tb and to compute cb2 */
  for (i=preio->rt_size-1; i >0 ; i--) {

    z = preio->reionization_table[i*preio->re_size+preio->index_re_z];

    class_call(background_functions(pba,1./(1.+z),pba->normal_info,pvecback),
         pba->error_message,
         pth->error_message);
    
    dz = (preio->reionization_table[i*preio->re_size+preio->index_re_z]-preio->reionization_table[(i-1)*preio->re_size+preio->index_re_z]);
    
    opacity = (1.+z) * (1.+z) * pth->n_e 
      * preio->reionization_table[i*preio->re_size+preio->index_re_xe] * _sigma_ * _Mpc_over_m_;

    mu = _m_H_/(1. + (1./_not4_ - 1.) * pth->YHe + preio->reionization_table[i*preio->re_size+preio->index_re_xe] * (1.-pth->YHe));

    /** - derivative of baryon temperature */

    dTdz=2./(1+z)*preio->reionization_table[i*preio->re_size+preio->index_re_Tb]
      -2.*mu/_m_e_*4.*pvecback[pba->index_bg_rho_g]/3./pvecback[pba->index_bg_rho_b]*opacity*(pba->T_cmb * (1.+z)-preio->reionization_table[i*preio->re_size+preio->index_re_Tb])/pvecback[pba->index_bg_H];

    /** - increment baryon temperature */

    preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb] =
      preio->reionization_table[i*preio->re_size+preio->index_re_Tb]-dTdz*dz;

    /** - get baryon sound speed */

    preio->reionization_table[(i-1)*preio->re_size+preio->index_re_cb2] = _k_B_/ ( _c_ * _c_ * mu)
      * preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb]
      *(1.+(1+z)/3.*dTdz/preio->reionization_table[(i-1)*preio->re_size+preio->index_re_Tb]);

  }

  /** - spline \f$ d tau / dz \f$ with respect to z in view of integrating for optical depth */
  class_call(array_spline(preio->reionization_table,
        preio->re_size,
        preio->rt_size,
        preio->index_re_z,
        preio->index_re_dkappadz,
        preio->index_re_d3kappadz3,
        _SPLINE_EST_DERIV_,
        pth->error_message),
       pth->error_message,
       pth->error_message);
  
  /** - integrate for optical depth */
  class_call(array_integrate_all_spline(preio->reionization_table,
          preio->re_size,
          preio->rt_size,
          preio->index_re_z,
          preio->index_re_dkappadz,
          preio->index_re_d3kappadz3,
          &(preio->reionization_optical_depth),
          pth->error_message),
       pth->error_message,
       pth->error_message);

  return _SUCCESS_;

}

/** 
 * Integrate thermodynamics with your favorite recombination code.
 *
 */

int thermodynamics_recombination(
         struct precision * ppr,
         struct background * pba,
         struct thermo * pth,
         struct recombination * preco,
         double * pvecback
         ) {

  if (pth->recombination==hyrec) {

      class_call(thermodynamics_recombination_with_hyrec(ppr,pba,pth,preco,pvecback),
       pth->error_message,
       pth->error_message);

  }
  
  if (pth->recombination==recfast) {

    class_call(thermodynamics_recombination_with_recfast(ppr,pba,pth,preco,pvecback),
         pth->error_message,
         pth->error_message);

  }

  return _SUCCESS_;

}

/** 
 * Integrate thermodynamics with HyRec.
 *
 * Integrate thermodynamics with HyRec, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_recombination(), from thermodynamics_init().
 *
 *************************************************************************************************
 *                 HYREC: Hydrogen and Helium Recombination Code                                 *
 *         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              *
 ************************************************************************************************* 
 *
 * @param ppr      Input: pointer to precision structure
 * @param pba      Input: pointer to background structure
 * @param pth      Input: pointer to thermodynamics structure
 * @param preco    Ouput: pointer to recombination structure
 * @param pvecback Input: pointer to an allocated (but empty) vector of background variables
 *
 */

int thermodynamics_recombination_with_hyrec(
              struct precision * ppr,
              struct background * pba,
              struct thermo * pth,
              struct recombination * preco,
              double * pvecback
              ) {

#ifdef HYREC

  REC_COSMOPARAMS param;
  HRATEEFF rate_table;
  TWO_PHOTON_PARAMS twog_params;
  double *xe_output, *Tm_output;
  int i,j,l,Nz,b;
  double z, xe, Tm, Hz;
  FILE *fA;
  FILE *fR;
  double L2s1s_current;
  void * buffer;
  int buf_size;

  /** - Fill hyrec parameter structure */
  
  param.T0 = pba->T_cmb;
  param.obh2 = pba->Omega0_b*pba->h*pba->h;
  param.omh2 = (pba->Omega0_b+pba->Omega0_cdm+pba->Omega0_ncdm_tot)*pba->h*pba->h;
  param.okh2 = pba->Omega0_k*pba->h*pba->h;
  param.odeh2 = (pba->Omega0_lambda+pba->Omega0_fld)*pba->h*pba->h;
  param.w0 = pba->w0_fld;
  param.wa = pba->wa_fld;
  param.Y = pth->YHe;
  param.Nnueff = pba->Neff;
  param.nH0 = 11.223846333047*param.obh2*(1.-param.Y);  /* number density of hudrogen today in m-3 */
  param.fHe = param.Y/(1-param.Y)/3.97153;              /* abundance of helium by number */
  param.zstart = ppr->recfast_z_initial; /* Redshift range */
  param.zend = 0.;
  param.dlna = 8.49e-5;
  param.nz = (long) floor(2+log((1.+param.zstart)/(1.+param.zend))/param.dlna);  

  /** - Build effective rate tables */

  /* allocate contiguous memory zone */

  buf_size = (2*NTR+NTM+2*NTR*NTM+2*param.nz)*sizeof(double) + 2*sizeof(double*);

  class_alloc(buffer,
        buf_size,
        pth->error_message);
  
  /* distribute addresses for each table */

  rate_table.logTR_tab = (double*)buffer;
  rate_table.TM_TR_tab = (double*)(rate_table.logTR_tab + NTR);
  rate_table.logAlpha_tab[0] = (double**)(rate_table.TM_TR_tab+NTM);
  rate_table.logAlpha_tab[1] = (double**)(rate_table.logAlpha_tab[0]+NTM);
  rate_table.logAlpha_tab[0][0] = (double*)(rate_table.logAlpha_tab[1]+NTM);
  for (j=1;j<NTM;j++) {
    rate_table.logAlpha_tab[0][j] = (double*)(rate_table.logAlpha_tab[0][j-1]+NTR);
  }
  rate_table.logAlpha_tab[1][0] = (double*)(rate_table.logAlpha_tab[0][NTM-1]+NTR);
  for (j=1;j<NTM;j++) {
    rate_table.logAlpha_tab[1][j] = (double*)(rate_table.logAlpha_tab[1][j-1]+NTR);
  }
  rate_table.logR2p2s_tab = (double*)(rate_table.logAlpha_tab[1][NTM-1]+NTR);

  xe_output = (double*)(rate_table.logR2p2s_tab+NTR);
  Tm_output = (double*)(rate_table.logR2p2s_tab+param.nz);

  /* store sampled values of temperatures */

  for (i = 0; i < NTR; i++)
    rate_table.logTR_tab[i] = log(TR_MIN) + i * (log(TR_MAX)-log(TR_MIN))/(NTR-1.);
  for (i = 0; i < NTM; i++)
    rate_table.TM_TR_tab[i] = TM_TR_MIN + i * (TM_TR_MAX-TM_TR_MIN)/(NTM-1.);

  rate_table.DlogTR = rate_table.logTR_tab[1] - rate_table.logTR_tab[0];
  rate_table.DTM_TR = rate_table.TM_TR_tab[1] - rate_table.TM_TR_tab[0];  

  /* read in file */

  class_open(fA,ppr->hyrec_Alpha_inf_file, "r",pth->error_message);
  class_open(fR,ppr->hyrec_R_inf_file, "r",pth->error_message);
  
  for (i = 0; i < NTR; i++) {  
    for (j = 0; j < NTM; j++) {  
      for (l = 0; l <= 1; l++) {
  if (fscanf(fA, "%le", &(rate_table.logAlpha_tab[l][j][i])) != 1)
    class_stop(pth->error_message,"Error reading hyrec data file %s",ppr->hyrec_Alpha_inf_file);
  rate_table.logAlpha_tab[l][j][i] = log(rate_table.logAlpha_tab[l][j][i]);
      }
    }
    
    if (fscanf(fR, "%le", &(rate_table.logR2p2s_tab[i])) !=1)
      class_stop(pth->error_message,"Error reading hyrec data file %s",ppr->hyrec_R_inf_file);
    rate_table.logR2p2s_tab[i] = log(rate_table.logR2p2s_tab[i]);
    
  }
  fclose(fA);
  fclose(fR);
 
  /* Read two-photon rate tables */
 
  class_open(fA,ppr->hyrec_two_photon_tables_file, "r",pth->error_message);  

  for (b = 0; b < NVIRT; b++) { 
    if ((fscanf(fA, "%le", &(twog_params.Eb_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A1s_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A2s_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A3s3d_tab[b])) != 1) ||
        (fscanf(fA, "%le", &(twog_params.A4s4d_tab[b])) != 1))
      class_stop(pth->error_message,"Error reading hyrec data file %s",ppr->hyrec_two_photon_tables_file);
   }  

   fclose(fA); 

   /** - Normalize 2s--1s differential decay rate to L2s1s (can be set by user in hydrogen.h) */
   L2s1s_current = 0.;
   for (b = 0; b < NSUBLYA; b++) L2s1s_current += twog_params.A2s_tab[b];
   for (b = 0; b < NSUBLYA; b++) twog_params.A2s_tab[b] *= L2s1s/L2s1s_current;

   /*  In CLASS, we have neutralized the switches for the various
       effects considered in Hirata (2008), keeping the full
       calculation as a default; but you could restore their
       functionality by copying a few lines from hyrec/hyrec.c to
       here */

   /** - Compute the recombination history by calling a function in hyrec (no CLASS-like error management here) */
 
   if (pth->thermodynamics_verbose > 0)
     printf(" -> calling HyRec version %s,\n",HYREC_VERSION);
 
   rec_build_history(&param, &rate_table, &twog_params, xe_output, Tm_output);

   if (pth->thermodynamics_verbose > 0)
     printf("    by Y. Ali-Haïmoud & C. Hirata\n");

  /** - fill a few parameters in preco and pth */
  
  Nz=ppr->recfast_Nz0;

  preco->rt_size = Nz;
  preco->H0 = pba->H0 * _c_ / _Mpc_over_m_; 
  /* preco->H0 in inverse seconds (while pba->H0 is [H0/c] in inverse Mpcs) */
  preco->Nnow = 3.*preco->H0*preco->H0*pba->Omega0_b*(1.-pth->YHe)/(8.*_PI_*_G_*_m_H_);

  pth->n_e=preco->Nnow;

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) and fill it */

  
  class_alloc(preco->recombination_table,preco->re_size*preco->rt_size*sizeof(double),pth->error_message);

  for(i=0; i <Nz; i++) {

    /** -> get redshift, corresponding results from hyrec, and background quantities */

    z = param.zstart * (1. - (double)(i+1) / (double)Nz);

    /* get (xe,Tm) by interpolating in pre-computed tables */

    class_call(array_interpolate_cubic_equal(-log(1.+param.zstart), 
               param.dlna, 
               xe_output, 
               param.nz, 
               -log(1.+z),
               &xe,
               pth->error_message),
         pth->error_message,
         pth->error_message);
 
    class_call(array_interpolate_cubic_equal(-log(1.+param.zstart), 
               param.dlna, 
               Tm_output, 
               param.nz, 
               -log(1.+z),
               &Tm,
               pth->error_message),
         pth->error_message,
         pth->error_message);
    
    class_call(background_functions(pba,pba->a_today/(1.+z),pba->short_info,pvecback),
         pba->error_message,
         pth->error_message);
  
    /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
    Hz=pvecback[pba->index_bg_H] * _c_ / _Mpc_over_m_;

    /** -> store the results in the table */

    /* results are obtained in order of decreasing z, and stored in order of growing z */

    /* redshift */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z)=z;

    /* ionization fraction */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe)=xe;

    /* Tb */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb)=Tm;
    
    /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz) 
       with (1+z)dlnTb/dz= - [dlnTb/dlna] */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2)
      = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * pth->YHe + xe * (1.-pth->YHe)) * Tm * (1. - rec_dTmdlna(xe, Tm, pba->T_cmb*(1.+z), Hz, param.fHe) / Tm / 3.);

    /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
      = (1.+z) * (1.+z) * preco->Nnow * xe * _sigma_ * _Mpc_over_m_;
    
  }
  
  /* Cleanup */

  free(buffer);
  
#else

  class_stop(pth->error_message,
       "you compiled without including the HyRec code, and now whish to use it. Either set the input parameter 'recombination' to something else than 'HyRec', or recompile after setting in the Makefile the appropriate path HYREC=... ");

#endif

  return _SUCCESS_;
}

/** 
 * Integrate thermodynamics with RECFAST.
 *
 * Integrate thermodynamics with RECFAST, allocate and fill the part
 * of the thermodynamics interpolation table (the rest is filled in
 * thermodynamics_init()). Called once by
 * thermodynamics_recombination, from thermodynamics_init().
 *
 *******************************************************************************
 * RECFAST is an integrator for Cosmic Recombination of Hydrogen and Helium,   *
 * developed by Douglas Scott (dscott@astro.ubc.ca)                            *
 * based on calculations in the paper Seager, Sasselov & Scott                 *
 * (ApJ, 523, L1, 1999).                                                       *
 * and "fudge" updates in Wong, Moss & Scott (2008).                           *
 *                                                                             *
 * Permission to use, copy, modify and distribute without fee or royalty at    *
 * any tier, this software and its documentation, for any purpose and without  *
 * fee or royalty is hereby granted, provided that you agree to comply with    *
 * the following copyright notice and statements, including the disclaimer,    *
 * and that the same appear on ALL copies of the software and documentation,   *
 * including modifications that you make for internal use or for distribution: *
 *                                                                             *
 * Copyright 1999-2010 by University of British Columbia.  All rights reserved.*
 *                                                                             *
 * THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO                      *
 * REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.                          *
 * BY WAY OF EXAMPLE, BUT NOT LIMITATION,                                      *
 * U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF                            *
 * MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT               *
 * THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE         *
 * ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.            *
 *******************************************************************************
 *
 * Version 1.5: includes extra fitting function from
 *              Rubino-Martin et al. arXiv:0910.4383v1 [astro-ph.CO]
 *
 * @param ppr      Input: pointer to precision structure
 * @param pba      Input: pointer to background structure
 * @param pth      Input: pointer to thermodynamics structure
 * @param preco    Ouput: pointer to recombination structure
 * @param pvecback Input: pointer to an allocated (but empty) vector of background variables
 * @return the error status
 */

int thermodynamics_recombination_with_recfast(
                struct precision * ppr,
                struct background * pba,
                struct thermo * pth,
                struct recombination * preco,
                double * pvecback
                ) {

  /** Summary: */

  /** - define local variables */

  /* vector of variables to be integrated: x_H, x_He, Tmat */
  double y[3],dy[3];

  /* other recfast variables */
  double OmegaB,Yp,zinitial,x_He0,x0;
  double x_H0=0.;
  double z,mu_H,n,Lalpha,Lalpha_He,DeltaB,DeltaB_He,mu_T;
  double zstart,zend,rhs;
  int i,Nz;

  /* introduced by JL for smoothing the various steps */
  double x0_previous,x0_new,s,weight;

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;

  /* contains all fixed parameters which should be passed to thermodynamics_derivs_with_recfast */
  struct thermodynamics_parameters_and_workspace tpaw;

  /** - allocate memory for thermodynamics interpolation tables (size known in advance) */
  preco->rt_size = ppr->recfast_Nz0;
  class_alloc(preco->recombination_table,preco->re_size*preco->rt_size*sizeof(double),pth->error_message);

  /** - initialize generic integrator with initialize_generic_integrator() */
  class_call(initialize_generic_integrator(_RECFAST_INTEG_SIZE_, &gi),
       gi.error_message,
       pth->error_message);
  
  /** - read a few precision/cosmological parameters */

  /* Nz */
  Nz=ppr->recfast_Nz0;

  /* preco->H0 is H0 in inverse seconds (while pba->H0 is [H0/c] in inverse Mpcs) */
  preco->H0 = pba->H0 * _c_ / _Mpc_over_m_;

  /* Omega_b */
  OmegaB = pba->Omega0_b;

  /* Yp */
  Yp = pth->YHe;

  /* Tnow */
  preco->Tnow = pba->T_cmb;

  /* z_initial */
  zinitial=ppr->recfast_z_initial;

  /* H_frac */ 
  preco->H_frac = ppr->recfast_H_frac;

  /* H fudging */
  class_test((ppr->recfast_Hswitch != _TRUE_) && (ppr->recfast_Hswitch != _FALSE_),
       pth->error_message,
       "RECFAST error: unknown H fudging scheme");
  preco->fu = ppr->recfast_fudge_H;
  if (ppr->recfast_Hswitch == _TRUE_) 
    preco->fu += ppr->recfast_delta_fudge_H;

  /* He fudging */
  class_test((ppr->recfast_Heswitch < 0) || (ppr->recfast_Heswitch > 6),
       pth->error_message,
       "RECFAST error: unknown He fudging scheme");

  /* related quantities */ 
  z=zinitial;
  mu_H = 1./(1.-Yp);
  mu_T = _not4_ /(_not4_ - (_not4_-1.)*Yp); /* recfast 1.4*/
  preco->fHe = Yp/(_not4_ *(1.-Yp)); /* recfast 1.4 */
  preco->Nnow = 3.*preco->H0*preco->H0*OmegaB/(8.*_PI_*_G_*mu_H*_m_H_);
  pth->n_e = preco->Nnow;

  /* quantities related to constants defined in thermodynamics.h */
  n = preco->Nnow * pow((1.+z),3);
  Lalpha = 1./_L_H_alpha_;
  Lalpha_He = 1./_L_He_2p_;
  DeltaB = _h_P_*_c_*(_L_H_ion_-_L_H_alpha_);
  preco->CDB = DeltaB/_k_B_;
  DeltaB_He = _h_P_*_c_*(_L_He1_ion_-_L_He_2s_);
  preco->CDB_He = DeltaB_He/_k_B_;
  preco->CB1 = _h_P_*_c_*_L_H_ion_/_k_B_;
  preco->CB1_He1 = _h_P_*_c_*_L_He1_ion_/_k_B_;
  preco->CB1_He2 = _h_P_*_c_*_L_He2_ion_/_k_B_;
  preco->CR = 2.*_PI_*(_m_e_/_h_P_)*(_k_B_/_h_P_);
  preco->CK = pow(Lalpha,3)/(8.*_PI_);
  preco->CK_He = pow(Lalpha_He,3)/(8.*_PI_);
  preco->CL = _c_*_h_P_/(_k_B_*Lalpha);
  preco->CL_He = _c_*_h_P_/(_k_B_/_L_He_2s_);
  preco->CT = (8./3.) * (_sigma_/(_m_e_*_c_)) * 
    (8.*pow(_PI_,5)*pow(_k_B_,4)/ 15./ pow(_h_P_,3)/pow(_c_,3));

  preco->Bfact = _h_P_*_c_*(_L_He_2p_-_L_He_2s_)/_k_B_;

  /** - define the fields of the 'thermodynamics parameter and workspace' structure */
  tpaw.pba = pba;
  tpaw.ppr = ppr;
  tpaw.preco = preco;
  tpaw.pvecback = pvecback;

  /** - impose initial conditions at early times */

  class_test(zinitial < ppr->recfast_z_He_3,
       pth->error_message,
       "increase zinitial, otherwise should get initial conditions from recfast's get_init routine (less precise anyway)");
  
  y[0] = 1.;
  y[1] = 1.;
  x0 = 1.+2.*preco->fHe;
  y[2] = preco->Tnow*(1.+z);

  /** - loop over redshift steps Nz; integrate over each step with
      generic_integrator(), store the results in the table using
      thermodynamics_derivs_with_recfast()*/

  for(i=0; i <Nz; i++) {

    zstart = zinitial * (1. - (double)i / (double)Nz);
    zend   = zinitial * (1. - (double)(i+1) / (double)Nz);
    z = zend;

    /** -> first approximation: H and Helium fully ionized */

    if (z > ppr->recfast_z_He_1+ppr->recfast_delta_z_He_1) {
      x_H0 = 1.;
      x_He0 = 1.;
      x0 = 1.+2.*preco->fHe;
      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** -> second approximation: first Helium recombination (analytic approximation) */

    else if (z > ppr->recfast_z_He_2+ppr->recfast_delta_z_He_2) {
      x_H0 = 1.;
      x_He0 = 1.;
      
      rhs = exp( 1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He2/(preco->Tnow*(1.+z)) ) / preco->Nnow;
      
      /* smoothed transition */
      if (z > ppr->recfast_z_He_1-ppr->recfast_delta_z_He_1) {
  x0_previous = 1.+2.*preco->fHe;
  x0_new = 0.5*(sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe));
  
  /* get s from -1 to 1 */
  s = (ppr->recfast_z_He_1-z)/ppr->recfast_delta_z_He_1;
  /* infer f1(s) = smooth function interpolating from 0 to 1 */
  weight = f1(s);
  
  x0 = weight*x0_new+(1.-weight)*x0_previous;
      }
      /* transition finished */
      else {
  x0 = 0.5*(sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe));
      }
      
      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** -> third approximation: first Helium recombination completed */

    else if (z > ppr->recfast_z_He_3+ppr->recfast_delta_z_He_3) {
      x_H0 = 1.;
      x_He0 = 1.;

      /* smoothed transition */
      if (z > ppr->recfast_z_He_2-ppr->recfast_delta_z_He_2) {
  rhs = exp( 1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He2/(preco->Tnow*(1.+z)) ) / preco->Nnow;
  x0_previous = 0.5*(sqrt(pow((rhs-1.-preco->fHe),2) + 4.*(1.+2.*preco->fHe)*rhs) - (rhs-1.-preco->fHe));
  x0_new = 1. + preco->fHe;
  /* get s from -1 to 1 */
  s = (ppr->recfast_z_He_2-z)/ppr->recfast_delta_z_He_2;
  /* infer f1(s) = smooth function interpolating from 0 to 1 */
  weight = f1(s);
  
  x0 = weight*x0_new+(1.-weight)*x0_previous;
  
      }
      /* transition finished */
      else {
  x0 = 1.+preco->fHe;
      }
      
      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }
  
    /** -> fourth approximation: second Helium recombination starts (analytic approximation) */

    else if (y[1] > ppr->recfast_x_He0_trigger) {
      x_H0 = 1.;

      rhs = 4.*exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He1/(preco->Tnow*(1.+z)))/preco->Nnow;
      x_He0 = 0.5*(sqrt(pow((rhs-1.),2) + 4.*(1.+preco->fHe)*rhs )- (rhs-1.));
      
      /* smoothed transition */
      if (z > ppr->recfast_z_He_3-ppr->recfast_delta_z_He_3) {
  x0_previous = 1. + preco->fHe;
  x0_new = x_He0;
  /* get s from -1 to 1 */
  s = (ppr->recfast_z_He_3-z)/ppr->recfast_delta_z_He_3;
  /* infer f1(x) = smooth function interpolating from 0 to 1 */
  weight = f1(s);
  
  x0 = weight*x0_new+(1.-weight)*x0_previous;
      }
      /* transition finished */
      else {
  x0 = x_He0;
      }
      
      x_He0 = (x0-1.)/preco->fHe;
      y[0] = x_H0;
      y[1] = x_He0;
      y[2] = preco->Tnow*(1.+z);
    }

    /** -> fifth approximation: second Helium recombination (full
           evolution for Helium), H recombination starts (analytic
           approximation) */

    else if (y[0] > ppr->recfast_x_H0_trigger) {

      rhs = exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1/(preco->Tnow*(1.+z)))/preco->Nnow;
      x_H0 = 0.5*(sqrt(pow(rhs,2)+4.*rhs) - rhs);
      
      class_call(generic_integrator(thermodynamics_derivs_with_recfast,
            zstart,
            zend,
            y,
            &tpaw,
            ppr->tol_thermo_integration,
            ppr->smallest_allowed_variation,
            &gi),
     gi.error_message,
     pth->error_message);
      
      y[0] = x_H0;
      
      /* smoothed transition */
      if (ppr->recfast_x_He0_trigger - y[1] < ppr->recfast_x_He0_trigger_delta) {
  rhs = 4.*exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1_He1/(preco->Tnow*(1.+z)))/preco->Nnow;
  x0_previous = 0.5*(sqrt(pow((rhs-1.),2) + 4.*(1.+preco->fHe)*rhs )- (rhs-1.));
  x0_new = y[0] + preco->fHe*y[1];
  /* get s from 0 to 1 */
  s = (ppr->recfast_x_He0_trigger - y[1])/ppr->recfast_x_He0_trigger_delta;
  /* infer f2(x) = smooth function interpolating from 0 to 1 */
  weight = f2(s);
  
  x0 = weight*x0_new+(1.-weight)*x0_previous;
      }
      /* transition finished */
      else {
  x0 = y[0] + preco->fHe*y[1];
      }
                
    }

    /** -> last case: full evolution for H and Helium */
      
    else {
        
      /* quantities used for smoothed transition */
      if (ppr->recfast_x_H0_trigger - y[0] < ppr->recfast_x_H0_trigger_delta) {
  rhs = exp(1.5*log(preco->CR*preco->Tnow/(1.+z)) - preco->CB1/(preco->Tnow*(1.+z)))/preco->Nnow;
  x_H0 = 0.5*(sqrt(pow(rhs,2)+4.*rhs) - rhs);
      }
      
      class_call(generic_integrator(thermodynamics_derivs_with_recfast,
            zstart,
            zend,
            y,
            &tpaw,
            ppr->tol_thermo_integration,
            ppr->smallest_allowed_variation,
            &gi),
     gi.error_message,
     pth->error_message);
      
      /* smoothed transition */
      if (ppr->recfast_x_H0_trigger - y[0] < ppr->recfast_x_H0_trigger_delta) {
  /* get s from 0 to 1 */
  s = (ppr->recfast_x_H0_trigger - y[0])/ppr->recfast_x_H0_trigger_delta; 
  /* infer f2(s) = smooth function interpolating from 0 to 1 */
  weight = f2(s);
  
  x0 = weight*y[0]+(1.-weight)*x_H0 + preco->fHe*y[1];
  
      }
      /* transition finished */ 
      else {
  x0 = y[0] + preco->fHe*y[1];
      }
    }
    
    /** -> store the results in the table */
    /* results are obtained in order of decreasing z, and stored in order of growing z */

    /* redshift */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z)=zend;

    /* ionization fraction */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe)=x0;

    /* Tb */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb)=y[2];

    /* get dTb/dz=dy[2] */
    class_call(thermodynamics_derivs_with_recfast(zend, y, dy, &tpaw,pth->error_message),
         pth->error_message,
         pth->error_message);

    /* cb2 = (k_B/mu) Tb (1-1/3 dlnTb/dlna) = (k_B/mu) Tb (1+1/3 (1+z) dlnTb/dz) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2)
      = _k_B_ / ( _c_ * _c_ * _m_H_ ) * (1. + (1./_not4_ - 1.) * Yp + x0 * (1.-Yp)) * y[2] * (1. + (1.+zend) * dy[2] / y[2] / 3.);

    /* dkappa/dtau = a n_e x_e sigma_T = a^{-2} n_e(today) x_e sigma_T (in units of 1/Mpc) */
    *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau)
      = (1.+zend) * (1.+zend) * preco->Nnow * x0 * _sigma_ * _Mpc_over_m_;


    // *** MY MODIFICATIONS ***

    // *** Rayleigh scattering ***
    if (pth->has_rayleigh_scattering == _TRUE_) {

      /* x_HI, number of neutral H atoms per H nucleus. It is obtained as x_HI = 1 - x_HII = 1 - y[0] */
      double x_HI = 1 - y[0];
      
      /* x_HeI = n_HeI / (n_HI + n_HII), number of neutral He atoms per H nuclues. It is obtained as
      x_HeI = x_He - x_HeII - x_HeIII, where x_He=fHe is the total number of He nuclei per H nuclues. Since there
      are no fully ionised Helium nuclei around recombination, we set x_HeIII to zero. The quantity
      evolved by CLASS is y[1] = n_HeII/n_He = x_HeII*(n_H/n_He) = x_HeII/fHe, so that
      x_HeI = x_He - fHe*y[1], but x_He = fHe so that x_HeI = fHe * (1 - y[1]).  */
      double x_HeI = preco->fHe * (1 - y[1]);
      
      /* Relative strength of Rayleigh scattering on He compared to H (arXiv:1307.8148) */
      double R_He = 0.1;
      
      /* Rayleigh scattering rate, without frequency dependence. This is the same as for Thomson scattering,
      but using the neutral atoms densities instead of the free electron one. */
      *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_rayleigh_dkappadtau)
        = (1.+zend) * (1.+zend) * preco->Nnow * (x_HI + R_He*x_HeI) * _sigma_ * _Mpc_over_m_;
      
      /* Include the frequency dependence (see eq. 2.1 of arXiv:1307.8148) in the Rayleigh scattering rate,
      making sure to include the scaling due to the redshifting of photons. */
      double nu_ratio = pth->rayleigh_frequency/_NU_EFF_;
      double nu_4 = pow((1.+zend)*nu_ratio,4);
      double nu_6 = pow((1.+zend)*nu_ratio,6) * 638/243.;
      double nu_8 = pow((1.+zend)*nu_ratio,8) * 1299667/239196.;
      *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_rayleigh_dkappadtau) *= nu_4 + nu_6 + nu_8;
      
      /* Uncomment to modify the full dkappa */
      *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau) += 
        *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_rayleigh_dkappadtau);

    }
    
    // *** END OF MY MODIFICATIONS ***
    
    /* fprintf(stdout,"%e %e %e %e %e %e\n", */
    /*      *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_z), */
    /*      *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_xe), */
    /*      *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_Tb), */
    /*      (1.+zend) * dy[2], */
    /*      *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_cb2), */
    /*      *(preco->recombination_table+(Nz-i-1)*preco->re_size+preco->index_re_dkappadtau) */
    /*      ); */

  }

  /** - cleanup generic integrator with cleanup_generic_integrator() */

  class_call(cleanup_generic_integrator(&gi),
       gi.error_message,
       pth->error_message);
  
  return _SUCCESS_;
}

/**
 * Subroutine evaluating the derivative with respect to redshift of
 * thermodynamical quantities (from RECFAST version 1.4).
 *
 * Computes derivatives of the three variables to integrate: \f$ d x_H
 * / dz, d x_{He} / dz, d T_{mat} / dz \f$.
 * 
 * This is one of the few functions in the code which are passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed parameters and workspaces are passed through a generic
 *   pointer. Here, this pointer contains the precision, background
 *   and recombination structures, plus a background vector, but
 *   generic_integrator() doesn't know its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 *   usual to pth->error_message, but to a generic error_message
 *   passed in the list of arguments.
 *
 * @param z                        Input : redshift
 * @param y                        Input : vector of variable to integrate
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input : pointer to fixed parameters (e.g. indices) and workspace (already allocated)
 * @param error_message            Output: error message
 */

int thermodynamics_derivs_with_recfast(
               double z,
               double * y,
               double * dy,
               void * parameters_and_workspace,
               ErrorMsg error_message
               ) {

  /** Summary: */

  /** - define local variables */

  double x,n,n_He,Trad,Tmat,x_H,x_He,Hz,dHdz,epsilon;
  double Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz;
  double timeTh,timeH;
  double sq_0,sq_1;

  /* new in recfast 1.4: */
  double Rdown_trip,Rup_trip,tauHe_s,pHe_s,Doppler,gamma_2Ps,pb,qb,AHcon;
  double tauHe_t,pHe_t,CL_PSt,gamma_2Pt;
  double CfHe_t=0.;
  int Heflag;

  struct thermodynamics_parameters_and_workspace * ptpaw;
  struct precision * ppr;
  struct background * pba;
  struct recombination * preco;
  double * pvecback;

  ptpaw = parameters_and_workspace;
  ppr = ptpaw->ppr;
  pba = ptpaw->pba;
  preco = ptpaw->preco;
  pvecback = ptpaw->pvecback;

  x_H = y[0];
  x_He = y[1];
  x = x_H + preco->fHe * x_He;
  Tmat = y[2];

  n = preco->Nnow * (1.+z) * (1.+z) * (1.+z);
  n_He = preco->fHe * n;
  Trad = preco->Tnow * (1.+z);

  class_call(background_functions(pba,1./(1.+z),pba->short_info,pvecback),
       pba->error_message,
       error_message);
  
  /* Hz is H in inverse seconds (while pvecback returns [H0/c] in inverse Mpcs) */
  Hz=pvecback[pba->index_bg_H]* _c_ / _Mpc_over_m_;

  Rdown=1.e-19*_a_PPB_*pow((Tmat/1.e4),_b_PPB_)/(1.+_c_PPB_*pow((Tmat/1.e4),_d_PPB_));
  Rup = Rdown * pow((preco->CR*Tmat),1.5)*exp(-preco->CDB/Tmat);

  sq_0 = sqrt(Tmat/_T_0_);
  sq_1 = sqrt(Tmat/_T_1_);
  Rdown_He = _a_VF_/(sq_0 * pow((1.+sq_0),(1.-_b_VF_)) * pow((1. + sq_1),(1. + _b_VF_)));
  Rup_He = 4.*Rdown_He*pow((preco->CR*Tmat),1.5)*exp(-preco->CDB_He/Tmat);
  K = preco->CK/Hz;

  /* following is from recfast 1.5 */

  if (ppr->recfast_Hswitch == _TRUE_ )
    K *= 1.
      + ppr->recfast_AGauss1*exp(-pow((log(1.+z)-ppr->recfast_zGauss1)/ppr->recfast_wGauss1,2)) 
      + ppr->recfast_AGauss2*exp(-pow((log(1.+z)-ppr->recfast_zGauss2)/ppr->recfast_wGauss2,2));

  /* end of new recfast 1.5 piece */

  /* following is from recfast 1.4 */

  Rdown_trip = _a_trip_/(sq_0*pow((1.+sq_0),(1.-_b_trip_)) * pow((1.+sq_1),(1.+_b_trip_)));
  Rup_trip = Rdown_trip*exp(-_h_P_*_c_*_L_He2St_ion_/(_k_B_*Tmat))*pow(preco->CR*Tmat,1.5)*4./3.;

  if ((x_He < 5.e-9) || (x_He > ppr->recfast_x_He0_trigger2)) 
    Heflag = 0;
  else
    Heflag = ppr->recfast_Heswitch;

  if (Heflag == 0)
    K_He = preco->CK_He/Hz;
  else {
    tauHe_s = _A2P_s_*preco->CK_He*3.*n_He*(1.-x_He)/Hz;
    pHe_s = (1.-exp(-tauHe_s))/tauHe_s;
    K_He = 1./(_A2P_s_*pHe_s*3.*n_He*(1.-x_He));

    /*    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.99999)) { */
    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.9999999)) { /* threshold changed by Antony Lewis in 2008 to get smoother Helium */

      Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_c_*_c_);
      Doppler = _c_*_L_He_2p_*sqrt(Doppler);
      gamma_2Ps = 3.*_A2P_s_*preco->fHe*(1.-x_He)*_c_*_c_
  /(sqrt(_PI_)*_sigma_He_2Ps_*8.*_PI_*Doppler*(1.-x_H))
  /pow(_c_*_L_He_2p_,2);
      pb = 0.36;
      qb = ppr->recfast_fudge_He;
      AHcon = _A2P_s_/(1.+pb*pow(gamma_2Ps,qb));
      K_He=1./((_A2P_s_*pHe_s+AHcon)*3.*n_He*(1.-x_He));
    }

    if (Heflag >= 3) {
      tauHe_t = _A2P_t_*n_He*(1.-x_He)*3./(8.*_PI_*Hz*pow(_L_He_2Pt_,3));
      pHe_t = (1. - exp(-tauHe_t))/tauHe_t;
      CL_PSt = _h_P_*_c_*(_L_He_2Pt_ - _L_He_2St_)/_k_B_;
      if ((Heflag == 3) || (Heflag == 5) || (x_H >= 0.99999)) {
  CfHe_t = _A2P_t_*pHe_t*exp(-CL_PSt/Tmat);
  CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
      else {
  Doppler = 2.*_k_B_*Tmat/(_m_H_*_not4_*_c_*_c_);
  Doppler = _c_*_L_He_2Pt_*sqrt(Doppler);
  gamma_2Pt = 3.*_A2P_t_*preco->fHe*(1.-x_He)*_c_*_c_
    /(sqrt(_PI_)*_sigma_He_2Pt_*8.*_PI_*Doppler*(1.-x_H))
    /pow(_c_*_L_He_2Pt_,2);
  pb = 0.66;
  qb = 0.9;
  AHcon = _A2P_t_/(1.+pb*pow(gamma_2Pt,qb))/3.;
  CfHe_t = (_A2P_t_*pHe_t+AHcon)*exp(-CL_PSt/Tmat);
  CfHe_t = CfHe_t/(Rup_trip+CfHe_t);
      }
    }
  }

  /* end of new recfast 1.4 piece */

  timeTh=(1./(preco->CT*pow(Trad,4)))*(1.+x+preco->fHe)/x;
  timeH=2./(3.*preco->H0*pow(1.+z,1.5)); 

  if (x_H > ppr->recfast_x_H0_trigger)
    dy[0] = 0.;
  else {
    if (x_H > ppr->recfast_x_H0_trigger2) {
      dy[0] = (x*x_H*n*Rdown - Rup*(1.-x_H)*exp(-preco->CL/Tmat))/ (Hz*(1.+z));
    }
    else {
      dy[0] = ((x*x_H*n*Rdown - Rup*(1.-x_H)*exp(-preco->CL/Tmat)) *(1. + K*_Lambda_*n*(1.-x_H))) /(Hz*(1.+z)*(1./preco->fu+K*_Lambda_*n*(1.-x)/preco->fu +K*Rup*n*(1.-x)));

    }
  }

  if (x_He < 1.e-15) 
    dy[1]=0.;
  else {

    if (preco->Bfact/Tmat < 680.) 
      He_Boltz=exp(preco->Bfact/Tmat);
    else 
      He_Boltz=exp(680.);

    dy[1] = ((x*x_He*n*Rdown_He - Rup_He*(1.-x_He)*exp(-preco->CL_He/Tmat)) 
       *(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz))
      /(Hz*(1+z)* (1. + K_He*(_Lambda_He_+Rup_He)*n_He*(1.-x_He)*He_Boltz));

    /* following is from recfast 1.4 */

    if (Heflag >= 3)
      dy[1] = dy[1] + 
  (x*x_He*n*Rdown_trip
   - (1.-x_He)*3.*Rup_trip*exp(-_h_P_*_c_*_L_He_2St_/(_k_B_*Tmat)))
  *CfHe_t/(Hz*(1.+z));

    /* end of new recfast 1.4 piece */

  }

  if (timeTh < preco->H_frac*timeH) {
  /*   dy[2]=Tmat/(1.+z); */
    /* v 1.5: like in camb, add here a smoothing term as suggested by Adam Moss */
    dHdz=-pvecback[pba->index_bg_H_prime]/pvecback[pba->index_bg_H]/pba->a_today* _c_ / _Mpc_over_m_;
    epsilon = Hz * (1.+x+preco->fHe) / (preco->CT*pow(Trad,3)*x);
    dy[2] = preco->Tnow + epsilon*((1.+preco->fHe)/(1.+preco->fHe+x))*((dy[0]+preco->fHe*dy[1])/x)
      - epsilon* dHdz/Hz + 3.*epsilon/(1.+z);
    
  }
  else {
    dy[2]= preco->CT * pow(Trad,4) * x / (1.+x+preco->fHe) * (Tmat-Trad) / (Hz*(1.+z)) + 2.*Tmat/(1.+z);
  }

  return _SUCCESS_;
}      

/** 
 * This routine merges the two tables 'recombination_table' and
 * 'reionization_table' inside the table 'thermodynamics_table', and
 * frees the temporary structures 'recombination' and 'reionization'.
 *
 * @param ppr   Input : pointer to precision structure
 * @param pth   Input/Output : pointer to thermo structure
 * @param preco Input : pointer to filled recombination structure
 * @param preio Input : pointer to reionization structure
 * @return the error status
 */

int thermodynamics_merge_reco_and_reio(
               struct precision * ppr,
               struct thermo * pth,
               struct recombination * preco,
               struct reionization * preio
               ) {
  /** Summary: */

  /** - define local variables */

  int i,index_th,index_re;

  /** - first, a little check that the two tables match each other and can be merged */

  if (pth->reio_parametrization != reio_none) {
    class_test(preco->recombination_table[preio->index_reco_when_reio_start*preco->re_size+preco->index_re_z] !=
         preio->reionization_table[(preio->rt_size -1)*preio->re_size+preio->index_re_z],
         pth->error_message,
         "mismatch which should never happen");
  }

  /** - find number of redshift in full table = number in reco + number in reio - overlap */

  pth->tt_size = ppr->recfast_Nz0 + preio->rt_size - preio->index_reco_when_reio_start - 1;


  /** - allocate arrays in thermo structure */

  class_alloc(pth->z_table,pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->thermodynamics_table,pth->th_size*pth->tt_size*sizeof(double),pth->error_message);
  class_alloc(pth->d2thermodynamics_dz2_table,pth->th_size*pth->tt_size*sizeof(double),pth->error_message);  
  
  /** - fill these arrays */

  for (i=0; i < preio->rt_size; i++) {
    pth->z_table[i]=
      preio->reionization_table[i*preio->re_size+preio->index_re_z];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_xe]=
      preio->reionization_table[i*preio->re_size+preio->index_re_xe];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_dkappa]=
      preio->reionization_table[i*preio->re_size+preio->index_re_dkappadtau];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_Tb]=
      preio->reionization_table[i*preio->re_size+preio->index_re_Tb];
    pth->thermodynamics_table[i*pth->th_size+pth->index_th_cb2]=
      preio->reionization_table[i*preio->re_size+preio->index_re_cb2];
  }
  for (i=0; i < ppr->recfast_Nz0 - preio->index_reco_when_reio_start - 1; i++) {
    index_th=i+preio->rt_size;
    index_re=i+preio->index_reco_when_reio_start+1;
    pth->z_table[index_th]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_z];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_xe]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_xe];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_dkappa]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_dkappadtau];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_Tb]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_Tb];
    pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_cb2]=
      preco->recombination_table[index_re*preco->re_size+preco->index_re_cb2];
    // *** MY MODIFICATIONS ***
    if (pth->has_rayleigh_scattering==_TRUE_) {
      pth->thermodynamics_table[index_th*pth->th_size+pth->index_th_rayleigh_dkappa]=
        preco->recombination_table[index_re*preco->re_size+preco->index_re_rayleigh_dkappadtau];
    }
    // *** END OF MY MODIFICATIONS ***
  }

  // *** MY MODIFICATIONS ***
  /* During and after reionisation (z~10), the Rayleigh scattering is completely insignificant for
  two reasons: (1) it scales as (1+z)^6 and (2) the process of reionisation removes neutral hydrogen
  and helium, thus making Rayleigh scattering less efficient. Therefore, we just set Rayleigh 
  scattering to vanish after reionisation. */
  if (pth->has_rayleigh_scattering==_TRUE_) {
    for (i=1; i < preio->rt_size; i++) {
      pth->thermodynamics_table[i*pth->th_size+pth->index_th_rayleigh_dkappa] = 0;
    }
  }
  // *** END OF MY MODIFICATIONS ***

  /** - free the temporary structures */

  free(preco->recombination_table);

  if (pth->reio_parametrization != reio_none)
    free(preio->reionization_table);

  return _SUCCESS_;
}