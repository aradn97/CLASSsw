/** @file nonlinear.h Documented includes for trg module */

#include "primordial.h"
#include "extrapolate_source.h"
#include "trigonometric_integrals.h"

#ifndef __NONLINEAR__
#define __NONLINEAR__

#define _M_EV_TOO_BIG_FOR_HALOFIT_ 10. /**< above which value of non-CDM mass (in eV) do we stop trusting halofit? */

#define _M_SUN_ 1.98847e30 /**< Solar mass in Kg */

enum non_linear_method {nl_none,nl_halofit,nl_HMcode};
enum halofit_integral_type {halofit_integral_one, halofit_integral_two, halofit_integral_three};
enum hmcode_baryonic_feedback_model {nl_emu_dmonly, nl_owls_dmonly, nl_owls_ref, nl_owls_agn, nl_owls_dblim, nl_user_defined};

/**
 * Structure containing all information on non-linear spectra.
 *
 * Once initialized by nonlinear_init(), contains a table for all two points correlation functions
 * and for all the ai,bj functions (containing the three points correlation functions), for each
 * time and wave-number.
 */

struct nonlinear {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these
      parameters and the content of the 'precision', 'background',
      'thermo', 'primordial' and 'spectra' structures) */

  //@{

  enum non_linear_method method; /**< method for computing non-linear corrections (none, Halogit, etc.) */

  enum hmcode_baryonic_feedback_model feedback; /** to choose between different baryonic feedback models
                                                in hmcode (dmonly, gas cooling, Agn or supernova feedback) */
  double c_min;      /** for HMcode: minimum concentration in Bullock 2001 mass-concentration relation */
  double eta_0;      /** for HMcode: halo bloating parameter */
  double z_infinity; /** for HMcode: z value at which Dark Energy correction is evaluated needs to be at early times (default */

  //@}

  /** @name - information on number of modes and pairs of initial conditions */

  //@{

  int ic_size;         /**< for a given mode, ic_size[index_md] = number of initial conditions included in computation */
  int ic_ic_size;      /**< for a given mode, ic_ic_size[index_md] = number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just N(N+1)/2  where N = ic_size[index_md] */
  short * is_non_zero; /**< for a given mode, is_non_zero[index_md][index_ic1_ic2] is set to true if the pair of initial conditions (index_ic1, index_ic2) are statistically correlated, or to false if they are uncorrelated */

  //@}

  /** @name - table non-linear corrections for matter density, sqrt(P_NL(k,z)/P_NL(k,z)) */

  //@{

  short has_pk_m;  /**< do we want nonlinear corrections for total matter? */
  short has_pk_cb; /**< do we want nonlinear corrections for cdm+baryons? */

  int index_pk_m;  /**< index of pk for matter */
  int index_pk_cb; /**< index of pk for cold dark matter plus baryons */
  int pk_size;     /**< k_size = total number of pk */

  int k_size;      /**< k_size = total number of k values */
  double * k;      /**< k[index_k] = list of k values */
  double * ln_k;   /**< ln_k[index_k] = list of log(k) values */

  int k_size_extra;/** total number of k values of extrapolated k array (high k)*/
  double * k_extra;/** list of k-values with extrapolated high k-values  */

  int tau_size;    /**< tau_size = number of values */
  double * tau;    /**< tau[index_tau] = list of time values */

  double ** nl_corr_density;   /**< nl_corr_density[index_pk][index_tau * ppt->k_size + index_k] */
  double ** k_nl;              /**< wavenumber at which non-linear corrections become important,
                                    defined differently by different non_linear_method's */
  int index_tau_min_nl;        /**< index of smallest value of tau at which nonlinear corrections have been computed
                                    (so, for tau<tau_min_nl, the array nl_corr_density only contains some factors 1 */

  //@}

  /** @name - parameters for the pk_eq method */

  //@{

  short has_pk_eq;               /**< flag: will we use the pk_eq method? */

  int index_pk_eq_w;                /**< index of w in table pk_eq_w_and_Omega */
  int index_pk_eq_Omega_m;          /**< index of Omega_m in table pk_eq_w_and_Omega */
  int pk_eq_size;                   /**< number of indices in table pk_eq_w_and_Omega */

  int pk_eq_tau_size;               /**< number of times (and raws in table pk_eq_w_and_Omega) */

  double * pk_eq_tau;               /**< table of time values */
  double * pk_eq_w_and_Omega;       /**< table of background quantites */
  double * pk_eq_ddw_and_ddOmega;   /**< table of second derivatives */

  //@}

  /** @name - technical parameters */

  //@{

  short nonlinear_verbose;  	/**< amount of information written in standard output */

  ErrorMsg error_message; 	/**< zone for writing error messages */

  //@}
};

/**
 * Structure containing variables used only internally in nonlinear module by various functions.
 *
 */

struct nonlinear_workspace {

  /** @name - quantitites used by HMcode */

  //@{

  double * rtab; /** List of R values */
  double * stab; /** List of Sigma Values */
  double * ddstab; /** Splined sigma */

  double * growtable;
  double * ztable;
  double * tautable;

  double ** sigma_8;
  double ** sigma_disp;
  double ** sigma_disp_100;
  double ** sigma_prime;

  double dark_energy_correction; /** this is the ratio [g_wcdm(z_infinity)/g_lcdm(z_infinity)]^1.5
                                  * (power comes from Dolag et al. (2004) correction)
                                  * it is 1, if has_fld == _FALSE_ */

  //@}

};

/********************************************************************************/

/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int nonlinear_k_nl_at_z(
                          struct background *pba,
                          struct nonlinear * pnl,
                          double z,
                          double * k_nl,
                          double * k_nl_cb
                          );

  int nonlinear_init(
                     struct precision *ppr,
                     struct background *pba,
                     struct thermo *pth,
                     struct perturbs *ppt,
                     struct primordial *ppm,
                     struct nonlinear *pnl
                     );

  int nonlinear_free(
                     struct nonlinear *pnl
                     );

  // new functions:
  int nonlinear_pk_linear(
                          struct background *pba,
                          struct perturbs *ppt,
                          struct primordial *ppm,
                          struct nonlinear *pnl,
                          double ** sources,
                          int k_size,
                          int index_tau,
                          double *ln_pk_m_ic_l,
                          double *ln_pk_m_l,
                          double *ln_pk_cb_ic_l,
                          double *ln_pk_cb_l
                          );
  // end new functions

  int nonlinear_pk_l(struct background *pba,
                     struct perturbs *ppt,
                     struct primordial *ppm,
                     struct nonlinear *pnl,
                     int index_pk,
                     int index_tau,
                     double *pk_l,
                     double *lnk,
                     double *lnpk,
                     double *ddlnpk
                     );

  int nonlinear_halofit(
                        struct precision *ppr,
                        struct background *pba,
                        struct perturbs *ppt,
                        struct primordial *ppm,
                        struct nonlinear *pnl,
                        int index_pk,
                        double tau,
                        double *pk_l,
                        double *pk_nl,
                        double *lnk_l,
                        double *lnpk_l,
                        double *ddlnpk_l,
                        double *k_nl,
                        short * halofit_found_k_max
                        );

  int nonlinear_halofit_integrate(
                                  struct nonlinear *pnl,
                                  double * integrand_array,
                                  int integrand_size,
                                  int ia_size,
                                  int index_ia_k,
                                  int index_ia_pk,
                                  int index_ia_sum,
                                  int index_ia_ddsum,
                                  double R,
                                  enum halofit_integral_type type,
                                  double * sum
                                  );

  int nonlinear_hmcode(
                       struct precision *ppr,
                       struct background *pba,
                       struct perturbs *ppt,
                       struct primordial *ppm,
                       struct nonlinear *pnl,
                       int index_pk,
                       int index_tau,
                       double tau,
                       double *pk_l,
                       double *pk_nl,
                       double **lnk_l,
                       double **lnpk_l,
                       double **ddlnpk_l,
                       double *k_nl,
                       short * halofit_found_k_max,
                       struct nonlinear_workspace * pnw
                       );

  int nonlinear_hmcode_sigma(
                             struct precision * ppr,
                             struct background * pba,
                             struct perturbs * ppt,
                             struct primordial * ppm,
                             struct nonlinear * pnl,
                             double R,
                             double *lnk_l,
                             double *lnpk_l,
                             double *ddlnpk_l,
                             double * sigma
                             );


  int nonlinear_hmcode_sigma_prime(
                                   struct precision * ppr,
                                   struct background * pba,
                                   struct perturbs * ppt,
                                   struct primordial * ppm,
                                   struct nonlinear * pnl,
                                   double R,
                                   double *lnk_l,
                                   double *lnpk_l,
                                   double *ddlnpk_l,
                                   double * sigma_prime
                                   );

  int nonlinear_hmcode_sigma_disp(
                                  struct precision * ppr,
                                  struct background * pba,
                                  struct perturbs * ppt,
                                  struct primordial * ppm,
                                  struct nonlinear * pnl,
                                  double R,
                                  double *lnk_l,
                                  double *lnpk_l,
                                  double *ddlnpk_l,
                                  double * sigma_disp
                                  );

  int nonlinear_hmcode_fill_sigtab(
                                   struct precision *ppr,
                                   struct background * pba,
                                   struct perturbs *ppt,
                                   struct primordial * ppm,
                                   struct nonlinear * pnl,
                                   int index_tau,
                                   double *lnk_l,
                                   double *lnpk_l,
                                   double *ddlnpk_l,
                                   struct nonlinear_workspace * pnw
                                   );

  int nonlinear_hmcode_fill_growtab(
                                    struct precision *ppr,
                                    struct background * pba,
                                    struct nonlinear * pnl,
                                    struct nonlinear_workspace * pnw
                                    );

  int nonlinear_hmcode_halomassfunction(
                                        double nu,
                                        double *hmf
                                        );

  int nonlinear_hmcode_window_nfw(
                                  struct nonlinear * pnl,
                                  double k,
                                  double rv,
                                  double c,
                                  double *window_nfw
                                  );

  int nonlinear_hmcode_growint(
                               struct precision *ppr,
                               struct background * pba,
                               struct nonlinear * pnl,
                               double a,
                               double w,
                               double wa,
                               double * growth
                               );

  int nonlinear_hmcode_sigma8_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * sigma_8,
                        double * sigma_8_cb,
                        struct nonlinear_workspace * pnw
                        );

  int nonlinear_hmcode_sigmadisp_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * sigma_disp,
                        double * sigma_disp_cb,
                        struct nonlinear_workspace * pnw
                        );

  int nonlinear_hmcode_sigmadisp100_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * sigma_disp_100,
                        double * sigma_disp_100_cb,
                        struct nonlinear_workspace * pnw
                        );

  int nonlinear_hmcode_sigmaprime_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * sigma_prime,
                        double * sigma_prime_cb,
                        struct nonlinear_workspace * pnw
                        );

#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
