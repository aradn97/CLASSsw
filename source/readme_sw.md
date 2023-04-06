This file explains the changes I make to the CLASS code to accommodate swwerves, in chronological order.

********************************************************************
The terminal output without any changes:


Reading input parameters
 -> matched budget equations by adjusting Omega_Lambda = 0.690026
[WARNING: input line not used: 'Gamma_dcdm=0.0']
Running CLASS version v3.2.0
Computing background
 -> age = 13.770598 Gyr
 -> conformal age = 14151.897989 Mpc
 -> N_eff = 3.044 (summed over all species that are non-relativistic at early times) 
 -> radiation/matter equality at z = 3405.751104
    corresponding to conformal time = 112.722902 Mpc
Computing thermodynamics using HyRec 2020
 -> with primordial helium mass fraction Y_He = 0.2454
 -> recombination (maximum of visibility function) at z = 1088.775232
    corresponding to conformal time = 280.650005 Mpc
    with comoving sound horizon = 144.510134 Mpc
    angular diameter distance = 12.728540 Mpc
    sound horizon angle 100*theta_s = 1.041796
    Thomson optical depth crosses one at z_* = 1085.129318
    giving an angle 100*theta_* = 1.044144
 -> baryon drag stops at z = 1059.921811
    corresponding to conformal time = 286.316116 Mpc
    with comoving sound horizon rs = 147.071382 Mpc
 -> reionization with optical depth = 0.054308
    corresponding to conformal time = 5127.540389 Mpc
Computing sources
Computing primordial spectra (analytic spectrum)
Computing linear Fourier spectra.
 -> sigma8=0.824398 for total matter (computed till k = 1.08688 h/Mpc)
No non-linear spectra requested. Nonlinear calculations skipped.
Computing transfers
Computing unlensed harmonic spectra
Computing lensed spectra (fast mode)
Writing output files in output/default_sw15_... 
***********************************************************************

**First modfi.:   there is a decaying dark matter component in CLASS that we can turn on to represent swerving dm. So our dark matter is represented by dcdm in the code.

1. Change in background_functions():
 /* dcdm */
  if (pba->has_dcdm == _TRUE_) {
    /* Pass value of rho_dcdm to output */
    pvecback[pba->index_bg_rho_dcdm] = pvecback_B[pba->index_bi_rho_dcdm];
    rho_tot += pvecback[pba->index_bg_rho_dcdm];
    //p_tot += 0.; //swflg
    w_sdm = pvecback_B[pba->index_bg_w_dcdm]; //swflg
    pvecback[pba->index_bg_w_dcdm] = w_sdm; //swflg : this is the redundant part where B is also saved in A quantities
    p_tot += w_sdm * pvecback[pba->index_bg_rho_dcdm]; //swflg
    rho_m += pvecback[pba->index_bg_rho_dcdm];
  }
2. In background_indices():
class_define_index(pba->index_bg_w_dcdm,pba->has_dcdm,index_bg,1);  //swflg
class_define_index(pba->index_bi_w_dcdm,pba->has_dcdm,index_bi,1);  //swflg
3. In background_derivs():
dy[pba->index_bi_rho_dcdm] = -3.*y[pba->index_bi_rho_dcdm]*(1.+y[pba->index_bi_w_dcdm]) + pba->Gamma_dcdm/H*y[pba->index_bi_rho_dcdm]*(1.-23.*y[pba->index_bi_w_dcdm]/2.);  //swflg
dy[pba->index_bi_w_dcdm] = -2.*y[pba->index_bi_w_dcdm] + 3.*pow(y[pba->index_bi_w_dcdm],2) + 2./3.* pba->Gamma_dcdm/H - 26./3. * pba->Gamma_dcdm/H * y[pba->index_bi_w_dcdm] + 23/2 * pba->Gamma_dcdm/H*pow(y[pba->index_bi_w_dcdm],2);  //swflg
4. In int background_output_titles():
class_store_columntitle(titles,"(.)w_dcdm",pba->has_dcdm); //swflg
5. In background_output_data():
class_store_double(dataptr,pvecback[pba->index_bg_w_dcdm],pba->has_dcdm,storeidx); //swflg
6. In background_initial_conditions():
pvecback_integration[pba->index_bi_w_dcdm]  = 0.;
7. In background.h:
int index_bg_w_dcdm;  //swflg
int index_bi_w_dcdm;  /**< {B} dcdm equation of state*/ //swflg
8. In your input file
omega_dcdmdr = 0.1201075          # Reduced decaying dark matter density (Omega*h^2)
Gamma_dcdm = 1.0 
and comment out cdm.



















Edits to be done:
!  w_ini shouldn't be 0
!  in the file input.c, there is a piece of code like this:
  /* array of corresponding parameters that must be adjusted in order to meet the target (= unknown parameters) */
  char * const unknown_namestrings[] = {"h",                        /* unknown param for target '100*theta_s' */
                                        "Omega_ini_dcdm",           /* unknown param for target 'Omega_dcdmd' */
                                        "Omega_ini_dcdm",           /* unknown param for target 'omega_dcdmdr' */
                                        "scf_shooting_parameter",   /* unknown param for target 'Omega_scf' */
                                        "Omega_dcdmdr",             /* unknown param for target 'Omega_ini_dcdm' */
                                        "omega_dcdmdr"};             /* unknown param for target 'omega_ini_dcdm' */

   The second Omega_ini_dcdm shouldn't be omega_ini_dcdm?







