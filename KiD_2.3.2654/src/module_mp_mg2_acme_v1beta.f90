module module_mp_mg2_acme_v1beta

!---------------------------------------------------------------------------------
! Purpose:
!   MG microphysics version 2.0 - Update of MG microphysics with
!                                 prognostic precipitation.
!
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Peter Caldwell, Xiaohong Liu and Steve Ghan
! Version 2 history: Sep 2011: Development begun.
!                    Feb 2013: Added of prognostic precipitation.
! invoked in CAM by specifying -microphys=mg2.0
!
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------
!
! NOTE: Modified to allow other microphysics packages (e.g. CARMA) to do ice
! microphysics in cooperation with the MG liquid microphysics. This is
! controlled by the do_cldice variable.
!
! If do_cldice is false, then MG microphysics should not update CLDICE or
! NUMICE; it is assumed that the other microphysics scheme will have updated
! CLDICE and NUMICE. The other microphysics should handle the following
! processes that would have been done by MG:
!   - Detrainment (liquid and ice)
!   - Homogeneous ice nucleation
!   - Heterogeneous ice nucleation
!   - Bergeron process
!   - Melting of ice
!   - Freezing of cloud drops
!   - Autoconversion (ice -> snow)
!   - Growth/Sublimation of ice
!   - Sedimentation of ice
!
! This option has not been updated since the introduction of prognostic
! precipitation, and probably should be adjusted to cover snow as well.
!
!---------------------------------------------------------------------------------
! Based on micro_mg (restructuring of former cldwat2m_micro)
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------
! Code comments added by HM, 093011
! General code structure:
!
! Code is divided into two main subroutines:
!   subroutine micro_mg_init --> initializes microphysics routine, should be called
!                                  once at start of simulation
!   subroutine micro_mg_tend --> main microphysics routine to be called each time step
!                                this also calls several smaller subroutines to calculate
!                                microphysical processes and other utilities
!
! List of external functions:
!   qsat_water --> for calculating saturation vapor pressure with respect to liquid water
!   qsat_ice --> for calculating saturation vapor pressure with respect to ice
!   gamma   --> standard mathematical gamma function
! .........................................................................
! List of inputs through use statement in fortran90:
! Variable Name                      Description                Units
! .........................................................................
! gravit          acceleration due to gravity                    m s-2
! rair            dry air gas constant for air                  J kg-1 K-1
! tmelt           temperature of melting point for water          K
! cpair           specific heat at constant pressure for dry air J kg-1 K-1
! rh2o            gas constant for water vapor                  J kg-1 K-1
! latvap          latent heat of vaporization                   J kg-1
! latice          latent heat of fusion                         J kg-1
! qsat_water      external function for calculating liquid water
!                 saturation vapor pressure/humidity              -
! qsat_ice        external function for calculating ice
!                 saturation vapor pressure/humidity              pa
! rhmini          relative humidity threshold parameter for
!                 nucleating ice                                  -
! .........................................................................
! NOTE: List of all inputs/outputs passed through the call/subroutine statement
!       for micro_mg_tend is given below at the start of subroutine micro_mg_tend.
!---------------------------------------------------------------------------------

! Procedures required:
! 1) An implementation of the gamma function (if not intrinsic).
! 2) saturation vapor pressure and specific humidity over water
! 3) svp over ice


!#ifndef HAVE_GAMMA_INTRINSICS
!use shr_spfn_mod, only: gamma => shr_spfn_gamma
!#endif

use mphys_stats

use wv_sat_methods, only: &
     qsat_water => wv_sat_qsat_water, &
     qsat_ice => wv_sat_qsat_ice

! Parameters from the utilities module.
use micro_mg2_acme_v1beta_utils, only: &
     r8, &
     pi, &
     omsm, &
     qsmall, &
     mincld, &
     rhosn, &
     rhoi, &
     rhow, &
     rhows, &
     ac, bc, &
     ai, bi, &
     ar, br, &
     as, bs, &
     mi0, &
     rising_factorial

use namelists, only: KiD_outdir, KiD_outfile, fileNameOut

implicit none
private
save

public :: &
     micro_mg2_acme_v1beta_init, &
     micro_mg2_acme_v1beta_get_cols, &
     micro_mg2_acme_v1beta_tend

! switch for specification rather than prediction of droplet and crystal number
! note: number will be adjusted as needed to keep mean size within bounds,
! even when specified droplet or ice number is used

! If constant cloud ice number is set (nicons = .true.),
! then all microphysical processes except mass transfer due to ice nucleation
! (mnuccd) are based on the fixed cloud ice number. Calculation of
! mnuccd follows from the prognosed ice crystal number ni.

! nccons = .true. to specify constant cloud droplet number
! nicons = .true. to specify constant cloud ice number

logical, parameter, public :: nccons = .true. ! false originally
logical, parameter, public :: nicons = .true. ! false originally

!=========================================================
! Private module parameters
!=========================================================

! call counter
integer :: mphys_calls = 0

! file ID and file name for output
integer        :: fid_nstep
character(100) :: fname

! parameters for specified ice and droplet number concentration
! note: these are local in-cloud values, not grid-mean
real(r8), parameter :: ncnst = 100.e6_r8    ! droplet num concentration when nccons=.true. (m-3)
real(r8), parameter :: ninst = 0.1e6_r8     ! ice num concentration when nicons=.true. (m-3)

!Range of cloudsat reflectivities (dBz) for analytic simulator
real(r8), parameter :: csmin = -30._r8
real(r8), parameter :: csmax = 26._r8
real(r8), parameter :: mindbz = -99._r8
real(r8), parameter :: minrefl = 1.26e-10_r8    ! minrefl = 10._r8**(mindbz/10._r8)

! autoconversion size threshold for cloud ice to snow (m)
real(r8) :: dcs

!!== KZ_DCS
logical :: dcs_tdep
!!== KZ_DCS

! minimum mass of new crystal due to freezing of cloud droplets done
! externally (kg)
real(r8), parameter :: mi0l_min = 4._r8/3._r8*pi*rhow*(4.e-6_r8)**3

!=========================================================
! Constants set in initialization
!=========================================================

! Set using arguments to micro_mg_init
real(r8) :: g           ! gravity
real(r8) :: r           ! dry air gas constant
real(r8) :: rv          ! water vapor gas constant
real(r8) :: cpp         ! specific heat of dry air
real(r8) :: tmelt       ! freezing point of water (K)
real(r8) :: prc_coef1 = huge(1.0_r8)
real(r8) :: prc_exp   = huge(1.0_r8)
real(r8) :: prc_exp1  = huge(1.0_r8)
real(r8) :: cld_sed   = huge(1.0_r8) !scale fac for cld sedimentation velocity

logical  :: mg_prc_coeff_fix = .false.
! latent heats of:
real(r8) :: xxlv        ! vaporization
real(r8) :: xlf         ! freezing
real(r8) :: xxls        ! sublimation

real(r8) :: rhmini      ! Minimum rh for ice cloud fraction > 0.

! flags
logical :: microp_uniform
logical :: do_cldice
logical :: use_hetfrz_classnuc

real(r8) :: rhosu       ! typical 850mn air density

real(r8) :: icenuct     ! ice nucleation temperature: currently -5 degrees C

real(r8) :: snowmelt    ! what temp to melt all snow: currently 2 degrees C
real(r8) :: rainfrze    ! what temp to freeze all rain: currently -5 degrees C

! additional constants to help speed up code
real(r8) :: gamma_br_plus1
real(r8) :: gamma_br_plus4
real(r8) :: gamma_bs_plus1
real(r8) :: gamma_bs_plus4
real(r8) :: gamma_bi_plus1
real(r8) :: gamma_bi_plus4
real(r8) :: xxlv_squared
real(r8) :: xxls_squared

character(len=16)  :: micro_mg_precip_frac_method  ! type of precipitation fraction method
real(r8)           :: micro_mg_berg_eff_factor     ! berg efficiency factor

logical  :: allow_sed_supersat ! Allow supersaturated conditions after sedimentation loop

! Generic interface for packing routines
interface pack_array
   module procedure pack_array_1Dr8
   module procedure pack_array_2Dr8
   module procedure pack_array_3Dr8
end interface

interface unpack_array
   module procedure unpack_array_1Dr8
   module procedure unpack_array_1Dr8_arrayfill
   module procedure unpack_array_2Dr8
   module procedure unpack_array_2Dr8_arrayfill
end interface

!===============================================================================
contains
!===============================================================================

subroutine micro_mg2_acme_v1beta_init( &
     kind, gravit, rair, rh2o, cpair,    &
     tmelt_in, latvap, latice,           &
!!== KZ_DCS
     rhmini_in, micro_mg_dcs, micro_mg_dcs_tdep, &
!!== KZ_DCS
     microp_uniform_in, do_cldice_in, use_hetfrz_classnuc_in, &
     micro_mg_precip_frac_method_in, micro_mg_berg_eff_factor_in, &
     allow_sed_supersat_in, ice_sed_ai, prc_coef1_in,prc_exp_in,  &
     prc_exp1_in, cld_sed_in, mg_prc_coeff_fix_in, errstring)

  use micro_mg2_acme_v1beta_utils, only: micro_mg2_acme_v1beta_utils_init

  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! initialize constants for MG microphysics
  !
  ! Author: Andrew Gettelman Dec 2005
  !
  !-----------------------------------------------------------------------

  integer,  intent(in)  :: kind         ! Kind used for reals
  real(r8), intent(in)  :: gravit
  real(r8), intent(in)  :: rair
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: tmelt_in     ! Freezing point of water (K)
  real(r8), intent(in)  :: latvap
  real(r8), intent(in)  :: latice
  real(r8), intent(in)  :: rhmini_in    ! Minimum rh for ice cloud fraction > 0.
  real(r8), intent(in)  :: micro_mg_dcs
  real(r8), intent(in)  :: ice_sed_ai   !Fall speed parameter for cloud ice
!!== KZ_DCS
  logical,  intent(in)  :: micro_mg_dcs_tdep
!!== KZ_DCS

  logical,  intent(in)  :: microp_uniform_in    ! .true. = configure uniform for sub-columns
                                            ! .false. = use w/o sub-columns (standard)
  logical,  intent(in)  :: do_cldice_in     ! .true. = do all processes (standard)
                                            ! .false. = skip all processes affecting
                                            !           cloud ice
  logical,  intent(in)  :: use_hetfrz_classnuc_in ! use heterogeneous freezing

  character(len=16),intent(in)  :: micro_mg_precip_frac_method_in  ! type of precipitation fraction method
  real(r8),         intent(in)  :: micro_mg_berg_eff_factor_in     ! berg efficiency factor
  logical,  intent(in)  ::  allow_sed_supersat_in ! allow supersaturated conditions after sedimentation loop
  real(r8), intent(in)  :: prc_coef1_in,prc_exp_in,prc_exp1_in, cld_sed_in
  logical, intent(in)   :: mg_prc_coeff_fix_in


  character(128), intent(out) :: errstring    ! Output status (non-blank for error return)

  !-----------------------------------------------------------------------

  dcs = micro_mg_dcs

!!== KZ_DCS
  dcs_tdep = micro_mg_dcs_tdep
!!== KZ_DCS

 prc_coef1 = prc_coef1_in
 prc_exp   = prc_exp_in
 prc_exp1  = prc_exp1_in
 cld_sed   = cld_sed_in

 mg_prc_coeff_fix = mg_prc_coeff_fix_in

  ! Initialize subordinate utilities module.
  call micro_mg2_acme_v1beta_utils_init(kind, rh2o, cpair, tmelt_in, latvap, latice, &
       dcs, ice_sed_ai, errstring)

  if (trim(errstring) /= "") return

  ! declarations for MG code (transforms variable names)

  g= gravit                 ! gravity
  r= rair                   ! dry air gas constant: note units(phys_constants are in J/K/kmol)
  rv= rh2o                  ! water vapor gas constant
  cpp = cpair               ! specific heat of dry air
  tmelt = tmelt_in
  rhmini = rhmini_in
  micro_mg_precip_frac_method = micro_mg_precip_frac_method_in
  micro_mg_berg_eff_factor    = micro_mg_berg_eff_factor_in
  allow_sed_supersat          = allow_sed_supersat_in

  ! latent heats

  xxlv = latvap         ! latent heat vaporization
  xlf  = latice         ! latent heat freezing
  xxls = xxlv + xlf     ! latent heat of sublimation

  ! flags
  microp_uniform = microp_uniform_in
  do_cldice  = do_cldice_in
  use_hetfrz_classnuc = use_hetfrz_classnuc_in

  ! typical air density at 850 mb

  rhosu = 85000._r8/(rair * tmelt)

  ! Maximum temperature at which snow is allowed to exist
  snowmelt = tmelt + 2._r8
  ! Minimum temperature at which rain is allowed to exist
  rainfrze = tmelt - 40._r8

  ! Ice nucleation temperature
  icenuct  = tmelt - 5._r8

  ! Define constants to help speed up code (this limits calls to gamma function)
  gamma_br_plus1=gamma(1._r8+br)
  gamma_br_plus4=gamma(4._r8+br)
  gamma_bs_plus1=gamma(1._r8+bs)
  gamma_bs_plus4=gamma(4._r8+bs)
  gamma_bi_plus1=gamma(1._r8+bi)
  gamma_bi_plus4=gamma(4._r8+bi)
  xxlv_squared=xxlv**2
  xxls_squared=xxls**2

end subroutine micro_mg2_acme_v1beta_init

!===============================================================================
!microphysics routine for each timestep goes here...

subroutine micro_mg2_acme_v1beta_tend ( &
     mgncol, mgcols,               nlev, top_lev,                &
     deltatin,                                                   &
     tn,                           qn,                           &
     qcn,                          qin,                          &
     ncn,                          nin,                          &
     qrn,                          qsn,                          &
     nrn,                          nsn,                          &
     relvarn,                      accre_enhann,                 &
     pn,                           pdeln,                        &
     cldn,                         liqcldf,             icecldf, &
     qcsinksum_rate1ordo,                                        & ! output
     naain,                        npccnin,                      &
     rndstn,                       naconin,                      &
! outputs
     tlato,                        qvlato,                       &
     qctendo,                      qitendo,                      &
     nctendo,                      nitendo,                      &
     qrtendo,                      qstendo,                      &
     nrtendo,                      nstendo,                      &
     effco,                        effco_fn,            effio,   &
     precto,                       precio,                       &
     nevapro,                      evapsnowo,                    &
     praino,                       prodsnowo,                    &
     cmeouto,                      deffio,                       &
     pgamrado,                     lamcrado,                     &
     qsouto,                       dsouto,                       &
     rflxo,                        sflxo,               qrouto,  &
     reff_raino,                   reff_snowo,                   &
     qcsevapo,                     qisevapo,            qvreso,  &
     cmeiout,                      vtrmco,              vtrmio,  &
     umso,                         umro,                         &
     qcsedteno,                    qisedteno,                    &
     qrsedteno,                    qssedteno,                    &
     prao,                         prco,                         &
     mnuccco,                      mnuccto,             msacwio, &
     psacwso,                      bergso,              bergo,   &
     melto,                        homoo,                        &
     qcreso,                       prcio,               praio,   &
     qireso,                       mnuccro,             pracso,  &
     meltsdto,                     frzrdto,             mnuccdo, &
     nrouto,                       nsouto,                       &
     reflo,                        areflo,              areflzo, &
     freflo,                       csrflo,              acsrflo, &
     fcsrflo,                      rercldo,                      &
     ncaio,                        ncalo,                        &
     qrouto2,                      qsouto2,                      &
     nrouto2,                      nsouto2,                      &
     drouto2,                      dsouto2,                      &
     freqso,                       freqro,                       &
     nficeo,                       qcrato,                       &
! Below arguments are "optional" (pass null pointers to omit).
     tnd_qsnown,                   tnd_nsnown,          re_icen, &
     errstring,                                                  &
     prer_evapo,                                                 &
     frzimmn,                      frzcntn,             frzdepn)

  ! Constituent properties.
  use micro_mg2_acme_v1beta_utils, only: &
       mg_liq_props, &
       mg_ice_props, &
       mg_rain_props, &
       mg_snow_props

  ! Size calculation functions.
  use micro_mg2_acme_v1beta_utils, only: &
       size_dist_param_liq, &
!!== KZ_DCS
       size_dist_param_ice, &
!!== KZ_DCS
       size_dist_param_basic, &
       avg_diameter

  ! Microphysical processes.
  use micro_mg2_acme_v1beta_utils, only: &
       ice_deposition_sublimation, &
       kk2000_liq_autoconversion, &
       ice_autoconversion, &
       immersion_freezing, &
       contact_freezing, &
       snow_self_aggregation, &
       accrete_cloud_water_snow, &
       secondary_ice_production, &
       accrete_rain_snow, &
       heterogeneous_rain_freezing, &
       accrete_cloud_water_rain, &
       self_collection_rain, &
       accrete_cloud_ice_snow, &
       evaporate_sublimate_precip, &
       bergeron_process_snow

  use Sedimentation, only: &
       sed_CalcCFL, &
       sed_CalcFlux, &
       sed_AdvanceOneStep, &
       MG_LIQUID, &
       MG_ICE, &
       MG_RAIN,&
       MG_SNOW, &
       CFL

   use gptl, only: gptlstart_handle, gptlstop_handle
   use gptl_kid

  !Authors: Hugh Morrison, Andrew Gettelman, NCAR, Peter Caldwell, LLNL
  ! e-mail: morrison@ucar.edu, andrew@ucar.edu

  ! input arguments
  integer,  intent(in) :: mgncol         ! number of microphysics columns
  integer,  intent(in) :: mgcols(:)      ! list of microphysics columns
  integer,  intent(in) :: nlev           ! number of layers
  integer,  intent(in) :: top_lev        ! top level to do microphysics
  real(r8), intent(in) :: deltatin       ! time step (s)
  real(r8), intent(in) :: tn(:,:)        ! input temperature (K)
  real(r8), intent(in) :: qn(:,:)        ! input h20 vapor mixing ratio (kg/kg)

  ! note: all input cloud variables are grid-averaged
  real(r8), intent(in) :: qcn(:,:)       ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(:,:)       ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: ncn(:,:)       ! cloud water number conc (1/kg)
  real(r8), intent(in) :: nin(:,:)       ! cloud ice number conc (1/kg)

  real(r8), intent(in) :: qrn(:,:)       ! rain mixing ratio (kg/kg)
  real(r8), intent(in) :: qsn(:,:)       ! snow mixing ratio (kg/kg)
  real(r8), intent(in) :: nrn(:,:)       ! rain number conc (1/kg)
  real(r8), intent(in) :: nsn(:,:)       ! snow number conc (1/kg)

  real(r8), intent(in) :: relvarn(:,:)       ! cloud water relative variance (-)
  real(r8), intent(in) :: accre_enhann(:,:)  ! optional accretion
                                             ! enhancement factor (-)

  real(r8), intent(in) :: pn(:,:)        ! air pressure (pa)
  real(r8), intent(in) :: pdeln(:,:)     ! pressure difference across level (pa)

  real(r8), intent(in) :: cldn(:,:)      ! cloud fraction (no units)
  real(r8), intent(in) :: liqcldf(:,:)   ! liquid cloud fraction (no units)
  real(r8), intent(in) :: icecldf(:,:)   ! ice cloud fraction (no units)
  ! used for scavenging
  ! Inputs for aerosol activation
  real(r8), intent(in) :: naain(:,:)     ! ice nucleation number (from microp_aero_ts) (1/kg)
  real(r8), intent(in) :: npccnin(:,:)   ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)

  ! Note that for these variables, the dust bin is assumed to be the last index.
  ! (For example, in CAM, the last dimension is always size 4.)
  real(r8), intent(in) :: rndstn(:,:,:)  ! radius of each dust bin, for contact freezing (from microp_aero_ts) (m)
  real(r8), intent(in) :: naconin(:,:,:) ! number in each dust bin, for contact freezing  (from microp_aero_ts) (1/m^3)

  ! output arguments

  real(r8), intent(out) :: qcsinksum_rate1ordo(:,:)    ! 1st order rate for direct cw to precip conversion
  real(r8), intent(out) :: tlato(:,:)         ! latent heating rate       (W/kg)
  real(r8), intent(out) :: qvlato(:,:)        ! microphysical tendency qv (1/s)
  real(r8), intent(out) :: qctendo(:,:)       ! microphysical tendency qc (1/s)
  real(r8), intent(out) :: qitendo(:,:)       ! microphysical tendency qi (1/s)
  real(r8), intent(out) :: nctendo(:,:)       ! microphysical tendency nc (1/(kg*s))
  real(r8), intent(out) :: nitendo(:,:)       ! microphysical tendency ni (1/(kg*s))

  real(r8), intent(out) :: qrtendo(:,:)       ! microphysical tendency qr (1/s)
  real(r8), intent(out) :: qstendo(:,:)       ! microphysical tendency qs (1/s)
  real(r8), intent(out) :: nrtendo(:,:)       ! microphysical tendency nr (1/(kg*s))
  real(r8), intent(out) :: nstendo(:,:)       ! microphysical tendency ns (1/(kg*s))

  real(r8), intent(out) :: effco(:,:)         ! droplet effective radius (micron)
  real(r8), intent(out) :: effco_fn(:,:)      ! droplet effective radius, assuming nc = 1.e8 kg-1
  real(r8), intent(out) :: effio(:,:)         ! cloud ice effective radius (micron)
  real(r8), intent(out) :: precto(:)          ! surface precip rate (m/s)
  real(r8), intent(out) :: precio(:)          ! cloud ice/snow precip rate (m/s)
  real(r8), intent(out) :: nevapro(:,:)       ! evaporation rate of rain + snow (1/s)
  real(r8), intent(out) :: evapsnowo(:,:)     ! sublimation rate of snow (1/s)
  real(r8), intent(out) :: praino(:,:)        ! production of rain + snow (1/s)
  real(r8), intent(out) :: prodsnowo(:,:)     ! production of snow (1/s)
  real(r8), intent(out) :: cmeouto(:,:)       ! evap/sub of cloud (1/s)
  real(r8), intent(out) :: deffio(:,:)        ! ice effective diameter for optics (radiation) (micron)
  real(r8), intent(out) :: pgamrado(:,:)      ! ice gamma parameter for optics (radiation) (no units)
  real(r8), intent(out) :: lamcrado(:,:)      ! slope of droplet distribution for optics (radiation) (1/m)
  real(r8), intent(out) :: qsouto(:,:)        ! snow mixing ratio (kg/kg)
  real(r8), intent(out) :: dsouto(:,:)        ! snow diameter (m)
  real(r8), intent(out) :: rflxo(:,:)         ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8), intent(out) :: sflxo(:,:)         ! grid-box average snow flux (kg m^-2 s^-1)
  real(r8), intent(out) :: qrouto(:,:)        ! grid-box average rain mixing ratio (kg/kg)
  real(r8), intent(out) :: reff_raino(:,:)    ! rain effective radius (micron)
  real(r8), intent(out) :: reff_snowo(:,:)    ! snow effective radius (micron)
  real(r8), intent(out) :: qcsevapo(:,:)      ! cloud water evaporation due to sedimentation (1/s)
  real(r8), intent(out) :: qisevapo(:,:)      ! cloud ice sublimation due to sublimation (1/s)
  real(r8), intent(out) :: qvreso(:,:)        ! residual condensation term to ensure RH < 100% (1/s)
  real(r8), intent(out) :: cmeiout(:,:)       ! grid-mean cloud ice sub/dep (1/s)
  real(r8), intent(out) :: vtrmco(:,:)        ! mass-weighted cloud water fallspeed (m/s)
  real(r8), intent(out) :: vtrmio(:,:)        ! mass-weighted cloud ice fallspeed (m/s)
  real(r8), intent(out) :: umso(:,:)          ! mass weighted snow fallspeed (m/s)
  real(r8), intent(out) :: umro(:,:)          ! mass weighted rain fallspeed (m/s)
  real(r8), intent(out) :: qcsedteno(:,:)     ! qc sedimentation tendency (1/s)
  real(r8), intent(out) :: qisedteno(:,:)     ! qi sedimentation tendency (1/s)
  real(r8), intent(out) :: qrsedteno(:,:)     ! qr sedimentation tendency (1/s)
  real(r8), intent(out) :: qssedteno(:,:)     ! qs sedimentation tendency (1/s)

  ! microphysical process rates for output (mixing ratio tendencies) (all have units of 1/s)
  real(r8), intent(out) :: prao(:,:)          ! accretion of cloud by rain
  real(r8), intent(out) :: prco(:,:)          ! autoconversion of cloud to rain
  real(r8), intent(out) :: mnuccco(:,:)       ! mixing ratio tend due to immersion freezing
  real(r8), intent(out) :: mnuccto(:,:)       ! mixing ratio tend due to contact freezing
  real(r8), intent(out) :: msacwio(:,:)       ! mixing ratio tend due to H-M splintering
  real(r8), intent(out) :: psacwso(:,:)       ! collection of cloud water by snow
  real(r8), intent(out) :: bergso(:,:)        ! bergeron process on snow
  real(r8), intent(out) :: bergo(:,:)         ! bergeron process on cloud ice
  real(r8), intent(out) :: melto(:,:)         ! melting of cloud ice
  real(r8), intent(out) :: homoo(:,:)         ! homogeneous freezing cloud water
  real(r8), intent(out) :: qcreso(:,:)        ! residual cloud condensation due to removal of excess supersat
  real(r8), intent(out) :: prcio(:,:)         ! autoconversion of cloud ice to snow
  real(r8), intent(out) :: praio(:,:)         ! accretion of cloud ice by snow
  real(r8), intent(out) :: qireso(:,:)        ! residual ice deposition due to removal of excess supersat
  real(r8), intent(out) :: mnuccro(:,:)       ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
  real(r8), intent(out) :: pracso(:,:)        ! mixing ratio tendency due to accretion of rain by snow (1/s)
  real(r8), intent(out) :: meltsdto(:,:)      ! latent heating rate due to melting of snow  (W/kg)
  real(r8), intent(out) :: frzrdto(:,:)       ! latent heating rate due to homogeneous freezing of rain (W/kg)
  real(r8), intent(out) :: mnuccdo(:,:)       ! mass tendency from ice nucleation
  real(r8), intent(out) :: nrouto(:,:)        ! rain number concentration (1/m3)
  real(r8), intent(out) :: nsouto(:,:)        ! snow number concentration (1/m3)
  real(r8), intent(out) :: reflo(:,:)         ! analytic radar reflectivity
  real(r8), intent(out) :: areflo(:,:)        ! average reflectivity will zero points outside valid range
  real(r8), intent(out) :: areflzo(:,:)       ! average reflectivity in z.
  real(r8), intent(out) :: freflo(:,:)        ! fractional occurrence of radar reflectivity
  real(r8), intent(out) :: csrflo(:,:)        ! cloudsat reflectivity
  real(r8), intent(out) :: acsrflo(:,:)       ! cloudsat average
  real(r8), intent(out) :: fcsrflo(:,:)       ! cloudsat fractional occurrence of radar reflectivity
  real(r8), intent(out) :: rercldo(:,:)       ! effective radius calculation for rain + cloud
  real(r8), intent(out) :: ncaio(:,:)         ! output number conc of ice nuclei available (1/m3)
  real(r8), intent(out) :: ncalo(:,:)         ! output number conc of CCN (1/m3)
  real(r8), intent(out) :: qrouto2(:,:)       ! copy of qrout as used to compute drout2
  real(r8), intent(out) :: qsouto2(:,:)       ! copy of qsout as used to compute dsout2
  real(r8), intent(out) :: nrouto2(:,:)       ! copy of nrout as used to compute drout2
  real(r8), intent(out) :: nsouto2(:,:)       ! copy of nsout as used to compute dsout2
  real(r8), intent(out) :: drouto2(:,:)       ! mean rain particle diameter (m)
  real(r8), intent(out) :: dsouto2(:,:)       ! mean snow particle diameter (m)
  real(r8), intent(out) :: freqso(:,:)        ! fractional occurrence of snow
  real(r8), intent(out) :: freqro(:,:)        ! fractional occurrence of rain
  real(r8), intent(out) :: nficeo(:,:)        ! fractional occurrence of ice
  real(r8), intent(out) :: qcrato(:,:)        ! limiter for qc process rates (1=no limit --> 0. no qc)

  real(r8), intent(out) :: prer_evapo(:,:)

  character(128),   intent(out) :: errstring  ! output status (non-blank for error return)

  ! Tendencies calculated by external schemes that can replace MG's native
  ! process tendencies.

  ! Used with CARMA cirrus microphysics
  ! (or similar external microphysics model)
  real(r8), intent(in) :: tnd_qsnown(:,:) ! snow mass tendency (kg/kg/s)
  real(r8), intent(in) :: tnd_nsnown(:,:) ! snow number tendency (#/kg/s)
  real(r8), intent(in) :: re_icen(:,:)    ! ice effective radius (m)

  ! From external ice nucleation.
  real(r8), intent(in) :: frzimmn(:,:) ! Number tendency due to immersion freezing (1/cm3)
  real(r8), intent(in) :: frzcntn(:,:) ! Number tendency due to contact freezing (1/cm3)
  real(r8), intent(in) :: frzdepn(:,:) ! Number tendency due to deposition nucleation (1/cm3)

  ! local workspace
  ! all units mks unless otherwise stated

  ! local copies of input variables
  real(r8) :: q(mgncol,nlev)       ! water vapor mixing ratio (kg/kg)
  real(r8) :: t(mgncol,nlev)       ! temperature (K)
  real(r8) :: qc(mgncol,nlev)      ! cloud liquid mixing ratio (kg/kg)
  real(r8) :: qi(mgncol,nlev)      ! cloud ice mixing ratio (kg/kg)
  real(r8) :: nc(mgncol,nlev)      ! cloud liquid number concentration (1/kg)
  real(r8) :: ni(mgncol,nlev)      ! cloud liquid number concentration (1/kg)
  real(r8) :: qr(mgncol,nlev)      ! rain mixing ratio (kg/kg)
  real(r8) :: qs(mgncol,nlev)      ! snow mixing ratio (kg/kg)
  real(r8) :: nr(mgncol,nlev)      ! rain number concentration (1/kg)
  real(r8) :: ns(mgncol,nlev)      ! snow number concentration (1/kg)

  real(r8) :: p(mgncol,nlev)       ! pressure (Pa)
  real(r8) :: pdel(mgncol,nlev)    ! pressure difference across level (Pa)

  real(r8) :: naai(mgncol,nlev)    ! ice nucleation number (from microp_aero_ts) (1/kg)
  real(r8) :: npccn(mgncol,nlev)   ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)

  real(r8) :: relvar(mgncol,nlev)      ! cloud water relative variance (-)
  real(r8) :: accre_enhan(mgncol,nlev) ! optional accretion enhancement factor (-)

  ! These were made allocatable because of a problem on PGI (possibly
  ! exceeding the stack size limit on some machines).
  real(r8), allocatable :: rndst(:,:,:)
  real(r8), allocatable :: nacon(:,:,:)

  real(r8) :: tnd_qsnow(mgncol,nlev) ! snow mass tendency (kg/kg/s)
  real(r8) :: tnd_nsnow(mgncol,nlev) ! snow number tendency (#/kg/s)
  real(r8) :: re_ice(mgncol,nlev)    ! ice effective radius (m)

  real(r8) :: frzimm(mgncol,nlev) ! Number tendency due to immersion freezing (1/cm3)
  real(r8) :: frzcnt(mgncol,nlev) ! Number tendency due to contact freezing (1/cm3)
  real(r8) :: frzdep(mgncol,nlev) ! Number tendency due to deposition nucleation (1/cm3)


  ! Packed copies of output variables
  real(r8) :: qcsinksum_rate1ord(mgncol,nlev) ! 1st order rate for direct cw to precip conversion
  real(r8) :: tlat(mgncol,nlev)      ! latent heating rate       (W/kg)
  real(r8) :: qvlat(mgncol,nlev)     ! microphysical tendency qv (1/s)
  real(r8) :: qctend(mgncol,nlev)    ! microphysical tendency qc (1/s)
  real(r8) :: qitend(mgncol,nlev)    ! microphysical tendency qi (1/s)
  real(r8) :: nctend(mgncol,nlev)    ! microphysical tendency nc (1/(kg*s))
  real(r8) :: nitend(mgncol,nlev)    ! microphysical tendency ni (1/(kg*s))

  real(r8) :: effc(mgncol,nlev)      ! droplet effective radius (micron)
  real(r8) :: effc_fn(mgncol,nlev)   ! droplet effective radius, assuming nc = 1.e8 kg-1
  real(r8) :: effi(mgncol,nlev)      ! cloud ice effective radius (micron)

  real(r8) :: prect(mgncol)          ! surface precip rate (m/s)
  real(r8) :: preci(mgncol)          ! cloud ice/snow precip rate (m/s)

  real(r8) :: nevapr(mgncol,nlev)    ! evaporation rate of rain + snow (1/s)
  real(r8) :: evapsnow(mgncol,nlev)  ! sublimation rate of snow (1/s)
  real(r8) :: prain(mgncol,nlev)     ! production of rain + snow (1/s)
  real(r8) :: prodsnow(mgncol,nlev)  ! production of snow (1/s)
  real(r8) :: cmeout(mgncol,nlev)    ! evap/sub of cloud (1/s)
  real(r8) :: deffi(mgncol,nlev)     ! ice effective diameter for optics (radiation) (micron)
  real(r8) :: pgamrad(mgncol,nlev)   ! ice gamma parameter for optics (radiation) (no units)
  real(r8) :: lamcrad(mgncol,nlev)   ! slope of droplet distribution for optics (radiation) (1/m)


  real(r8) :: qsout(mgncol,nlev)     ! snow mixing ratio (kg/kg)
  real(r8) :: qsout2(mgncol,nlev)    ! copy of qsout as used to compute dsout2
  real(r8) :: nsout(mgncol,nlev)     ! snow number concentration (1/m3)
  real(r8) :: nsout2(mgncol,nlev)    ! copy of nsout as used to compute dsout2
  real(r8) :: dsout(mgncol,nlev)     ! snow diameter (m)
  real(r8) :: dsout2(mgncol,nlev)    ! mean snow particle diameter (m)

  real(r8) :: qrout(mgncol,nlev)     ! grid-box average rain mixing ratio (kg/kg)
  real(r8) :: qrout2(mgncol,nlev)    ! copy of qrout as used to compute drout2
  real(r8) :: nrout(mgncol,nlev)     ! rain number concentration (1/m3)
  real(r8) :: nrout2(mgncol,nlev)    ! copy of nrout as used to compute drout2
  real(r8) :: drout2(mgncol,nlev)    ! mean rain particle diameter (m)

  real(r8) :: reff_rain(mgncol,nlev) ! rain effective radius (micron)
  real(r8) :: reff_snow(mgncol,nlev) ! snow effective radius (micron)

  real(r8) :: freqs(mgncol,nlev)     ! fractional occurrence of snow
  real(r8) :: freqr(mgncol,nlev)     ! fractional occurrence of rain

  real(r8) :: rflx(mgncol,nlev+1)    ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8) :: sflx(mgncol,nlev+1)    ! grid-box average snow flux (kg m^-2 s^-1)

  real(r8) :: qcsevap(mgncol,nlev)   ! cloud water evaporation due to sedimentation (1/s)
  real(r8) :: qisevap(mgncol,nlev)   ! cloud ice sublimation due to sublimation (1/s)
  real(r8) :: qvres(mgncol,nlev)     ! residual condensation term to ensure RH < 100% (1/s)
  real(r8) :: cmeitot(mgncol,nlev)   ! grid-mean cloud ice sub/dep (1/s)
  real(r8) :: vtrmc(mgncol,nlev)     ! mass-weighted cloud water fallspeed (m/s)
  real(r8) :: vtrmi(mgncol,nlev)     ! mass-weighted cloud ice fallspeed (m/s)
  real(r8) :: qcsedten(mgncol,nlev)  ! qc sedimentation tendency (1/s)
  real(r8) :: qisedten(mgncol,nlev)  ! qi sedimentation tendency (1/s)
  real(r8) :: qrsedten(mgncol,nlev)  ! qr sedimentation tendency (1/s)
  real(r8) :: qssedten(mgncol,nlev)  ! qs sedimentation tendency (1/s)

  real(r8) :: pratot(mgncol,nlev)    ! accretion of cloud by rain
  real(r8) :: prctot(mgncol,nlev)    ! autoconversion of cloud to rain
  real(r8) :: mnuccctot(mgncol,nlev) ! mixing ratio tend due to immersion freezing
  real(r8) :: mnuccttot(mgncol,nlev) ! mixing ratio tend due to contact freezing
  real(r8) :: msacwitot(mgncol,nlev) ! mixing ratio tend due to H-M splintering
  real(r8) :: psacwstot(mgncol,nlev) ! collection of cloud water by snow
  real(r8) :: bergstot(mgncol,nlev)  ! bergeron process on snow
  real(r8) :: bergtot(mgncol,nlev)   ! bergeron process on cloud ice
  real(r8) :: melttot(mgncol,nlev)   ! melting of cloud ice
  real(r8) :: homotot(mgncol,nlev)   ! homogeneous freezing cloud water
  real(r8) :: qcrestot(mgncol,nlev)  ! residual cloud condensation due to removal of excess supersat
  real(r8) :: prcitot(mgncol,nlev)   ! autoconversion of cloud ice to snow
  real(r8) :: praitot(mgncol,nlev)   ! accretion of cloud ice by snow
  real(r8) :: qirestot(mgncol,nlev)  ! residual ice deposition due to removal of excess supersat
  real(r8) :: mnuccrtot(mgncol,nlev) ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
  real(r8) :: pracstot(mgncol,nlev)  ! mixing ratio tendency due to accretion of rain by snow (1/s)
  real(r8) :: mnuccdtot(mgncol,nlev) ! mass tendency from ice nucleation
  real(r8) :: meltsdttot(mgncol,nlev)! latent heating rate due to melting of snow  (W/kg)
  real(r8) :: frzrdttot(mgncol,nlev) ! latent heating rate due to homogeneous freezing of rain (W/kg)

  real(r8) :: refl(mgncol,nlev)      ! analytic radar reflectivity
  real(r8) :: arefl(mgncol,nlev)     ! average reflectivity will zero points outside valid range
  real(r8) :: areflz(mgncol,nlev)    ! average reflectivity in z.
  real(r8) :: frefl(mgncol,nlev)     ! fractional occurrence of radar reflectivity
  real(r8) :: csrfl(mgncol,nlev)     ! cloudsat reflectivity
  real(r8) :: acsrfl(mgncol,nlev)    ! cloudsat average
  real(r8) :: fcsrfl(mgncol,nlev)    ! cloudsat fractional occurrence of radar reflectivity

  real(r8) :: rercld(mgncol,nlev)    ! effective radius calculation for rain + cloud

  real(r8) :: nfice(mgncol,nlev)     ! fractional occurrence of ice

  real(r8) :: ncai(mgncol,nlev)      ! output number conc of ice nuclei available (1/m3)
  real(r8) :: ncal(mgncol,nlev)      ! output number conc of CCN (1/m3)

  ! general purpose variables
  real(r8) :: deltat            ! sub-time step (s)
  real(r8) :: mtime             ! the assumed ice nucleation timescale

  ! physical properties of the air at a given point
  real(r8) :: rho(mgncol,nlev)    ! density (kg m-3)
  real(r8) :: dv(mgncol,nlev)     ! diffusivity of water vapor
  real(r8) :: mu(mgncol,nlev)     ! viscosity
  real(r8) :: sc(mgncol,nlev)     ! schmidt number
  real(r8) :: rhof(mgncol,nlev)   ! density correction factor for fallspeed

  ! cloud fractions
  real(r8) :: precip_frac(mgncol,nlev) ! precip fraction assuming maximum overlap
  real(r8) :: cldm(mgncol,nlev)   ! cloud fraction
  real(r8) :: icldm(mgncol,nlev)  ! ice cloud fraction
  real(r8) :: lcldm(mgncol,nlev)  ! liq cloud fraction

  ! mass mixing ratios
  real(r8) :: qcic(mgncol,nlev)   ! in-cloud cloud liquid
  real(r8) :: qiic(mgncol,nlev)   ! in-cloud cloud ice
  real(r8) :: qsic(mgncol,nlev)   ! in-precip snow
  real(r8) :: qric(mgncol,nlev)   ! in-precip rain

  ! number concentrations
  real(r8) :: ncic(mgncol,nlev)   ! in-cloud droplet
  real(r8) :: niic(mgncol,nlev)   ! in-cloud cloud ice
  real(r8) :: nsic(mgncol,nlev)   ! in-precip snow
  real(r8) :: nric(mgncol,nlev)   ! in-precip rain
  ! maximum allowed ni value
  real(r8) :: nimax(mgncol,nlev)

  ! Size distribution parameters for:
  ! cloud ice
  real(r8) :: lami(mgncol,nlev)   ! slope
  real(r8) :: n0i(mgncol,nlev)    ! intercept
  ! cloud liquid
  real(r8) :: lamc(mgncol,nlev)   ! slope
  real(r8) :: pgam(mgncol,nlev)   ! spectral width parameter
  ! snow
  real(r8) :: lams(mgncol,nlev)   ! slope
  real(r8) :: n0s(mgncol,nlev)    ! intercept
  ! rain
  real(r8) :: lamr(mgncol,nlev)   ! slope
  real(r8) :: n0r(mgncol,nlev)    ! intercept

  ! Rates/tendencies due to:

  ! Instantaneous snow melting
  real(r8) :: minstsm(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: ninstsm(mgncol,nlev)    ! number concentration
  ! Instantaneous rain freezing
  real(r8) :: minstrf(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: ninstrf(mgncol,nlev)    ! number concentration

  ! deposition of cloud ice
  real(r8) :: vap_dep(mgncol,nlev)    ! deposition from vapor to ice PMC 12/3/12
  ! sublimation of cloud ice
  real(r8) :: ice_sublim(mgncol,nlev) ! sublimation from ice to vapor PMC 12/3/12
  ! ice nucleation
  real(r8) :: nnuccd(mgncol,nlev) ! number rate from deposition/cond.-freezing
  real(r8) :: mnuccd(mgncol,nlev) ! mass mixing ratio
  ! freezing of cloud water
  real(r8) :: mnuccc(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnuccc(mgncol,nlev) ! number concentration
  ! contact freezing of cloud water
  real(r8) :: mnucct(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnucct(mgncol,nlev) ! number concentration
  ! deposition nucleation in mixed-phase clouds (from external scheme)
  real(r8) :: mnudep(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnudep(mgncol,nlev) ! number concentration
  ! ice multiplication
  real(r8) :: msacwi(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nsacwi(mgncol,nlev) ! number concentration
  ! autoconversion of cloud droplets
  real(r8) :: prc(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: nprc(mgncol,nlev)   ! number concentration (rain)
  real(r8) :: nprc1(mgncol,nlev)  ! number concentration (cloud droplets)
  ! self-aggregation of snow
  real(r8) :: nsagg(mgncol,nlev)  ! number concentration
  ! self-collection of rain
  real(r8) :: nragg(mgncol,nlev)  ! number concentration
  ! collection of droplets by snow
  real(r8) :: psacws(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: npsacws(mgncol,nlev)    ! number concentration
  ! collection of rain by snow
  real(r8) :: pracs(mgncol,nlev)  ! mass mixing ratio
  real(r8) :: npracs(mgncol,nlev) ! number concentration
  ! freezing of rain
  real(r8) :: mnuccr(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnuccr(mgncol,nlev) ! number concentration
  ! freezing of rain to form ice (mg add 4/26/13)
  real(r8) :: mnuccri(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: nnuccri(mgncol,nlev)    ! number concentration
  ! accretion of droplets by rain
  real(r8) :: pra(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: npra(mgncol,nlev)   ! number concentration
  ! autoconversion of cloud ice to snow
  real(r8) :: prci(mgncol,nlev)   ! mass mixing ratio
  real(r8) :: nprci(mgncol,nlev)  ! number concentration
  ! accretion of cloud ice by snow
  real(r8) :: prai(mgncol,nlev)   ! mass mixing ratio
  real(r8) :: nprai(mgncol,nlev)  ! number concentration
  ! evaporation of rain
  real(r8) :: pre(mgncol,nlev)    ! mass mixing ratio
  ! sublimation of snow
  real(r8) :: prds(mgncol,nlev)   ! mass mixing ratio
  ! number evaporation
  real(r8) :: nsubi(mgncol,nlev)  ! cloud ice
  real(r8) :: nsubc(mgncol,nlev)  ! droplet
  real(r8) :: nsubs(mgncol,nlev)  ! snow
  real(r8) :: nsubr(mgncol,nlev)  ! rain
  ! bergeron process
  real(r8) :: berg(mgncol,nlev)   ! mass mixing ratio (cloud ice)
  real(r8) :: bergs(mgncol,nlev)  ! mass mixing ratio (snow)

  real(r8) :: qcrat(mgncol,nlev)  ! ratio for qc limiter

  real(r8) :: prer_evap(mgncol,nlev)

  ! fallspeeds
  ! number-weighted
  real(r8) :: uns(mgncol,nlev)    ! snow
  real(r8) :: unr(mgncol,nlev)    ! rain
  ! mass-weighted
  real(r8) :: ums(mgncol,nlev)    ! snow
  real(r8) :: umr(mgncol,nlev)    ! rain
  ! air density corrected fallspeed parameters
  real(r8) :: arn(mgncol,nlev)    ! rain
  real(r8) :: asn(mgncol,nlev)    ! snow
  real(r8) :: acn(mgncol,nlev)    ! cloud droplet
  real(r8) :: ain(mgncol,nlev)    ! cloud ice

  ! Mass of liquid droplets used with external heterogeneous freezing.
  real(r8) :: mi0l(mgncol)

  ! saturation vapor pressures
  real(r8) :: esl(mgncol,nlev)    ! liquid
  real(r8) :: esi(mgncol,nlev)    ! ice
  real(r8) :: esn                 ! checking for RH after rain evap

  ! saturation vapor mixing ratios
  real(r8) :: qvl(mgncol,nlev)    ! liquid
  real(r8) :: qvi(mgncol,nlev)    ! ice
  real(r8) :: qvn                 ! checking for RH after rain evap

  ! relative humidity
  real(r8) :: relhum(mgncol,nlev)

!!== KZ_DCS
  real(r8) :: dcst(mgncol,nlev)        ! t-dependent dcs
!!== KZ_DCS

  ! parameters for cloud water and cloud ice sedimentation calculations
  real(r8) :: fc(0:nlev,2), fnc(0:nlev,2)
  real(r8) :: fi(0:nlev,2), fni(0:nlev,2)
  real(r8) :: fr(0:nlev,2), fnr(0:nlev,2)
  real(r8) :: fs(0:nlev,2), fns(0:nlev,2)
  real(r8) :: alphaq(0:nlev), alphan(0:nlev)
  real(r8) :: s1(nlev),s2(nlev), w1(2,nlev), w2(2,nlev)

  ! sum of source/sink terms for diagnostic precip
  real(r8) :: qstend(mgncol,nlev)     ! snow mixing ratio
  real(r8) :: nstend(mgncol,nlev)     ! snow number concentration
  real(r8) :: qrtend(mgncol,nlev)     ! rain mixing ratio
  real(r8) :: nrtend(mgncol,nlev)     ! rain number concentration

  real(r8) :: rainrt(mgncol,nlev)     ! rain rate for reflectivity calculation

  ! dummy variables
  real(r8) :: dum
  real(r8) :: dum1
  real(r8) :: dum2
  ! dummies for checking RH
  real(r8) :: qtmp
  real(r8) :: ttmp
  ! dummies for conservation check
  real(r8) :: ratio
  real(r8) :: tmpfrz
  ! dummies for in-cloud variables
  real(r8) :: dumc(mgncol,nlev)   ! qc
  real(r8) :: dumnc(mgncol,nlev)  ! nc
  real(r8) :: dumi(mgncol,nlev)   ! qi
  real(r8) :: dumni(mgncol,nlev)  ! ni
  real(r8) :: dumr(mgncol,nlev)   ! rain mixing ratio
  real(r8) :: dumnr(mgncol,nlev)  ! rain number concentration
  real(r8) :: dums(mgncol,nlev)   ! snow mixing ratio
  real(r8) :: dumns(mgncol,nlev)  ! snow number concentration
  ! Array dummy variable
  real(r8) :: dum_2D(mgncol,nlev)

  ! sedimentation substep loop variables
  real(r8) :: deltat_sed, time_sed
  ! number of sub-steps for loops over "n" (for sedimentation)
  integer nstep

  ! variable for timing loops
  integer :: n


  ! loop array variables
  ! "i" and "k" are column/level iterators for internal (MG) variables
  integer i, k


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! default return error message
  errstring = ' '

  ! Process inputs

  ! assign variable deltat to deltatin
  deltat = deltatin

  call pack_array(     qn, mgcols, top_lev,     q)
  call pack_array(     tn, mgcols, top_lev,     t)
  call pack_array(    qcn, mgcols, top_lev,    qc)
  call pack_array(    qin, mgcols, top_lev,    qi)
  call pack_array(    ncn, mgcols, top_lev,    nc)
  call pack_array(    nin, mgcols, top_lev,    ni)
  call pack_array(    qrn, mgcols, top_lev,    qr)
  call pack_array(    qsn, mgcols, top_lev,    qs)
  call pack_array(    nrn, mgcols, top_lev,    nr)
  call pack_array(    nsn, mgcols, top_lev,    ns)

  ! Relative variance used only if microp_uniform = .false.
  call pack_array(relvarn, mgcols, top_lev,relvar)
  call pack_array(accre_enhann, mgcols, top_lev, accre_enhan)

  call pack_array(     pn, mgcols, top_lev,     p)
  call pack_array(  pdeln, mgcols, top_lev,  pdel)

  call pack_array(  naain, mgcols, top_lev,  naai)
  call pack_array(npccnin, mgcols, top_lev, npccn)

  ! These are allocated instead of used as automatic arrays
  ! purely to work around a PGI bug.
  allocate(rndst(mgncol,nlev,size(rndstn,3)))
  allocate(nacon(mgncol,nlev,size(rndstn,3)))
  call pack_array( rndstn, mgcols, top_lev, rndst)
  call pack_array(naconin, mgcols, top_lev, nacon)

  if (.not. do_cldice) then
     call pack_array(tnd_qsnown, mgcols, top_lev, tnd_qsnow)
     call pack_array(tnd_nsnown, mgcols, top_lev, tnd_nsnow)
     call pack_array(   re_icen, mgcols, top_lev,    re_ice)
  end if

  if (use_hetfrz_classnuc) then
     call pack_array(frzimmn, mgcols, top_lev, frzimm)
     call pack_array(frzcntn, mgcols, top_lev, frzcnt)
     call pack_array(frzdepn, mgcols, top_lev, frzdep)
  end if

  ! if (.not. (do_cldice .or. &
  !      (associated(tnd_qsnow) .and. associated(tnd_nsnow) .and. associated(re_ice)))) then
  !    errstring = "MG's native cloud ice processes are disabled, but &
  !         &no replacement values were passed in."
  ! end if

  ! if (use_hetfrz_classnuc .and. (.not. &
  !      (associated(frzimm) .and. associated(frzcnt) .and. associated(frzdep)))) then
  !    errstring = "External heterogeneous freezing is enabled, but the &
  !         &required tendencies were not all passed in."
  ! end if

  ! Copies of input concentrations that may be changed internally.
  ! qc = qcn
  ! nc = ncn
  ! qi = qin
  ! ni = nin
  ! qr = qrn
  ! nr = nrn
  ! qs = qsn
  ! ns = nsn

!!== KZ_DCS
  call get_dcst(mgncol,nlev,t,dcst)
!!== KZ_DCS


  ! cldn: used to set cldm, unused for subcolumns
  ! liqcldf: used to set lcldm, unused for subcolumns
  ! icecldf: used to set icldm, unused for subcolumns

  if (microp_uniform) then
     ! subcolumns, set cloud fraction variables to one
     ! if cloud water or ice is present, if not present
     ! set to mincld (mincld used instead of zero, to prevent
     ! possible division by zero errors).

     where (qc >= qsmall)
        lcldm = 1._r8
     elsewhere
        lcldm = mincld
     end where

     where (qi >= qsmall)
        icldm = 1._r8
     elsewhere
        icldm = mincld
     end where

     cldm = max(icldm, lcldm)

  else
     ! get cloud fraction, check for minimum
     cldm = max(cldn,mincld)
     lcldm = max(liqcldf,mincld)
     icldm = max(icecldf,mincld)
  end if

  ! Initialize local variables

  ! local physical properties
  rho = p/(r*t)
  dv = 8.794E-5_r8 * t**1.81_r8 / p
  mu = 1.496E-6_r8 * t**1.5_r8 / (t + 120._r8)
  sc = mu/(rho*dv)

  ! air density adjustment for fallspeed parameters
  ! includes air density correction factor to the
  ! power of 0.54 following Heymsfield and Bansemer 2007

  rhof=(rhosu/rho)**0.54_r8

  arn=ar*rhof
  asn=as*rhof
  acn=g*rhow/(18._r8*mu) * cld_sed
  ain=ai*(rhosu/rho)**0.35_r8

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Get humidity and saturation vapor pressures

  do k=1,nlev
     do i=1,mgncol

        call qsat_water(t(i,k), p(i,k), esl(i,k), qvl(i,k))

        ! make sure when above freezing that esi=esl, not active yet
        if (t(i,k) >= tmelt) then
           esi(i,k)=esl(i,k)
           qvi(i,k)=qvl(i,k)
        else
           call qsat_ice(t(i,k), p(i,k), esi(i,k), qvi(i,k))
        end if

     end do
  end do

  relhum = q / max(qvl, qsmall)

  !===============================================

  ! set mtime here to avoid answer-changing
  mtime=deltat

  ! reset limiter stats
  call zero_stats()

  ! initialize microphysics output
  qcsevap=0._r8
  qisevap=0._r8
  qvres  =0._r8
  cmeitot =0._r8
  vtrmc =0._r8
  vtrmi =0._r8
  qcsedten =0._r8
  qisedten =0._r8
  qrsedten =0._r8
  qssedten =0._r8

  pratot=0._r8
  prctot=0._r8
  mnuccctot=0._r8
  mnuccttot=0._r8
  msacwitot=0._r8
  psacwstot=0._r8
  bergstot=0._r8
  bergtot=0._r8
  melttot=0._r8
  homotot=0._r8
  qcrestot=0._r8
  prcitot=0._r8
  praitot=0._r8
  qirestot=0._r8
  mnuccrtot=0._r8
  pracstot=0._r8
  meltsdttot=0._r8
  frzrdttot=0._r8
  mnuccdtot=0._r8

  rflx=0._r8
  sflx=0._r8

  ! initialize precip output

  qrout=0._r8
  qsout=0._r8
  nrout=0._r8
  nsout=0._r8

  ! for refl calc
  rainrt = 0._r8

  ! initialize rain size
  rercld=0._r8

  qcsinksum_rate1ord = 0._r8

  ! initialize variables for trop_mozart
  nevapr = 0._r8
  prer_evap = 0._r8
  evapsnow = 0._r8
  prain = 0._r8
  prodsnow = 0._r8
  cmeout = 0._r8

  precip_frac = mincld

  lamc=0._r8

  ! initialize microphysical tendencies

  tlat=0._r8
  qvlat=0._r8
  qctend=0._r8
  qitend=0._r8
  qstend = 0._r8
  qrtend = 0._r8
  nctend=0._r8
  nitend=0._r8
  nrtend = 0._r8
  nstend = 0._r8

  ! initialize in-cloud and in-precip quantities to zero
  qcic  = 0._r8
  qiic  = 0._r8
  qsic  = 0._r8
  qric  = 0._r8

  ncic  = 0._r8
  niic  = 0._r8
  nsic  = 0._r8
  nric  = 0._r8

  ! initialize precip at surface

  prect = 0._r8
  preci = 0._r8

  ! initialize precip fallspeeds to zero
  ums = 0._r8
  uns = 0._r8
  umr = 0._r8
  unr = 0._r8

  ! initialize limiter for output
  qcrat = 1._r8

  ! Many outputs have to be initialized here at the top to work around
  ! ifort problems, even if they are always overwritten later.
  effc = 10._r8
  lamcrad = 0._r8
  pgamrad = 0._r8
  effc_fn = 10._r8
  effi = 25._r8
  deffi = 50._r8

  qrout2 = 0._r8
  nrout2 = 0._r8
  drout2 = 0._r8
  qsout2 = 0._r8
  nsout2 = 0._r8
  dsout = 0._r8
  dsout2 = 0._r8

  freqr = 0._r8
  freqs = 0._r8

  reff_rain = 0._r8
  reff_snow = 0._r8

  refl = -9999._r8
  arefl = 0._r8
  areflz = 0._r8
  frefl = 0._r8
  csrfl = 0._r8
  acsrfl = 0._r8
  fcsrfl = 0._r8

  ncal = 0._r8
  ncai = 0._r8

  nfice = 0._r8

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! droplet activation
  ! get provisional droplet number after activation. This is used for
  ! all microphysical process calculations, for consistency with update of
  ! droplet mass before microphysics

  ! calculate potential for droplet activation if cloud water is present
  ! tendency from activation (npccn) is read in from companion routine

  ! output activated liquid and ice (convert from #/kg -> #/m3)
  !--------------------------------------------------
  where (qc >= qsmall)
     nc = max(nc + npccn*deltat, 0._r8)
     ncal = nc*rho/lcldm ! sghan minimum in #/cm3
  elsewhere
     ncal = 0._r8
  end where

  where (t < icenuct)
     ncai = naai*rho
  elsewhere
     ncai = 0._r8
  end where

  !===============================================

  ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
  !-------------------------------------------------------

  if (do_cldice) then
     where (naai > 0._r8 .and. t < icenuct .and. &
          relhum*esl/esi > rhmini+0.05_r8)

        !if NAAI > 0. then set numice = naai (as before)
        !note: this is gridbox averaged
        nnuccd = (naai-ni/icldm)/mtime*icldm
        nnuccd = max(nnuccd,0._r8)
        nimax = naai*icldm

        !Calc mass of new particles using new crystal mass...
        !also this will be multiplied by mtime as nnuccd is...

        mnuccd = nnuccd * mi0

     elsewhere
        nnuccd = 0._r8
        nimax = 0._r8
        mnuccd = 0._r8
     end where

  end if


  !=============================================================================
  pre_vert_loop: do k=1,nlev

     pre_col_loop: do i=1,mgncol

        ! calculate instantaneous precip processes (melting and homogeneous freezing)

        ! melting of snow at +2 C

        if (t(i,k) > snowmelt) then
           if (qs(i,k) > 0._r8) then

              ! make sure melting snow doesn't reduce temperature below threshold
              dum = -xlf/cpp*qs(i,k)
              if (t(i,k)+dum < snowmelt) then
                 dum = (t(i,k)-snowmelt)*cpp/xlf
                 dum = dum/qs(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              minstsm(i,k) = dum*qs(i,k)
              ninstsm(i,k) = dum*ns(i,k)

              dum1=-xlf*minstsm(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1
              meltsdttot(i,k)=meltsdttot(i,k) + dum1

              qs(i,k) = max(qs(i,k) - minstsm(i,k), 0._r8)
              ns(i,k) = max(ns(i,k) - ninstsm(i,k), 0._r8)
              qr(i,k) = max(qr(i,k) + minstsm(i,k), 0._r8)
              nr(i,k) = max(nr(i,k) + ninstsm(i,k), 0._r8)
           end if
        end if

        ! freezing of rain at -5 C

        if (t(i,k) < rainfrze) then

           if (qr(i,k) > 0._r8) then

              ! make sure freezing rain doesn't increase temperature above threshold
              dum = xlf/cpp*qr(i,k)
              if (t(i,k)+dum > rainfrze) then
                 dum = -(t(i,k)-rainfrze)*cpp/xlf
                 dum = dum/qr(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              minstrf(i,k) = dum*qr(i,k)
              ninstrf(i,k) = dum*nr(i,k)

              ! heating tendency
              dum1 = xlf*minstrf(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1
              frzrdttot(i,k)=frzrdttot(i,k) + dum1

              qr(i,k) = max(qr(i,k) - minstrf(i,k), 0._r8)
              nr(i,k) = max(nr(i,k) - ninstrf(i,k), 0._r8)
              qs(i,k) = max(qs(i,k) + minstrf(i,k), 0._r8)
              ns(i,k) = max(ns(i,k) + ninstrf(i,k), 0._r8)

           end if
        end if

        ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
        !-------------------------------------------------------
        ! for microphysical process calculations
        ! units are kg/kg for mixing ratio, 1/kg for number conc

        if (qc(i,k).ge.qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           qcic(i,k)=min(qc(i,k)/lcldm(i,k),5.e-3_r8)
           ncic(i,k)=max(nc(i,k)/lcldm(i,k),0._r8)

           ! specify droplet concentration
           if (nccons) then
              ncic(i,k)=ncnst/rho(i,k)
           end if
        else
           qcic(i,k)=0._r8
           ncic(i,k)=0._r8
        end if

        if (qi(i,k).ge.qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           qiic(i,k)=min(qi(i,k)/icldm(i,k),5.e-3_r8)
           niic(i,k)=max(ni(i,k)/icldm(i,k),0._r8)

           ! switch for specification of cloud ice number
           if (nicons) then
              niic(i,k)=ninst/rho(i,k)
           end if
        else
           qiic(i,k)=0._r8
           niic(i,k)=0._r8
        end if

     end do pre_col_loop
  end do pre_vert_loop

  !========================================================================

  ! for sub-columns cldm has already been set to 1 if cloud
  ! water or ice is present, so precip_frac will be correctly set below
  ! and nothing extra needs to be done here

  precip_frac = cldm

  micro_vert_loop: do k=1,nlev

     if (trim(micro_mg_precip_frac_method) == 'in_cloud') then

        if (k /= 1) then
           where (qc(:,k) < qsmall .and. qi(:,k) < qsmall)
              precip_frac(:,k) = precip_frac(:,k-1)
           end where
        endif

     else if (trim(micro_mg_precip_frac_method) == 'max_overlap') then

        ! calculate precip fraction based on maximum overlap assumption

        ! if rain or snow mix ratios are smaller than threshold,
        ! then leave precip_frac as cloud fraction at current level
        if (k /= 1) then
           where (qr(:,k-1) >= qsmall .or. qs(:,k-1) >= qsmall)
              precip_frac(:,k)=max(precip_frac(:,k-1),precip_frac(:,k))
           end where
        end if

     endif

     do i = 1, mgncol

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! get size distribution parameters based on in-cloud cloud water
        ! these calculations also ensure consistency between number and mixing ratio
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ! cloud liquid
        !-------------------------------------------

        call size_dist_param_liq(mg_liq_props, qcic(i,k), ncic(i,k), rho(i,k), &
             pgam(i,k), lamc(i,k))

     end do

     !========================================================================
     ! autoconversion of cloud liquid water to rain
     ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
     ! minimum qc of 1 x 10^-8 prevents floating point error

     call kk2000_liq_autoconversion(microp_uniform, qcic(:,k), &
          ncic(:,k), rho(:,k), relvar(:,k),mg_prc_coeff_fix,prc_coef1,prc_exp,prc_exp1, prc(:,k), nprc(:,k), nprc1(:,k))

     ! assign qric based on prognostic qr, using assumed precip fraction
     ! note: this could be moved above for consistency with qcic and qiic calculations
     qric(:,k) = qr(:,k)/precip_frac(:,k)
     nric(:,k) = nr(:,k)/precip_frac(:,k)

     ! LIMITER: flag and save change due to limiter
     where (qric(:,k) > 0.01_r8)
        qric_limiter(:,k)     = 1
        qric_limiter_mag(:,k) = (qric(:,k) - 0.01_r8)/qric(:,k)
     elsewhere
        qric_limiter(:,k)     = 0
        qric_limiter_mag(:,k) = 0.0d0
     end where

     ! limit in-precip mixing ratios to 10 g/kg at most
     qric(:,k)=min(qric(:,k),0.01_r8)

     ! add autoconversion to precip from above to get provisional rain mixing ratio
     ! and number concentration (qric and nric)

     ! LIMITER: flag and save change due to limiter
     where (qric(:,k) < qsmall)
        qric_qsmall(:,k) = 1
        nric_qsmall(:,k) = 1

        qric_qsmall_mag(:,k) = abs(qric(:,k))
        nric_qsmall_mag(:,k) = abs(nric(:,k))
     elsewhere
        qric_qsmall(:,k) = 0
        nric_qsmall(:,k) = 0

        qric_qsmall_mag(:,k) = 0.0d0
        nric_qsmall_mag(:,k) = 0.0d0
     end where

     where (qric(:,k).lt.qsmall)
        qric(:,k)=0._r8
        nric(:,k)=0._r8
     end where

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

     ! LIMITER: flag and save change due to limiter
     where (nric(:,k) < 0.0_r8)
        nric_neg(:,k)     = 1
        nric_neg_mag(:,k) = abs(nric(:,k))
     elsewhere
        nric_neg(:,k)     = 0
        nric_neg_mag(:,k) = 0.0d0
     end where

     nric(:,k)=max(nric(:,k),0._r8)

!!== KZ_DCS
     if(dcs_tdep) then
        ! Get size distribution parameters for cloud ice
        call size_dist_param_ice(mg_ice_props, dcst(:,k), qiic(:,k), niic(:,k), &
             lami(:,k), n0i(:,k))
     else
        ! Get size distribution parameters for cloud ice
        call size_dist_param_basic(mg_ice_props, qiic(:,k), niic(:,k), &
             lami(:,k), n0i(:,k))
     end if
!!== KZ_DCS

     !.......................................................................
     ! Autoconversion of cloud ice to snow
     ! similar to Ferrier (1994)

     if (do_cldice) then
        call ice_autoconversion(t(:,k), qiic(:,k), lami(:,k), n0i(:,k), &
!!== KZ_DCS
             dcs, dcst(:,k), dcs_tdep, prci(:,k), nprci(:,k))
!!== KZ_DCS
     else
        ! Add in the particles that we have already converted to snow, and
        ! don't do any further autoconversion of ice.
        prci(:,k)  = tnd_qsnow(:,k) / cldm(:,k)
        nprci(:,k) = tnd_nsnow(:,k) / cldm(:,k)
     end if

     ! note, currently we don't have this
     ! inside the do_cldice block, should be changed later
     ! assign qsic based on prognostic qs, using assumed precip fraction
     qsic(:,k) = qs(:,k)/precip_frac(:,k)
     nsic(:,k) = ns(:,k)/precip_frac(:,k)

     ! LIMITER: flag and save change due to limiter
     where (qsic(:,k) > 0.01_r8)
        qsic_limiter(:,k)     = 1
        qsic_limiter_mag(:,k) = (qsic(:,k) - 0.01_r8)/qsic(:,k)
     elsewhere
        qsic_limiter(:,k)     = 0
        qsic_limiter_mag(:,k) = 0.0d0
     end where

     ! limit in-precip mixing ratios to 10 g/kg at most (>>> limiter <<<)
     qsic(:,k)=min(qsic(:,k),0.01_r8)

     ! if precip mix ratio is zero so should number concentration

     ! LIMITER: flag and save change due to limiter
     where (qsic(:,k) < qsmall)
        qsic_qsmall(:,k) = 1
        nsic_qsmall(:,k) = 1

        qsic_qsmall_mag(:,k) = abs(qsic(:,k))
        nsic_qsmall_mag(:,k) = abs(nsic(:,k))
     elsewhere
        qsic_qsmall(:,k) = 0
        nsic_qsmall(:,k) = 0

        qsic_qsmall_mag(:,k) = 0.0d0
        nsic_qsmall_mag(:,k) = 0.0d0
     end where

     where (qsic(:,k) < qsmall)
        qsic(:,k)=0._r8
        nsic(:,k)=0._r8
     end where

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

     ! LIMITER: flag and save change due to limiter
     where (nsic(:,k) < 0.0_r8)
        nsic_neg(:,k)     = 1
        nsic_neg_mag(:,k) = abs(nsic(:,k))
     elsewhere
        nsic_neg(:,k)     = 0
        nsic_neg_mag(:,k) = 0.0d0
     end where

     nsic(:,k)=max(nsic(:,k),0._r8)

     !.......................................................................
     ! get size distribution parameters for precip
     !......................................................................
     ! rain

     call size_dist_param_basic(mg_rain_props, qric(:,k), nric(:,k), &
          lamr(:,k), n0r(:,k))

     where (lamr(:,k) >= qsmall)

        ! provisional rain number and mass weighted mean fallspeed (m/s)

        unr(:,k) = min(arn(:,k)*gamma_br_plus1/lamr(:,k)**br,9.1_r8*rhof(:,k))
        umr(:,k) = min(arn(:,k)*gamma_br_plus4/(6._r8*lamr(:,k)**br),9.1_r8*rhof(:,k))

     elsewhere
        umr(:,k) = 0._r8
        unr(:,k) = 0._r8
     end where

     !......................................................................
     ! snow

     call size_dist_param_basic(mg_snow_props, qsic(:,k), nsic(:,k), &
          lams(:,k), n0s(:,k))

     where (lams(:,k) > 0._r8)

        ! provisional snow number and mass weighted mean fallspeed (m/s)

        ums(:,k) = min(asn(:,k)*gamma_bs_plus4/(6._r8*lams(:,k)**bs),1.2_r8*rhof(:,k))
        uns(:,k) = min(asn(:,k)*gamma_bs_plus1/lams(:,k)**bs,1.2_r8*rhof(:,k))

     elsewhere
        ums(:,k) = 0._r8
        uns(:,k) = 0._r8
     end where

     if (do_cldice) then
        if (.not. use_hetfrz_classnuc) then

           ! heterogeneous freezing of cloud water
           !----------------------------------------------

           call immersion_freezing(microp_uniform, t(:,k), pgam(:,k), lamc(:,k), &
                qcic(:,k), ncic(:,k), relvar(:,k), mnuccc(:,k), nnuccc(:,k))

           ! make sure number of droplets frozen does not exceed available ice nuclei concentration
           ! this prevents 'runaway' droplet freezing (>>> limiter? <<<)

           where (qcic(:,k).ge.qsmall .and. t(:,k).lt.269.15_r8)
              where (nnuccc(:,k)*lcldm(:,k).gt.nnuccd(:,k))
                 ! scale mixing ratio of droplet freezing with limit
                 mnuccc(:,k)=mnuccc(:,k)*(nnuccd(:,k)/(nnuccc(:,k)*lcldm(:,k)))
                 nnuccc(:,k)=nnuccd(:,k)/lcldm(:,k)
              end where
           end where

           call contact_freezing(microp_uniform, t(:,k), p(:,k), rndst(:,k,:), &
                nacon(:,k,:), pgam(:,k), lamc(:,k), qcic(:,k), ncic(:,k), &
                relvar(:,k), mnucct(:,k), nnucct(:,k))

           mnudep(:,k)=0._r8
           nnudep(:,k)=0._r8

        else

           ! Mass of droplets frozen is the average droplet mass, except
           ! with two limiters: concentration must be at least 1/cm^3, and
           ! mass must be at least the minimum defined above. (>>> limiter? <<<)
           mi0l = qcic(:,k)/max(ncic(:,k), 1.0e6_r8/rho(:,k))
           mi0l = max(mi0l_min, mi0l)

           where (qcic(:,k) >= qsmall)
              nnuccc(:,k) = frzimm(:,k)*1.0e6_r8/rho(:,k)
              mnuccc(:,k) = nnuccc(:,k)*mi0l

              nnucct(:,k) = frzcnt(:,k)*1.0e6_r8/rho(:,k)
              mnucct(:,k) = nnucct(:,k)*mi0l

              nnudep(:,k) = frzdep(:,k)*1.0e6_r8/rho(:,k)
              mnudep(:,k) = nnudep(:,k)*mi0
           elsewhere
              nnuccc(:,k) = 0._r8
              mnuccc(:,k) = 0._r8

              nnucct(:,k) = 0._r8
              mnucct(:,k) = 0._r8

              nnudep(:,k) = 0._r8
              mnudep(:,k) = 0._r8
           end where

        end if

     else
        mnuccc(:,k)=0._r8
        nnuccc(:,k)=0._r8
        mnucct(:,k)=0._r8
        nnucct(:,k)=0._r8
        mnudep(:,k)=0._r8
        nnudep(:,k)=0._r8
     end if

     call snow_self_aggregation(t(:,k), rho(:,k), asn(:,k), rhosn, qsic(:,k), nsic(:,k), &
          nsagg(:,k))

     call accrete_cloud_water_snow(t(:,k), rho(:,k), asn(:,k), uns(:,k), mu(:,k), &
          qcic(:,k), ncic(:,k), qsic(:,k), pgam(:,k), lamc(:,k), lams(:,k), n0s(:,k), &
          psacws(:,k), npsacws(:,k))

     if (do_cldice) then
        call secondary_ice_production(t(:,k), psacws(:,k), msacwi(:,k), nsacwi(:,k))
     else
        nsacwi(:,k) = 0.0_r8
        msacwi(:,k) = 0.0_r8
     end if

     call accrete_rain_snow(t(:,k), rho(:,k), umr(:,k), ums(:,k), unr(:,k), uns(:,k), &
          qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), &
          pracs(:,k), npracs(:,k))

     call heterogeneous_rain_freezing(t(:,k), qric(:,k), nric(:,k), lamr(:,k), &
          mnuccr(:,k), nnuccr(:,k))

     call accrete_cloud_water_rain(microp_uniform, qric(:,k), qcic(:,k), &
          ncic(:,k), relvar(:,k), accre_enhan(:,k), pra(:,k), npra(:,k))

     call self_collection_rain(rho(:,k), qric(:,k), nric(:,k), nragg(:,k))

     if (do_cldice) then
        call accrete_cloud_ice_snow(t(:,k), rho(:,k), asn(:,k), qiic(:,k), niic(:,k), &
             qsic(:,k), lams(:,k), n0s(:,k), prai(:,k), nprai(:,k))
     else
        prai(:,k) = 0._r8
        nprai(:,k) = 0._r8
     end if

     call evaporate_sublimate_precip(t(:,k), rho(:,k), &
          dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k), &
          lcldm(:,k), precip_frac(:,k), arn(:,k), asn(:,k), qcic(:,k), qiic(:,k), &
          qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), &
          pre(:,k), prds(:,k))

     call bergeron_process_snow(t(:,k), rho(:,k), dv(:,k), mu(:,k), sc(:,k), &
          qvl(:,k), qvi(:,k), asn(:,k), qcic(:,k), qsic(:,k), lams(:,k), n0s(:,k), &
          bergs(:,k))

     bergs(:,k)=bergs(:,k)*micro_mg_berg_eff_factor

     !+++PMC 12/3/12 - NEW VAPOR DEP/SUBLIMATION GOES HERE!!!
     if (do_cldice) then

        call ice_deposition_sublimation(t(:,k), q(:,k), qi(:,k), ni(:,k), &
!!== KZ_DCS
             icldm(:,k), rho(:,k), dv(:,k), qvl(:,k), qvi(:,k), dcst(:,k), dcs_tdep, &
!!== KZ_DCS
             berg(:,k), vap_dep(:,k), ice_sublim(:,k))

        berg(:,k)=berg(:,k)*micro_mg_berg_eff_factor

        where (vap_dep(:,k) < 0._r8 .and. qi(:,k) > qsmall .and. icldm(:,k) > mincld)
           nsubi(:,k) = vap_dep(:,k) / qi(:,k) * ni(:,k) / icldm(:,k)
        elsewhere
           nsubi(:,k) = 0._r8
        end where

        ! bergeron process should not reduce nc unless
        ! all ql is removed (which is handled elsewhere)
        !in fact, nothing in this entire file makes nsubc nonzero.
        nsubc(:,k) = 0._r8

     end if !do_cldice
     !---PMC 12/3/12

     do i=1,mgncol

        ! conservation to ensure no negative values of cloud water/precipitation
        ! in case microphysical process rates are large
        !===================================================================

        ! note: for check on conservation, processes are multiplied by omsm
        ! to prevent problems due to round off error

        ! conservation of qc (>>> limiter <<<)
        !-------------------------------------------------------------------

        dum = ((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+ &
             psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k))*deltat

        if (dum.gt.qc(i,k)) then
           ratio = qc(i,k)/deltat/((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+ &
                msacwi(i,k)+psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k))*omsm
           prc(i,k) = prc(i,k)*ratio
           pra(i,k) = pra(i,k)*ratio
           mnuccc(i,k) = mnuccc(i,k)*ratio
           mnucct(i,k) = mnucct(i,k)*ratio
           msacwi(i,k) = msacwi(i,k)*ratio
           psacws(i,k) = psacws(i,k)*ratio
           bergs(i,k) = bergs(i,k)*ratio
           berg(i,k) = berg(i,k)*ratio
           qcrat(i,k) = ratio

           ! LIMITER: flag and save change due to limiter
           qc_conservation(i,k)     = 1
           qc_conservation_mag(i,k) = deltat*(ratio-1.d0)* &
                (prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+ &
                    msacwi(i,k)+psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k)

        else
           qcrat(i,k) = 1._r8

           ! LIMITER: flag and save change due to limiter
           qc_conservation(i,k)     = 0
           qc_conservation_mag(i,k) = 0.0d0
        end if

        !PMC 12/3/12: ratio is also frac of step w/ liquid.
        !thus we apply berg for "ratio" of timestep and vapor
        !deposition for the remaining frac of the timestep.
        if (qc(i,k) >= qsmall) then
           vap_dep(i,k) = vap_dep(i,k)*(1._r8-qcrat(i,k))
        end if

     end do

     do i=1,mgncol

        !=================================================================
        ! apply limiter to ensure that ice/snow sublimation and rain evap
        ! don't push conditions into supersaturation, and ice deposition/nucleation don't
        ! push conditions into sub-saturation
        ! note this is done after qc conservation since we don't know how large
        ! vap_dep is before then
        ! estimates are only approximate since other process terms haven't been limited
        ! for conservation yet

        ! first limit ice deposition/nucleation vap_dep + mnuccd
        dum1 = vap_dep(i,k) + mnuccd(i,k)
        if (dum1 > 1.e-20_r8) then
           dum = (q(i,k)-qvi(i,k))/(1._r8 + xxls_squared*qvi(i,k)/(cpp*rv*t(i,k)**2))/deltat
           dum = max(dum,0._r8)
           if (dum1 > dum) then
              ! Allocate the limited "dum" tendency to mnuccd and vap_dep
              ! processes. Don't divide by cloud fraction; these are grid-
              ! mean rates.
              dum1 = mnuccd(i,k) / (vap_dep(i,k)+mnuccd(i,k))

              ! LIMITER: flag and save change due to limiter
              ice_nucleation_limiter(i,k)     = 1
              ice_nucleation_limiter_mag(i,k) = abs(mnuccd(i,k) - dum*dum1)

              mnuccd(i,k) = dum*dum1

              ! LIMITER: flag and save change due to limiter
              ice_deposition_limiter(i,k)     = 1
              ice_deposition_limiter_mag(i,k) = abs(vap_dep(i,k) - dum + mnuccd(i,k))

              vap_dep(i,k) = dum - mnuccd(i,k)

           else
              ! LIMITER: flag and save change due to limiter
              ice_nucleation_limiter(i,k)     = 0
              ice_nucleation_limiter_mag(i,k) = 0.0d0

              ! LIMITER: flag and save change due to limiter
              ice_deposition_limiter(i,k)     = 0
              ice_deposition_limiter_mag(i,k) = 0.0d0
           end if
        end if

     end do

     do i=1,mgncol

        !===================================================================
        ! conservation of nc
        !-------------------------------------------------------------------
        dum = (nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+ &
             npsacws(i,k)-nsubc(i,k))*lcldm(i,k)*deltat

        if (dum.gt.nc(i,k)) then
           ratio = nc(i,k)/deltat/((nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+&
                npsacws(i,k)-nsubc(i,k))*lcldm(i,k))*omsm

           nprc1(i,k) = nprc1(i,k)*ratio
           npra(i,k) = npra(i,k)*ratio
           nnuccc(i,k) = nnuccc(i,k)*ratio
           nnucct(i,k) = nnucct(i,k)*ratio
           npsacws(i,k) = npsacws(i,k)*ratio
           nsubc(i,k)=nsubc(i,k)*ratio

           ! LIMITER: flag and save change due to limiter
           nc_conservation(i,k)     = 1
           nc_conservation_mag(i,k) = deltat*(ratio-1.d0)* &
              (nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+&
                  npsacws(i,k)-nsubc(i,k))*lcldm(i,k)
        else
           ! LIMITER: flag and save change due to limiter
           nc_conservation(i,k)     = 0
           nc_conservation_mag(i,k) = 0.0d0
        end if

        mnuccri(i,k)=0._r8
        nnuccri(i,k)=0._r8

        if (do_cldice) then

           ! freezing of rain to produce ice if mean rain size is smaller than Dcs (???)
           if (lamr(i,k) > qsmall .and. 1._r8/lamr(i,k) < Dcs) then
              mnuccri(i,k)=mnuccr(i,k)
              nnuccri(i,k)=nnuccr(i,k)
              mnuccr(i,k)=0._r8
              nnuccr(i,k)=0._r8
           end if
        end if

     end do

     do i=1,mgncol

        ! conservation of rain mixing ratio
        !-------------------------------------------------------------------
        dum = ((-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k))*precip_frac(i,k)- &
             (pra(i,k)+prc(i,k))*lcldm(i,k))*deltat

        ! note that qrtend is included below because of instantaneous freezing/melt
        if (dum.gt.qr(i,k).and. &
             (-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k)).ge.qsmall) then
           ratio = (qr(i,k)/deltat+(pra(i,k)+prc(i,k))*lcldm(i,k))/   &
                precip_frac(i,k)/(-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k))*omsm
           pre(i,k)=pre(i,k)*ratio
           pracs(i,k)=pracs(i,k)*ratio
           mnuccr(i,k)=mnuccr(i,k)*ratio
           mnuccri(i,k)=mnuccri(i,k)*ratio

           ! LIMITER: flag and save change due to limiter
           qr_conservation(i,k)     = 1
           qr_conservation_mag(i,k) = deltat*(ratio-1.d0)* &
              precip_frac(i,k)*(-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k))
        else
           ! LIMITER: flag and save change due to limiter
           qr_conservation(i,k)     = 0
           qr_conservation_mag(i,k) = 0.0d0
        end if

     end do

     do i=1,mgncol

        ! conservation of rain number (>>> limiter? <<< not really)
        !-------------------------------------------------------------------

        ! Add evaporation of rain number.
        if (pre(i,k) < 0._r8) then
           dum = pre(i,k)*deltat/qr(i,k)
           dum = max(-1._r8,dum)
           nsubr(i,k) = dum*nr(i,k)/deltat
        else
           nsubr(i,k) = 0._r8
        end if

     end do

     do i=1,mgncol

        dum = ((-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k))*precip_frac(i,k)- &
             nprc(i,k)*lcldm(i,k))*deltat

        if (dum.gt.nr(i,k)) then
           ratio = (nr(i,k)/deltat+nprc(i,k)*lcldm(i,k))/precip_frac(i,k)/ &
                (-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k))*omsm

           nragg(i,k)=nragg(i,k)*ratio
           npracs(i,k)=npracs(i,k)*ratio
           nnuccr(i,k)=nnuccr(i,k)*ratio
           nsubr(i,k)=nsubr(i,k)*ratio
           nnuccri(i,k)=nnuccri(i,k)*ratio

           ! LIMITER: flag and save change due to limiter
           nr_conservation(i,k)     = 1
           nr_conservation_mag(i,k) = deltat*(ratio-1.d0)* &
              (-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k))

        else
           ! LIMITER: flag and save change due to limiter
           nr_conservation(i,k)     = 0
           nr_conservation_mag(i,k) = 0.0d0
        end if

     end do

     if (do_cldice) then

        do i=1,mgncol

           ! conservation of qi
           !-------------------------------------------------------------------

           dum = ((-mnuccc(i,k)-mnucct(i,k)-mnudep(i,k)-msacwi(i,k))*lcldm(i,k)+(prci(i,k)+ &
                prai(i,k))*icldm(i,k)-mnuccri(i,k)*precip_frac(i,k) &
                -ice_sublim(i,k)-vap_dep(i,k)-berg(i,k)-mnuccd(i,k))*deltat

           if (dum.gt.qi(i,k)) then
              ratio = (qi(i,k)/deltat+vap_dep(i,k)+berg(i,k)+mnuccd(i,k)+ &
                   (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k))*lcldm(i,k)+ &
                   mnuccri(i,k)*precip_frac(i,k))/ &
                   ((prci(i,k)+prai(i,k))*icldm(i,k)-ice_sublim(i,k))*omsm
              prci(i,k) = prci(i,k)*ratio
              prai(i,k) = prai(i,k)*ratio
              ice_sublim(i,k) = ice_sublim(i,k)*ratio

              ! LIMITER: flag and save change due to limiter
              qi_conservation(i,k)     = 1
              qi_conservation_mag(i,k) = ratio
           else
              ! LIMITER: flag and save change due to limiter
              qi_conservation(i,k)     = 0
              qi_conservation_mag(i,k) = 1.0d0
           end if

        end do

     end if

     if (do_cldice) then

        do i=1,mgncol

           ! conservation of ni
           !-------------------------------------------------------------------
           if (use_hetfrz_classnuc) then
              tmpfrz = nnuccc(i,k)
           else
              tmpfrz = 0._r8
           end if
           dum = ((-nnucct(i,k)-tmpfrz-nnudep(i,k)-nsacwi(i,k))*lcldm(i,k)+(nprci(i,k)+ &
                nprai(i,k)-nsubi(i,k))*icldm(i,k)-nnuccri(i,k)*precip_frac(i,k)- &
                nnuccd(i,k))*deltat

           if (dum.gt.ni(i,k)) then
              ratio = (ni(i,k)/deltat+nnuccd(i,k)+ &
                   (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k))*lcldm(i,k)+ &
                   nnuccri(i,k)*precip_frac(i,k))/ &
                   ((nprci(i,k)+nprai(i,k)-nsubi(i,k))*icldm(i,k))*omsm
              nprci(i,k) = nprci(i,k)*ratio
              nprai(i,k) = nprai(i,k)*ratio
              nsubi(i,k) = nsubi(i,k)*ratio

              ! LIMITER: flag and save change due to limiter
              ni_conservation(i,k)     = 1
              ni_conservation_mag(i,k) = ratio
           else
              ! LIMITER: flag and save change due to limiter
              ni_conservation(i,k)     = 0
              ni_conservation_mag(i,k) = 1.0d0
           end if

        end do

     end if

     do i=1,mgncol

        ! conservation of snow mixing ratio
        !-------------------------------------------------------------------
        dum = (-(prds(i,k)+pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)-(prai(i,k)+prci(i,k))*icldm(i,k) &
             -(bergs(i,k)+psacws(i,k))*lcldm(i,k))*deltat

        if (dum.gt.qs(i,k).and.-prds(i,k).ge.qsmall) then
           ratio = (qs(i,k)/deltat+(prai(i,k)+prci(i,k))*icldm(i,k)+ &
                (bergs(i,k)+psacws(i,k))*lcldm(i,k)+(pracs(i,k)+mnuccr(i,k))*precip_frac(i,k))/ &
                precip_frac(i,k)/(-prds(i,k))*omsm
           prds(i,k)=prds(i,k)*ratio

           ! LIMITER: flag and save change due to limiter
           qs_conservation(i,k)     = 1
           qs_conservation_mag(i,k) = ratio
        else
           ! LIMITER: flag and save change due to limiter
           qs_conservation(i,k)     = 0
           qs_conservation_mag(i,k) = 1.0d0
        end if

     end do

     do i=1,mgncol

        ! conservation of snow number
        !-------------------------------------------------------------------
        ! calculate loss of number due to sublimation
        ! for now neglect sublimation of ns
        nsubs(i,k)=0._r8

        dum = ((-nsagg(i,k)-nsubs(i,k)-nnuccr(i,k))*precip_frac(i,k)-nprci(i,k)*icldm(i,k))*deltat

        if (dum.gt.ns(i,k)) then
           ratio = (ns(i,k)/deltat+nnuccr(i,k)* &
                precip_frac(i,k)+nprci(i,k)*icldm(i,k))/precip_frac(i,k)/ &
                (-nsubs(i,k)-nsagg(i,k))*omsm

           ! BUG: if nsubs and nsagg are zero then the ratio is infinite
           ! write(*,*) '------'
           ! write(*,*) dum, ns(i,k), deltat, nnuccr(i,k), precip_frac(i,k), &
           !      nprci(i,k), icldm(i,k), precip_frac(i,k), &
           !      nsubs(i,k), nsagg(i,k), omsm, ratio

           nsubs(i,k)=nsubs(i,k)*ratio
           nsagg(i,k)=nsagg(i,k)*ratio

           ! BUG: if nsubs and nsagg are zero then the ratio is infinite
           ! write(*,*) nsubs(i,k), nsagg(i,k)

           ! LIMITER: flag and save change due to limiter
           ns_conservation(i,k)     = 1
           ns_conservation_mag(i,k) = ratio
        else
           ! LIMITER: flag and save change due to limiter
           ns_conservation(i,k)     = 0
           ns_conservation_mag(i,k) = 1.0d0
        end if

     end do

     do i=1,mgncol

        ! DJG NOTE: These limiters appeared after qc limiter in older MG2

        ! next limit ice and snow sublimation and rain evaporation
        ! get estimate of q and t at end of time step
        ! don't include other microphysical processes since they haven't
        ! been limited via conservation checks yet

        if ((pre(i,k)+prds(i,k))*precip_frac(i,k)+ice_sublim(i,k) < -1.e-20_r8) then

           qtmp=q(i,k)-(ice_sublim(i,k)+vap_dep(i,k)+mnuccd(i,k)+ &
                (pre(i,k)+prds(i,k))*precip_frac(i,k))*deltat
           ttmp=t(i,k)+((pre(i,k)*precip_frac(i,k))*xxlv+ &
                (prds(i,k)*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k))*xxls)*deltat/cpp

           ! use rhw to allow ice supersaturation
           call qsat_water(ttmp, p(i,k), esn, qvn)

           ! modify ice/precip evaporation rate if q > qsat
           if (qtmp > qvn) then

              ! LIMITER: flag limiter, see changes due to limiter below
              qiqs_sublimation_qr_evaporation_limiter(i,k) = 1

              dum1=pre(i,k)*precip_frac(i,k)/((pre(i,k)+prds(i,k))*precip_frac(i,k)+ice_sublim(i,k))
              dum2=prds(i,k)*precip_frac(i,k)/((pre(i,k)+prds(i,k))*precip_frac(i,k)+ice_sublim(i,k))
              ! recalculate q and t after vap_dep and mnuccd but without evap or sublim
              qtmp=q(i,k)-(vap_dep(i,k)+mnuccd(i,k))*deltat
              ttmp=t(i,k)+((vap_dep(i,k)+mnuccd(i,k))*xxls)*deltat/cpp

              ! use rhw to allow ice supersaturation
              call qsat_water(ttmp, p(i,k), esn, qvn)

              dum=(qtmp-qvn)/(1._r8 + xxlv_squared*qvn/(cpp*rv*ttmp**2))
              dum=min(dum,0._r8)

              ! LIMITER: save change due to limiter
              rain_evaporation_limiter_mag(i,k) = &
                   abs(pre(i,k) - dum*dum1/deltat/precip_frac(i,k))

              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              pre(i,k)=dum*dum1/deltat/precip_frac(i,k)

              ! do separately using RHI for prds and ice_sublim
              call qsat_ice(ttmp, p(i,k), esn, qvn)

              dum=(qtmp-qvn)/(1._r8 + xxls_squared*qvn/(cpp*rv*ttmp**2))
              dum=min(dum,0._r8)

              ! LIMITER: save change due to limiter
              snow_sublimation_limiter_mag(i,k) = &
                   abs(prds(i,k) - dum*dum2/deltat/precip_frac(i,k))

              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              prds(i,k) = dum*dum2/deltat/precip_frac(i,k)

              ! LIMITER: save change due to limiter
              ice_sublimation_limiter_mag(i,k) = &
                   abs(ice_sublim(i,k) - dum*(1._r8-dum1-dum2)/deltat)

              ! don't divide ice_sublim by cloud fraction since it is grid-averaged
              dum1 = (1._r8-dum1-dum2)
              ice_sublim(i,k) = dum*dum1/deltat
           else

              ! LIMITER: flag and save change due to limiter
              qiqs_sublimation_qr_evaporation_limiter(i,k) = 0
              rain_evaporation_limiter_mag(i,k) = 0.0d0
              snow_sublimation_limiter_mag(i,k) = 0.0d0
              ice_sublimation_limiter_mag(i,k)  = 0.0d0

           end if
        end if

     end do

     ! End Big "administration" loop enforces conservation, updates variables
     ! that accumulate over substeps, and sets output variables.

     do i=1,mgncol

        ! get tendencies due to microphysical conversion processes
        !==========================================================
        ! note: tendencies are multiplied by appropriate cloud/precip
        ! fraction to get grid-scale values
        ! note: vap_dep is already grid-average values

        ! The net tendencies need to be added to rather than overwritten,
        ! because they may have a value already set for instantaneous
        ! melting/freezing.

        qvlat(i,k) = qvlat(i,k)-(pre(i,k)+prds(i,k))*precip_frac(i,k)-&
             vap_dep(i,k)-ice_sublim(i,k)-mnuccd(i,k)-mnudep(i,k)*lcldm(i,k)

        tlat(i,k) = tlat(i,k)+((pre(i,k)*precip_frac(i,k)) &
             *xxlv+(prds(i,k)*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k)+mnudep(i,k)*lcldm(i,k))*xxls+ &
             ((bergs(i,k)+psacws(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k))*lcldm(i,k)+(mnuccr(i,k)+ &
             pracs(i,k)+mnuccri(i,k))*precip_frac(i,k)+berg(i,k))*xlf)

        qctend(i,k) = qctend(i,k)+ &
             (-pra(i,k)-prc(i,k)-mnuccc(i,k)-mnucct(i,k)-msacwi(i,k)- &
             psacws(i,k)-bergs(i,k))*lcldm(i,k)-berg(i,k)

        if (do_cldice) then
           qitend(i,k) = qitend(i,k)+ &
                (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k))*lcldm(i,k)+(-prci(i,k)- &
                prai(i,k))*icldm(i,k)+vap_dep(i,k)+berg(i,k)+ice_sublim(i,k)+ &
                mnuccd(i,k)+mnuccri(i,k)*precip_frac(i,k)
        end if

        qrtend(i,k) = qrtend(i,k)+ &
             (pra(i,k)+prc(i,k))*lcldm(i,k)+(pre(i,k)-pracs(i,k)- &
             mnuccr(i,k)-mnuccri(i,k))*precip_frac(i,k)

        qstend(i,k) = qstend(i,k)+ &
             (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(prds(i,k)+ &
             pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)


        cmeout(i,k) = vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)

        ! add output for cmei (accumulate)
        cmeitot(i,k) = vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)

        ! assign variables for trop_mozart, these are grid-average
        !-------------------------------------------------------------------
        ! evaporation/sublimation is stored here as positive term

        evapsnow(i,k) = -prds(i,k)*precip_frac(i,k)
        nevapr(i,k) = -pre(i,k)*precip_frac(i,k)
        prer_evap(i,k) = -pre(i,k)*precip_frac(i,k)

        ! change to make sure prain is positive: do not remove snow from
        ! prain used for wet deposition
        prain(i,k) = (pra(i,k)+prc(i,k))*lcldm(i,k)+(-pracs(i,k)- &
             mnuccr(i,k)-mnuccri(i,k))*precip_frac(i,k)
        prodsnow(i,k) = (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(&
             pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)

        ! following are used to calculate 1st order conversion rate of cloud water
        !    to rain and snow (1/s), for later use in aerosol wet removal routine
        ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
        !    used to calculate pra, prc, ... in this routine
        ! qcsinksum_rate1ord = { rate of direct transfer of cloud water to rain & snow }
        !                      (no cloud ice or bergeron terms)
        qcsinksum_rate1ord(i,k) = (pra(i,k)+prc(i,k)+psacws(i,k))*lcldm(i,k)
        ! Avoid zero/near-zero division.
        qcsinksum_rate1ord(i,k) = qcsinksum_rate1ord(i,k) / &
             max(qc(i,k),1.0e-30_r8)


        ! microphysics output, note this is grid-averaged
        pratot(i,k) = pra(i,k)*lcldm(i,k)
        prctot(i,k) = prc(i,k)*lcldm(i,k)
        mnuccctot(i,k) = mnuccc(i,k)*lcldm(i,k)
        mnuccttot(i,k) = mnucct(i,k)*lcldm(i,k)
        msacwitot(i,k) = msacwi(i,k)*lcldm(i,k)
        psacwstot(i,k) = psacws(i,k)*lcldm(i,k)
        bergstot(i,k) = bergs(i,k)*lcldm(i,k)
        bergtot(i,k) = berg(i,k)
        prcitot(i,k) = prci(i,k)*icldm(i,k)
        praitot(i,k) = prai(i,k)*icldm(i,k)
        mnuccdtot(i,k) = mnuccd(i,k)*icldm(i,k)

        pracstot(i,k) = pracs(i,k)*precip_frac(i,k)
        mnuccrtot(i,k) = mnuccr(i,k)*precip_frac(i,k)


        nctend(i,k) = nctend(i,k)+&
             (-nnuccc(i,k)-nnucct(i,k)-npsacws(i,k)+nsubc(i,k) &
             -npra(i,k)-nprc1(i,k))*lcldm(i,k)

        if (do_cldice) then
           if (use_hetfrz_classnuc) then
              tmpfrz = nnuccc(i,k)
           else
              tmpfrz = 0._r8
           end if
           nitend(i,k) = nitend(i,k)+ nnuccd(i,k)+ &
                (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k))*lcldm(i,k)+(nsubi(i,k)-nprci(i,k)- &
                nprai(i,k))*icldm(i,k)+nnuccri(i,k)*precip_frac(i,k)
        end if

        nstend(i,k) = nstend(i,k)+(nsubs(i,k)+ &
             nsagg(i,k)+nnuccr(i,k))*precip_frac(i,k)+nprci(i,k)*icldm(i,k)

        nrtend(i,k) = nrtend(i,k)+ &
             nprc(i,k)*lcldm(i,k)+(nsubr(i,k)-npracs(i,k)-nnuccr(i,k) &
             -nnuccri(i,k)+nragg(i,k))*precip_frac(i,k)

        ! make sure that ni at advanced time step does not exceed
        ! maximum (existing N + source terms*dt), which is possible if mtime < deltat
        ! note that currently mtime = deltat
        !================================================================

        if (do_cldice .and. nitend(i,k).gt.0._r8.and.ni(i,k)+nitend(i,k)*deltat.gt.nimax(i,k)) then
           ! LIMITER: flag and save change due to limiter
           ni_tendency_limiter(i,k)     = 1
           ni_tendency_limiter_mag(i,k) = &
                abs((nitend(i,k) - max(0._r8,(nimax(i,k)-ni(i,k))/deltat))/nitend(i,k))

           nitend(i,k)=max(0._r8,(nimax(i,k)-ni(i,k))/deltat)
        else
           ! LIMITER: flag and save change due to limiter
           ni_tendency_limiter(i,k)     = 0
           ni_tendency_limiter_mag(i,k) = 0.0d0
        end if

     end do

     ! End of "administration" loop

  end do micro_vert_loop ! end k loop

  !-----------------------------------------------------
  ! convert rain/snow q and N for output to history, note,
  ! output is for gridbox average

  qrout = qr
  nrout = nr * rho
  qsout = qs
  nsout = ns * rho

  ! calculate precip fluxes
  ! calculate the precip flux (kg/m2/s) as mixingratio(kg/kg)*airdensity(kg/m3)*massweightedfallspeed(m/s)
  ! ---------------------------------------------------------------------

  rflx(:,2:) = rflx(:,2:) + (qric*rho*umr*precip_frac)
  sflx(:,2:) = sflx(:,2:) + (qsic*rho*ums*precip_frac)

  ! calculate n0r and lamr from rain mass and number
  ! divide by precip fraction to get in-precip (local) values of
  ! rain mass and number, divide by rhow to get rain number in kg^-1

  call size_dist_param_basic(mg_rain_props, qric, nric, lamr, n0r)

  ! Calculate rercld

  ! calculate mean size of combined rain and cloud water

  call calc_rercld(lamr, n0r, lamc, pgam, qric, qcic, ncic, &
       rercld)


  ! Assign variables back to start-of-timestep values
  ! Some state variables are changed before the main microphysics loop
  ! to make "instantaneous" adjustments. Afterward, we must move those changes
  ! back into the tendencies.
  ! These processes:
  !  - Droplet activation (npccn, impacts nc)
  !  - Instantaneous snow melting  (minstsm/ninstsm, impacts qr/qs/nr/ns)
  !  - Instantaneous rain freezing (minstfr/ninstrf, impacts qr/qs/nr/ns)
  !================================================================================

  ! Re-apply droplet activation tendency
  call pack_array(ncn, mgcols, top_lev, nc)
  nctend = nctend + npccn

  ! Re-apply rain freezing and snow melting.
  dum_2D = qs
  call pack_array(qsn, mgcols, top_lev, qs)
  qstend = qstend + (dum_2D-qs)/deltat

  dum_2D = ns
  call pack_array(nsn, mgcols, top_lev, ns)
  nstend = nstend + (dum_2D-ns)/deltat

  dum_2D = qr
  call pack_array(qrn, mgcols, top_lev, qr)
  qrtend = qrtend + (dum_2D-qr)/deltat

  dum_2D = nr
  call pack_array(nrn, mgcols, top_lev, nr)
  nrtend = nrtend + (dum_2D-nr)/deltat

  !.............................................................................

  !================================================================================

  ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
  nevapr = nevapr + evapsnow
  prain = prain + prodsnow

  sed_col_loop: do i=1,mgncol

    do k=1,nlev

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation

        dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
        dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat),0._r8)
        dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
        dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat),0._r8)
        dumr(i,k) = (qr(i,k)+qrtend(i,k)*deltat)
        dumnr(i,k) = max((nr(i,k)+nrtend(i,k)*deltat),0._r8)
        dums(i,k) = (qs(i,k)+qstend(i,k)*deltat)
        dumns(i,k) = max((ns(i,k)+nstend(i,k)*deltat),0._r8)

     end do       !!! vertical loop

     where (dumc(i,:) .lt. qsmall) dumc(i,:) = 0._r8
     where (dumi(i,:) .lt. qsmall) dumi(i,:) = 0._r8
     where (dumr(i,:) .lt. qsmall) dumr(i,:) = 0._r8
     where (dums(i,:) .lt. qsmall) dums(i,:) = 0._r8

     ! initialize nstep for sedimentation sub-steps
     fid_nstep = 102
     fname = trim(fileNameOut)//'.nstep.txt'
     if (mphys_calls == 0) then
        open(unit = fid_nstep, file = fname, status="new", action="write")
        write(fid_nstep,'(A12,4X,A12,4X,A12,4X,A12,4X,A12)') "mphys_call", &
             "nstep_qi", "nstep_qc", "nstep_qr","nstep_qs"
     else
        open(unit = fid_nstep, file = fname, status="old", &
             position="append", action="write")
     end if

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     if (do_cldice) then
       deltat_sed = deltat
       time_sed = deltat
     else
       time_sed = 0._r8
     end if
     nstep = 0

     ! subcycle
     do while (time_sed > qsmall)
!       gptl_ret = gptlstart_handle('Sedimentation (ice)',gptl_sed_ice)
       ! obtain CFL number
#ifndef SED_UPDATECFL
       if (nstep == 0) then
#endif
       gptl_ret = gptlstart_handle('CFL Calculation (ice)',gptl_cfl_ice)
       do n=1,gptl_loop
         call sed_CalcCFL(dumi,dumni,icldm,rho,pdel,nlev,i,&
           MG_ICE,deltat_sed,g,ain,rhof,alphaq,alphan,s1,s2,w1,w2,dum,&
           ncons=nicons,nnst=ninst,&
           gamma_b_plus1=gamma_bi_plus1,gamma_b_plus4=gamma_bi_plus4)
       end do
       gptl_ret = gptlstop_handle('CFL Calculation (ice)',gptl_cfl_ice)
#ifndef SED_UPDATECFL
      else
        dum = CFL
      end if
#endif
       ! update deltat_sed for target CFL number
       deltat_sed = min(deltat_sed*CFL/dum,time_sed)
       ! advance cloud ice sedimentation
       time_sed = time_sed - deltat_sed
       nstep = nstep+1
       gptl_ret = gptlstart_handle('Flux Calculation (ice)',gptl_flux_ice)
       do n=1,gptl_loop
         call sed_CalcFlux(dumi,dumni,alphaq,alphan,s1,s2,w1,w2,nlev,i,MG_ICE,fi,fni)
       end do
       gptl_ret = gptlstop_handle('Flux Calculation (ice)',gptl_flux_ice)
       gptl_ret = gptlstart_handle('Advance Solution (ice)',gptl_advance_sol_ice)
       call sed_AdvanceOneStep(dumi,fi,dumni,fni,pdel,deltat,deltat_sed,nlev,i,MG_ICE,g, &
         qitend,nitend,preci,qisedten,&
         cloud_frac=icldm,qvlat=qvlat,tlat=tlat,xxl=xxls,preci=preci,qsevap=qisevap)
       gptl_ret = gptlstop_handle('Advance Solution (ice)',gptl_advance_sol_ice)

!       gptl_ret = gptlstop_handle('Sedimentation (ice)',gptl_sed_ice)
     end do

     write(fid_nstep,'(I12,4X,I12,4X)',advance="no") mphys_calls, nstep

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     deltat_sed = deltat
     time_sed = deltat
     nstep = 0

     ! subcycle
     do while (time_sed > qsmall)
!       gptl_ret = gptlstart_handle('Sedimentation (cloud)',gptl_sed_cloud)
       ! obtain fall speeds
#ifndef SED_UPDATECFL
       if (nstep == 0) then
#endif
       gptl_ret = gptlstart_handle('CFL Calculation (cloud)',gptl_cfl_cloud)
       do n=1,gptl_loop
         call sed_CalcCFL(dumc,dumnc,lcldm,rho,pdel,nlev,i,&
           MG_LIQUID,deltat_sed,g,acn,rhof,alphaq,alphan,s1,s2,w1,w2,dum,&
           ncons=nccons,nnst=ncnst)
       end do
       gptl_ret = gptlstop_handle('CFL Calculation (cloud)',gptl_cfl_cloud)
#ifndef SED_UPDATECFL
       else
         dum = CFL
       end if
#endif
       ! update deltat_sed for target CFL number
       deltat_sed = min(deltat_sed*CFL/dum,time_sed)
       ! advance cloud liquid sedimentation
       time_sed = time_sed - deltat_sed
       nstep = nstep + 1
       gptl_ret = gptlstart_handle('Flux Calculation (cloud)',gptl_flux_cloud)
       do n=1,gptl_loop
         call sed_CalcFlux(dumc,dumnc,alphaq,alphan,s1,s2,w1,w2,nlev,i,MG_LIQUID,fc,fnc)
       end do
       gptl_ret = gptlstop_handle('Flux Calculation (cloud)',gptl_flux_cloud)
       gptl_ret = gptlstart_handle('Advance Solution (cloud)',gptl_advance_sol_cloud)
       call sed_AdvanceOneStep(dumc,fc,dumnc,fnc,pdel,deltat,deltat_sed,nlev,i,MG_LIQUID,g, &
         qctend,nctend,prect,qcsedten,&
         cloud_frac=lcldm,qvlat=qvlat,tlat=tlat,xxl=xxlv,qsevap=qcsevap)
       gptl_ret = gptlstop_handle('Advance Solution (cloud)',gptl_advance_sol_cloud)
!       gptl_ret = gptlstop_handle('Sedimentation (cloud)',gptl_sed_cloud)
     end do

     write(fid_nstep,'(I12,4X)',advance="no") nstep
     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     deltat_sed = deltat
     time_sed = deltat
     nstep = 0

     ! subcycle
     do while (time_sed > qsmall)
!       gptl_ret = gptlstart_handle('Sedimentation (rain)',gptl_sed_rain)
       ! obtain fall speeds
#ifndef SED_UPDATECFL
       if (nstep == 0) then
#endif
       gptl_ret = gptlstart_handle('CFL Calculation (rain)',gptl_cfl_rain)
       do n=1,gptl_loop
         call sed_CalcCFL(dumr,dumnr,precip_frac,rho,pdel,nlev,i,&
           MG_RAIN,deltat_sed,g,arn,rhof,alphaq,alphan,s1,s2,w1,w2,dum,&
           gamma_b_plus1=gamma_br_plus1,gamma_b_plus4=gamma_br_plus4)
       end do
       gptl_ret = gptlstop_handle('CFL Calculation (rain)',gptl_cfl_rain)
#ifndef SED_UPDATECFL
       else
         dum = CFL
       end if
#endif
       ! update deltat_sed for target CFL number
       deltat_sed = min(deltat_sed*CFL/dum,time_sed)
       ! advance rain sedimentation
       time_sed = time_sed - deltat_sed
       nstep = nstep + 1
       gptl_ret = gptlstart_handle('Flux Calculation (rain)',gptl_flux_rain)
       do n=1,gptl_loop
         call sed_CalcFlux(dumr,dumnr,alphaq,alphan,s1,s2,w1,w2,nlev,i,MG_RAIN,fr,fnr)
       end do
       gptl_ret = gptlstop_handle('Flux Calculation (rain)',gptl_flux_rain)
       gptl_ret = gptlstart_handle('Advance Solution (rain)',gptl_advance_sol_rain)
       call sed_AdvanceOneStep(dumr,fr,dumnr,fnr,pdel,deltat,deltat_sed,nlev,i,MG_RAIN,g, &
         qrtend,nrtend,prect,qrsedten)
       gptl_ret = gptlstop_handle('Advance Solution (rain)',gptl_advance_sol_rain)
!       gptl_ret = gptlstop_handle('Sedimentation (rain)',gptl_sed_rain)
     end do

     write(fid_nstep,'(I12,4X)',advance="no") nstep

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     deltat_sed = deltat
     time_sed = deltat
     nstep = 0

     ! subcycle
     do while (time_sed > qsmall)
!       gptl_ret = gptlstart_handle('Sedimentation (snow)',gptl_sed_snow)
       ! obtain fall speeds
#ifndef SED_UPDATECFL
      if (nstep == 0) then
#endif
      gptl_ret = gptlstart_handle('CFL Calculation (snow)',gptl_cfl_snow)
      do n=1,gptl_loop
        call sed_CalcCFL(dums,dumns,precip_frac,rho,pdel,nlev,i,&
           MG_SNOW,deltat_sed,g,asn,rhof,alphaq,alphan,s1,s2,w1,w2,dum,&
           gamma_b_plus1=gamma_bs_plus1,gamma_b_plus4=gamma_bs_plus4)
      end do
      gptl_ret = gptlstop_handle('CFL Calculation (snow)',gptl_cfl_snow)
#ifndef SED_UPDATECFL
      else
        dum = CFL
      end if
#endif
       ! update deltat_sed for target CFL number
       deltat_sed = min(deltat_sed*CFL/dum,time_sed)
       ! advance snow sedimentation
       time_sed = time_sed - deltat_sed
       nstep = nstep + 1
       gptl_ret = gptlstart_handle('Flux Calculation (snow)',gptl_flux_snow)
       do n=1,gptl_loop
         call sed_CalcFlux(dums,dumns,alphaq,alphan,s1,s2,w1,w2,nlev,i,MG_SNOW,fs,fns)
       end do
       gptl_ret = gptlstop_handle('Flux Calculation (snow)',gptl_flux_snow)
       gptl_ret = gptlstart_handle('Advance Solution (snow)',gptl_advance_sol_snow)
       call sed_AdvanceOneStep(dums,fs,dumns,fns,pdel,deltat,deltat_sed,nlev,i,MG_SNOW,g, &
         qstend,nstend,prect,qssedten,preci=preci)
       gptl_ret = gptlstop_handle('Advance Solution (snow)',gptl_advance_sol_snow)
!       gptl_ret = gptlstop_handle('Sedimentation (snow)',gptl_sed_snow)
     end do

     write(fid_nstep,'(I12,4X)') nstep

     ! close sedimentation step file
     close(fid_nstep)

     ! end sedimentation
     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     ! get new update for variables that includes sedimentation tendency
     ! note : here dum variables are grid-average, NOT in-cloud

     do k=1,nlev

        dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)
        dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)
        dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)
        dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)

        dumr(i,k) = max(qr(i,k)+qrtend(i,k)*deltat,0._r8)
        dumnr(i,k) = max(nr(i,k)+nrtend(i,k)*deltat,0._r8)
        dums(i,k) = max(qs(i,k)+qstend(i,k)*deltat,0._r8)
        dumns(i,k) = max(ns(i,k)+nstend(i,k)*deltat,0._r8)

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)*lcldm(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)*icldm(i,k)
        end if

        if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
        if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
        if (dumr(i,k).lt.qsmall) dumnr(i,k)=0._r8
        if (dums(i,k).lt.qsmall) dumns(i,k)=0._r8

        ! calculate instantaneous processes (melting, homogeneous freezing)
        !====================================================================

        ! melting of snow at +2 C

        if (t(i,k)+tlat(i,k)/cpp*deltat > snowmelt) then
           if (dums(i,k) > 0._r8) then

              ! make sure melting snow doesn't reduce temperature below threshold
              dum = -xlf/cpp*dums(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt. snowmelt) then
                 dum = (t(i,k)+tlat(i,k)/cpp*deltat-snowmelt)*cpp/xlf
                 dum = dum/dums(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qstend(i,k)=qstend(i,k)-dum*dums(i,k)/deltat
              nstend(i,k)=nstend(i,k)-dum*dumns(i,k)/deltat
              qrtend(i,k)=qrtend(i,k)+dum*dums(i,k)/deltat
              nrtend(i,k)=nrtend(i,k)+dum*dumns(i,k)/deltat

              dum1=-xlf*dum*dums(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1
              meltsdttot(i,k)=meltsdttot(i,k) + dum1
           end if
        end if

        ! freezing of rain at -5 C

        if (t(i,k)+tlat(i,k)/cpp*deltat < rainfrze) then

           if (dumr(i,k) > 0._r8) then

              ! make sure freezing rain doesn't increase temperature above threshold
              dum = xlf/cpp*dumr(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.rainfrze) then
                 dum = -(t(i,k)+tlat(i,k)/cpp*deltat-rainfrze)*cpp/xlf
                 dum = dum/dumr(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qrtend(i,k)=qrtend(i,k)-dum*dumr(i,k)/deltat
              nrtend(i,k)=nrtend(i,k)-dum*dumnr(i,k)/deltat

              ! get mean size of rain = 1/lamr, add frozen rain to either snow or cloud ice
              ! depending on mean rain size

              call size_dist_param_basic(mg_rain_props, dumr(i,k), dumnr(i,k), &
                   lamr(i,k))

              if (lamr(i,k) < 1._r8/Dcs) then
                 qstend(i,k)=qstend(i,k)+dum*dumr(i,k)/deltat
                 nstend(i,k)=nstend(i,k)+dum*dumnr(i,k)/deltat
              else
                 qitend(i,k)=qitend(i,k)+dum*dumr(i,k)/deltat
                 nitend(i,k)=nitend(i,k)+dum*dumnr(i,k)/deltat
              end if

              ! heating tendency
              dum1 = xlf*dum*dumr(i,k)/deltat
              frzrdttot(i,k)=frzrdttot(i,k) + dum1
              tlat(i,k)=tlat(i,k)+dum1

           end if
        end if


        if (do_cldice) then
           if (t(i,k)+tlat(i,k)/cpp*deltat > tmelt) then
              if (dumi(i,k) > 0._r8) then

                 ! limit so that melting does not push temperature below freezing
                 !-----------------------------------------------------------------
                 dum = -dumi(i,k)*xlf/cpp
                 if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.tmelt) then
                    dum = (t(i,k)+tlat(i,k)/cpp*deltat-tmelt)*cpp/xlf
                    dum = dum/dumi(i,k)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qctend(i,k)=qctend(i,k)+dum*dumi(i,k)/deltat

                 ! for output
                 melttot(i,k)=dum*dumi(i,k)/deltat

                 ! assume melting ice produces droplet
                 ! mean volume radius of 8 micron

                 nctend(i,k)=nctend(i,k)+3._r8*dum*dumi(i,k)/deltat/ &
                      (4._r8*pi*5.12e-16_r8*rhow)

                 qitend(i,k)=((1._r8-dum)*dumi(i,k)-qi(i,k))/deltat
                 nitend(i,k)=((1._r8-dum)*dumni(i,k)-ni(i,k))/deltat
                 tlat(i,k)=tlat(i,k)-xlf*dum*dumi(i,k)/deltat
              end if
           end if

           ! homogeneously freeze droplets at -40 C
           !-----------------------------------------------------------------

           if (t(i,k)+tlat(i,k)/cpp*deltat < 233.15_r8) then
              if (dumc(i,k) > 0._r8) then

                 ! limit so that freezing does not push temperature above threshold
                 dum = dumc(i,k)*xlf/cpp
                 if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.233.15_r8) then
                    dum = -(t(i,k)+tlat(i,k)/cpp*deltat-233.15_r8)*cpp/xlf
                    dum = dum/dumc(i,k)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qitend(i,k)=qitend(i,k)+dum*dumc(i,k)/deltat
                 ! for output
                 homotot(i,k)=dum*dumc(i,k)/deltat

                 ! assume 25 micron mean volume radius of homogeneously frozen droplets
                 ! consistent with size of detrained ice in stratiform.F90
                 nitend(i,k)=nitend(i,k)+dum*3._r8*dumc(i,k)/(4._r8*3.14_r8*1.563e-14_r8* &
                      500._r8)/deltat
                 qctend(i,k)=((1._r8-dum)*dumc(i,k)-qc(i,k))/deltat
                 nctend(i,k)=((1._r8-dum)*dumnc(i,k)-nc(i,k))/deltat
                 tlat(i,k)=tlat(i,k)+xlf*dum*dumc(i,k)/deltat
              end if
           end if

           ! remove any excess over-saturation, which is possible due to non-linearity when adding
           ! together all microphysical processes
           !-----------------------------------------------------------------
           ! follow code similar to old CAM scheme

           qtmp=q(i,k)+qvlat(i,k)*deltat
           ttmp=t(i,k)+tlat(i,k)/cpp*deltat

           ! use rhw to allow ice supersaturation
           call qsat_water(ttmp, p(i,k), esn, qvn)

           if (qtmp > qvn .and. qvn > 0 .and. allow_sed_supersat) then
              ! expression below is approximate since there may be ice deposition
              dum = (qtmp-qvn)/(1._r8+xxlv_squared*qvn/(cpp*rv*ttmp**2))/deltat
              ! add to output cme
              cmeout(i,k) = cmeout(i,k)+dum
              ! now add to tendencies, partition between liquid and ice based on temperature
              if (ttmp > 268.15_r8) then
                 dum1=0.0_r8
                 ! now add to tendencies, partition between liquid and ice based on te
                 !-------------------------------------------------------
              else if (ttmp < 238.15_r8) then
                 dum1=1.0_r8
              else
                 dum1=(268.15_r8-ttmp)/30._r8
              end if

              dum = (qtmp-qvn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                   *qvn/(cpp*rv*ttmp**2))/deltat
              qctend(i,k)=qctend(i,k)+dum*(1._r8-dum1)
              ! for output
              qcrestot(i,k)=dum*(1._r8-dum1)
              qitend(i,k)=qitend(i,k)+dum*dum1
              qirestot(i,k)=dum*dum1
              qvlat(i,k)=qvlat(i,k)-dum
              ! for output
              qvres(i,k)=-dum
              tlat(i,k)=tlat(i,k)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
           end if
        end if

        ! calculate effective radius for pass to radiation code
        !=========================================================
        ! if no cloud water, default value is 10 micron for droplets,
        ! 25 micron for cloud ice

        ! update cloud variables after instantaneous processes to get effective radius
        ! variables are in-cloud to calculate size dist parameters

        dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)/lcldm(i,k)
        dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)/icldm(i,k)
        dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)/lcldm(i,k)
        dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)/icldm(i,k)

        dumr(i,k) = max(qr(i,k)+qrtend(i,k)*deltat,0._r8)/precip_frac(i,k)
        dumnr(i,k) = max(nr(i,k)+nrtend(i,k)*deltat,0._r8)/precip_frac(i,k)
        dums(i,k) = max(qs(i,k)+qstend(i,k)*deltat,0._r8)/precip_frac(i,k)
        dumns(i,k) = max(ns(i,k)+nstend(i,k)*deltat,0._r8)/precip_frac(i,k)

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)
        end if

        ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
        dumc(i,k)=min(dumc(i,k),5.e-3_r8)
        dumi(i,k)=min(dumi(i,k),5.e-3_r8)
        ! limit in-precip mixing ratios
        dumr(i,k)=min(dumr(i,k),10.e-3_r8)
        dums(i,k)=min(dums(i,k),10.e-3_r8)

        ! cloud ice effective radius
        !-----------------------------------------------------------------

        if (do_cldice) then
           if (dumi(i,k).ge.qsmall) then

              dum_2D(i,k) = dumni(i,k)
              call size_dist_param_basic(mg_ice_props, dumi(i,k), dumni(i,k), &
                   lami(i,k))

              if (dumni(i,k) /=dum_2D(i,k)) then
                 ! adjust number conc if needed to keep mean size in reasonable range
                 nitend(i,k)=(dumni(i,k)*icldm(i,k)-ni(i,k))/deltat
              end if

              effi(i,k) = 1.5_r8/lami(i,k)*1.e6_r8

           else
              effi(i,k) = 25._r8
           end if

           ! ice effective diameter for david mitchell's optics
           deffi(i,k)=effi(i,k)*rhoi/rhows*2._r8
        else
           ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
           ! radius has already been determined from the size distribution.
           effi(i,k) = re_ice(i,k) * 1.e6_r8      ! m -> um
           deffi(i,k)=effi(i,k) * 2._r8
        end if

        ! cloud droplet effective radius
        !-----------------------------------------------------------------
        if (dumc(i,k).ge.qsmall) then


           ! switch for specification of droplet and crystal number
           if (nccons) then
              ! make sure nc is consistence with the constant N by adjusting tendency, need
              ! to multiply by cloud fraction
              ! note that nctend may be further adjusted below if mean droplet size is
              ! out of bounds

              nctend(i,k)=(ncnst/rho(i,k)*lcldm(i,k)-nc(i,k))/deltat

           end if

           dum = dumnc(i,k)

           call size_dist_param_liq(mg_liq_props, dumc(i,k), dumnc(i,k), rho(i,k), &
                pgam(i,k), lamc(i,k))

           if (dum /= dumnc(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nctend(i,k)=(dumnc(i,k)*lcldm(i,k)-nc(i,k))/deltat
           end if

           effc(i,k) = (pgam(i,k)+3._r8)/lamc(i,k)/2._r8*1.e6_r8
           !assign output fields for shape here
           lamcrad(i,k)=lamc(i,k)
           pgamrad(i,k)=pgam(i,k)


           ! recalculate effective radius for constant number, in order to separate
           ! first and second indirect effects
           !======================================
           ! assume constant number of 10^8 kg-1

           dumnc(i,k)=1.e8_r8

           ! Pass in "false" adjust flag to prevent number from being changed within
           ! size distribution subroutine.
           call size_dist_param_liq(mg_liq_props, dumc(i,k), dumnc(i,k), rho(i,k), &
                pgam(i,k), lamc(i,k))

           effc_fn(i,k) = (pgam(i,k)+3._r8)/lamc(i,k)/2._r8*1.e6_r8

        else
           effc(i,k) = 10._r8
           lamcrad(i,k)=0._r8
           pgamrad(i,k)=0._r8
           effc_fn(i,k) = 10._r8
        end if

        ! recalculate 'final' rain size distribution parameters
        ! to ensure that rain size is in bounds, adjust rain number if needed

        if (dumr(i,k).ge.qsmall) then

           dum = dumnr(i,k)

           call size_dist_param_basic(mg_rain_props, dumr(i,k), dumnr(i,k), &
                lamr(i,k))

           if (dum /= dumnr(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nrtend(i,k)=(dumnr(i,k)*precip_frac(i,k)-nr(i,k))/deltat
           end if

        end if

        ! recalculate 'final' snow size distribution parameters
        ! to ensure that snow size is in bounds, adjust snow number if needed

        if (dums(i,k).ge.qsmall) then

           dum = dumns(i,k)

           call size_dist_param_basic(mg_snow_props, dums(i,k), dumns(i,k), &
                lams(i,k))

           if (dum /= dumns(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nstend(i,k)=(dumns(i,k)*precip_frac(i,k)-ns(i,k))/deltat
           end if

        end if


     end do ! vertical k loop

     do k=1,nlev
        ! if updated q (after microphysics) is zero, then ensure updated n is also zero
        !=================================================================================
        if (qc(i,k)+qctend(i,k)*deltat.lt.qsmall) nctend(i,k)=-nc(i,k)/deltat
        if (do_cldice .and. qi(i,k)+qitend(i,k)*deltat.lt.qsmall) nitend(i,k)=-ni(i,k)/deltat
        if (qr(i,k)+qrtend(i,k)*deltat.lt.qsmall) nrtend(i,k)=-nr(i,k)/deltat
        if (qs(i,k)+qstend(i,k)*deltat.lt.qsmall) nstend(i,k)=-ns(i,k)/deltat

     end do

  end do sed_col_loop! i loop

  ! DO STUFF FOR OUTPUT:
  !==================================================

  ! qc and qi are only used for output calculations past here,
  ! so add qctend and qitend back in one more time
  qc = qc + qctend*deltat
  qi = qi + qitend*deltat

  ! averaging for snow and rain number and diameter
  !--------------------------------------------------

  ! drout2/dsout2:
  ! diameter of rain and snow
  ! dsout:
  ! scaled diameter of snow (passed to radiation in CAM)
  ! reff_rain/reff_snow:
  ! calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual

  where (qrout .gt. 1.e-7_r8 &
       .and. nrout.gt.0._r8)
     qrout2 = qrout * precip_frac
     nrout2 = nrout * precip_frac
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just drout2 times constants.
     drout2 = avg_diameter(qrout, nrout, rho, rhow)
     freqr = precip_frac

     reff_rain=1.5_r8*drout2*1.e6_r8
  elsewhere
     qrout2 = 0._r8
     nrout2 = 0._r8
     drout2 = 0._r8
     freqr = 0._r8
     reff_rain = 0._r8
  end where

  where (qsout .gt. 1.e-7_r8 &
       .and. nsout.gt.0._r8)
     qsout2 = qsout * precip_frac
     nsout2 = nsout * precip_frac
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just dsout2 times constants.
     dsout2 = avg_diameter(qsout, nsout, rho, rhosn)
     freqs = precip_frac

     dsout=3._r8*rhosn/rhows*dsout2

     reff_snow=1.5_r8*dsout2*1.e6_r8
  elsewhere
     dsout  = 0._r8
     qsout2 = 0._r8
     nsout2 = 0._r8
     dsout2 = 0._r8
     freqs  = 0._r8
     reff_snow=0._r8
  end where

  ! analytic radar reflectivity
  !--------------------------------------------------
  ! formulas from Matthew Shupe, NOAA/CERES
  ! *****note: radar reflectivity is local (in-precip average)
  ! units of mm^6/m^3

  do i = 1,mgncol
     do k=1,nlev
        if (qc(i,k).ge.qsmall .and. (nc(i,k)+nctend(i,k)*deltat).gt.10._r8) then
           dum=(qc(i,k)/lcldm(i,k)*rho(i,k)*1000._r8)**2 &
                /(0.109_r8*(nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)/1.e6_r8)*lcldm(i,k)/precip_frac(i,k)
        else
           dum=0._r8
        end if
        if (qi(i,k).ge.qsmall) then
           dum1=(qi(i,k)*rho(i,k)/icldm(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*icldm(i,k)/precip_frac(i,k)
        else
           dum1=0._r8
        end if

        if (qsout(i,k).ge.qsmall) then
           dum1=dum1+(qsout(i,k)*rho(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)
        end if

        refl(i,k)=dum+dum1

        ! add rain rate, but for 37 GHz formulation instead of 94 GHz
        ! formula approximated from data of Matrasov (2007)
        ! rainrt is the rain rate in mm/hr
        ! reflectivity (dum) is in DBz

        if (rainrt(i,k).ge.0.001_r8) then
           dum=log10(rainrt(i,k)**6._r8)+16._r8

           ! convert from DBz to mm^6/m^3

           dum = 10._r8**(dum/10._r8)
        else
           ! don't include rain rate in R calculation for values less than 0.001 mm/hr
           dum=0._r8
        end if

        ! add to refl

        refl(i,k)=refl(i,k)+dum

        !output reflectivity in Z.
        areflz(i,k)=refl(i,k) * precip_frac(i,k)

        ! convert back to DBz

        if (refl(i,k).gt.minrefl) then
           refl(i,k)=10._r8*log10(refl(i,k))
        else
           refl(i,k)=-9999._r8
        end if

        !set averaging flag
        if (refl(i,k).gt.mindbz) then
           arefl(i,k)=refl(i,k) * precip_frac(i,k)
           frefl(i,k)=precip_frac(i,k)
        else
           arefl(i,k)=0._r8
           areflz(i,k)=0._r8
           frefl(i,k)=0._r8
        end if

        ! bound cloudsat reflectivity

        csrfl(i,k)=min(csmax,refl(i,k))

        !set averaging flag
        if (csrfl(i,k).gt.csmin) then
           acsrfl(i,k)=refl(i,k) * precip_frac(i,k)
           fcsrfl(i,k)=precip_frac(i,k)
        else
           acsrfl(i,k)=0._r8
           fcsrfl(i,k)=0._r8
        end if

     end do
  end do

  !redefine fice here....
  dum_2D = qsout + qrout + qc + qi
  dumi = qsout + qi
  where (dumi .gt. qsmall .and. dum_2D .gt. qsmall)
     nfice=min(dumi/dum_2D,1._r8)
  elsewhere
     nfice=0._r8
  end where


  ! Unpack all outputs

  ! Avoid zero/near-zero division.
  ! qcsinksum_rate1ord = qcsinksum_rate1ord/max(qcsum_rate1ord,1.0e-30_r8)

  ! call unpack_array(qcsinksum_rate1ord, mgcols, top_lev, 0._r8, rate1ord_cw2pr_st)

  call unpack_array(qcsinksum_rate1ord, mgcols, top_lev, 0._r8, qcsinksum_rate1ordo)

  call unpack_array(tlat, mgcols, top_lev, 0._r8, tlato)
  call unpack_array(qvlat, mgcols, top_lev, 0._r8, qvlato)

  call unpack_array(qctend, mgcols, top_lev, 0._r8, qctendo)
  call unpack_array(qitend, mgcols, top_lev, 0._r8, qitendo)

  ! hm mg2******
  call unpack_array(qrtend, mgcols, top_lev, 0._r8, qrtendo)
  call unpack_array(qstend, mgcols, top_lev, 0._r8, qstendo)

  ! Note that where there is no water, we set nctend and nitend to remove number
  ! concentration as well.
  call unpack_array(nctend, mgcols, top_lev, -ncn/deltat, nctendo)
  if (do_cldice) then
     call unpack_array(nitend, mgcols, top_lev, -nin/deltat, nitendo)
  else
     call unpack_array(nitend, mgcols, top_lev, 0._r8, nitendo)
  end if

  ! hm mg2 ****note we should apply a similar logic to unpacking of precip number
  ! as employed for cloud water and ice above
  call unpack_array(nrtend, mgcols, top_lev, -nrn/deltat, nrtendo)
  call unpack_array(nstend, mgcols, top_lev, -nsn/deltat, nstendo)

  call unpack_array(effc,    mgcols, top_lev, 10._r8, effco)
  call unpack_array(effc_fn, mgcols, top_lev, 10._r8, effco_fn)
  call unpack_array(effi,    mgcols, top_lev, 25._r8, effio)

  call unpack_array(prect, mgcols, 0._r8, precto)
  call unpack_array(preci, mgcols, 0._r8, precio)

  call unpack_array(nevapr,   mgcols, top_lev, 0._r8, nevapro)
  call unpack_array(evapsnow, mgcols, top_lev, 0._r8, evapsnowo)
  call unpack_array(prain,    mgcols, top_lev, 0._r8, praino)
  call unpack_array(prodsnow, mgcols, top_lev, 0._r8, prodsnowo)
  call unpack_array(cmeout,   mgcols, top_lev, 0._r8, cmeouto)

  call unpack_array(lamcrad, mgcols, top_lev, 0._r8, lamcrado)
  call unpack_array(pgamrad, mgcols, top_lev, 0._r8, pgamrado)
  call unpack_array(deffi,   mgcols, top_lev, 0._r8, deffio)

  call unpack_array(qsout,  mgcols, top_lev, 0._r8, qsouto)
  call unpack_array(qsout2, mgcols, top_lev, 0._r8, qsouto2)
  call unpack_array(nsout,  mgcols, top_lev, 0._r8, nsouto)
  call unpack_array(nsout2, mgcols, top_lev, 0._r8, nsouto2)
  call unpack_array(dsout,  mgcols, top_lev, 0._r8, dsouto)
  call unpack_array(dsout2, mgcols, top_lev, 0._r8, dsouto2)

  call unpack_array(qrout,  mgcols, top_lev, 0._r8, qrouto)
  call unpack_array(qrout2, mgcols, top_lev, 0._r8, qrouto2)
  call unpack_array(nrout,  mgcols, top_lev, 0._r8, nrouto)
  call unpack_array(nrout2, mgcols, top_lev, 0._r8, nrouto2)
  call unpack_array(drout2, mgcols, top_lev, 0._r8, drouto2)

  call unpack_array(reff_rain, mgcols, top_lev, 0._r8, reff_raino)
  call unpack_array(reff_snow, mgcols, top_lev, 0._r8, reff_snowo)

  call unpack_array(freqs, mgcols, top_lev, 0._r8, freqso)
  call unpack_array(freqr, mgcols, top_lev, 0._r8, freqro)

  call unpack_array(rflx, mgcols, top_lev, 0._r8, rflxo)
  call unpack_array(sflx, mgcols, top_lev, 0._r8, sflxo)

  call unpack_array(qcsevap, mgcols, top_lev, 0._r8, qcsevapo)
  call unpack_array(qisevap, mgcols, top_lev, 0._r8, qisevapo)
  call unpack_array(qvres,   mgcols, top_lev, 0._r8, qvreso)
  call unpack_array(cmeitot, mgcols, top_lev, 0._r8, cmeiout)
  call unpack_array(vtrmc,   mgcols, top_lev, 0._r8, vtrmco)
  call unpack_array(vtrmi,   mgcols, top_lev, 0._r8, vtrmio)
  !++ag
  call unpack_array(ums, mgcols, top_lev, 0._r8, umso)
  call unpack_array(umr, mgcols, top_lev, 0._r8, umro)
  !--ag
  call unpack_array(qcsedten, mgcols, top_lev, 0._r8, qcsedteno)
  call unpack_array(qisedten, mgcols, top_lev, 0._r8, qisedteno)
  call unpack_array(qrsedten, mgcols, top_lev, 0._r8, qrsedteno)
  call unpack_array(qssedten, mgcols, top_lev, 0._r8, qssedteno)

  call unpack_array(pratot,     mgcols,top_lev, 0._r8, prao)
  call unpack_array(prctot,     mgcols,top_lev, 0._r8, prco)
  call unpack_array(mnuccctot,  mgcols,top_lev, 0._r8, mnuccco)
  call unpack_array(mnuccttot,  mgcols,top_lev, 0._r8, mnuccto)
  call unpack_array(msacwitot,  mgcols,top_lev, 0._r8, msacwio)
  call unpack_array(psacwstot,  mgcols,top_lev, 0._r8, psacwso)
  call unpack_array(bergstot,   mgcols,top_lev, 0._r8, bergso)
  call unpack_array(bergtot,    mgcols,top_lev, 0._r8, bergo)
  call unpack_array(melttot,    mgcols,top_lev, 0._r8, melto)
  call unpack_array(homotot,    mgcols,top_lev, 0._r8, homoo)
  call unpack_array(qcrestot,   mgcols,top_lev, 0._r8, qcreso)
  call unpack_array(prcitot,    mgcols,top_lev, 0._r8, prcio)
  call unpack_array(praitot,    mgcols,top_lev, 0._r8, praio)
  call unpack_array(qirestot,   mgcols,top_lev, 0._r8, qireso)
  call unpack_array(mnuccrtot,  mgcols,top_lev, 0._r8, mnuccro)
  call unpack_array(pracstot,   mgcols,top_lev, 0._r8, pracso)
  call unpack_array(mnuccdtot,  mgcols,top_lev, 0._r8, mnuccdo)
  call unpack_array(meltsdttot, mgcols,top_lev, 0._r8, meltsdto)
  call unpack_array(frzrdttot,  mgcols,top_lev, 0._r8, frzrdto)

  call unpack_array(refl,   mgcols, top_lev, -9999._r8, reflo)
  call unpack_array(arefl,  mgcols, top_lev, 0._r8, areflo)
  call unpack_array(areflz, mgcols, top_lev, 0._r8, areflzo)
  call unpack_array(frefl,  mgcols, top_lev, 0._r8, freflo)
  call unpack_array(csrfl,  mgcols, top_lev, -9999._r8, csrflo)
  call unpack_array(acsrfl, mgcols, top_lev, 0._r8, acsrflo)
  call unpack_array(fcsrfl, mgcols, top_lev, 0._r8, fcsrflo)

  call unpack_array(rercld, mgcols, top_lev, 0._r8, rercldo)

  call unpack_array(nfice, mgcols, top_lev, 0._r8, nficeo)

  call unpack_array(qcrat, mgcols, top_lev, 1._r8, qcrato)

  call unpack_array(ncai, mgcols, top_lev, 0._r8, ncaio)
  call unpack_array(ncal, mgcols, top_lev, 0._r8, ncalo)

  call unpack_array(prer_evap, mgcols, top_lev, 1._r8, prer_evapo)

  ! update call count
  mphys_calls = mphys_calls + 1

end subroutine micro_mg2_acme_v1beta_tend

!========================================================================
!OUTPUT CALCULATIONS
!========================================================================

elemental subroutine calc_rercld(lamr, n0r, lamc, pgam, qric, qcic, ncic, &
     rercld)
  real(r8), intent(in) :: lamr          ! rain size parameter (slope)
  real(r8), intent(in) :: n0r           ! rain size parameter (intercept)
  real(r8), intent(in) :: lamc          ! size distribution parameter (slope)
  real(r8), intent(in) :: pgam          ! droplet size parameter
  real(r8), intent(in) :: qric          ! in-cloud rain mass mixing ratio
  real(r8), intent(in) :: qcic          ! in-cloud cloud liquid
  real(r8), intent(in) :: ncic          ! in-cloud droplet number concentration

  real(r8), intent(inout) :: rercld     ! effective radius calculation for rain + cloud

  ! combined size of precip & cloud drops
  real(r8) :: Atmp

  ! Rain drops
  if (lamr > 0._r8) then
     Atmp = n0r * pi / (2._r8 * lamr**3._r8)
  else
     Atmp = 0._r8
  end if

  ! Add cloud drops
  if (lamc > 0._r8) then
     Atmp = Atmp + &
          ncic * pi * rising_factorial(pgam+1._r8, 2)/(4._r8 * lamc**2._r8)
  end if

  if (Atmp > 0._r8) then
     rercld = rercld + 3._r8 *(qric + qcic) / (4._r8 * rhow * Atmp)
  end if

end subroutine calc_rercld

!========================================================================
!UTILITIES
!========================================================================

pure subroutine micro_mg2_acme_v1beta_get_cols(ncol, nlev, top_lev, qcn, qin, &
     qrn, qsn, mgncol, mgcols)

  ! Determines which columns microphysics should operate over by
  ! checking for non-zero cloud water/ice.

  integer, intent(in) :: ncol      ! Number of columns with meaningful data
  integer, intent(in) :: nlev      ! Number of levels to use
  integer, intent(in) :: top_lev   ! Top level for microphysics

  real(r8), intent(in) :: qcn(:,:) ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(:,:) ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: qrn(:,:) ! rain mixing ratio (kg/kg)
  real(r8), intent(in) :: qsn(:,:) ! snow mixing ratio (kg/kg)

  integer, intent(out) :: mgncol   ! Number of columns MG will use
  integer, allocatable, intent(out) :: mgcols(:) ! column indices

  integer :: lev_offset  ! top_lev - 1 (defined here for consistency)
  logical :: ltrue(ncol) ! store tests for each column

  integer :: i, ii ! column indices

  if (allocated(mgcols)) deallocate(mgcols)

  lev_offset = top_lev - 1

  ! Using "any" along dimension 2 collapses across levels, but
  ! not columns, so we know if water is present at any level
  ! in each column.

  ltrue = any(qcn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qin(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qrn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qsn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)

  ! Scan for true values to get a usable list of indices.

  mgncol = count(ltrue)
  allocate(mgcols(mgncol))
  i = 0
  do ii = 1,ncol
     if (ltrue(ii)) then
        i = i + 1
        mgcols(i) = ii
     end if
  end do

end subroutine micro_mg2_acme_v1beta_get_cols

! Subroutines to pack arrays into smaller, contiguous pieces
!========================================================================
! Rank 1 array of reals, columns only
pure subroutine pack_array_1Dr8(old_array, cols, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:)   ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include

  ! Output
  real(r8), intent(out) :: new_array(:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = old_array(cols)
  end if

end subroutine pack_array_1Dr8

! Rank 2 array of reals, columns and levels
pure subroutine pack_array_2Dr8(old_array, cols, top_lev, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:,:) ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include
  integer,  intent(in)  :: top_lev        ! First level to use

  ! Output
  real(r8), intent(out) :: new_array(:,:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = old_array(cols, top_lev:)
  end if

end subroutine pack_array_2Dr8

! Rank 3 array of reals, assume last index is extra
pure subroutine pack_array_3Dr8(old_array, cols, top_lev, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:,:,:) ! Array to be packed
  integer,  intent(in)  :: cols(:)          ! List of columns to include
  integer,  intent(in)  :: top_lev          ! First level to use

  ! Output
  real(r8), intent(out) :: new_array(:,:,:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = old_array(cols, top_lev:,:)
  end if

end subroutine pack_array_3Dr8

! Subroutines to unpack arrays for output
!========================================================================
! Rank 1 array of reals, columns only
pure subroutine unpack_array_1Dr8(old_array, cols, fill, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:)   ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include
  real(r8), intent(in)  :: fill           ! Value with which to fill unused
                                          ! sections of new_array.

  ! Output
  real(r8), intent(out) :: new_array(:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = fill

     new_array(cols) = old_array
  end if

end subroutine unpack_array_1Dr8

! Rank 1 array of reals, columns only, "fill" value is an array
pure subroutine unpack_array_1Dr8_arrayfill(old_array, cols, fill, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:)   ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include
  real(r8), intent(in)  :: fill(:)        ! Value with which to fill unused
                                          ! sections of new_array.

  ! Output
  real(r8), intent(out) :: new_array(:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = fill

     new_array(cols) = old_array
  end if

end subroutine unpack_array_1Dr8_arrayfill

! Rank 2 array of reals, columns and levels
pure subroutine unpack_array_2Dr8(old_array, cols, top_lev, fill, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:,:) ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include
  integer,  intent(in)  :: top_lev        ! First level to use
  real(r8), intent(in)  :: fill           ! Value with which to fill unused
                                          ! sections of new_array.

  ! Output
  real(r8), intent(out) :: new_array(:,:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = fill

     new_array(cols, top_lev:) = old_array
  end if

end subroutine unpack_array_2Dr8

! Rank 2 array of reals, columns and levels, "fill" value is an array
pure subroutine unpack_array_2Dr8_arrayfill(old_array, cols, top_lev, fill, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:,:) ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include
  integer,  intent(in)  :: top_lev        ! First level to use
  real(r8), intent(in)  :: fill(:,:)      ! Value with which to fill unused
                                          ! sections of new_array.

  ! Output
  real(r8), intent(out) :: new_array(:,:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = fill

     new_array(cols, top_lev:) = old_array
  end if

end subroutine unpack_array_2Dr8_arrayfill



!!== KZ_DCS
subroutine get_dcst(ncol,pver,temp,dcst)

implicit none

integer,  intent(in) :: ncol
integer,  intent(in) :: pver                 ! number of layers in columns
real(r8), intent(in) :: temp(ncol,pver)       ! input temperature (K)
real(r8), intent(out) :: dcst(ncol,pver)      ! temperature dependent dcs

integer :: i,k
real(r8) :: st


dcst = 400.e-6_r8

do k=1,pver
   do i=1,ncol
      st = temp(i,k) - 273.15
      if(st.le.-70.) then
         dcst(i,k) = 100.e-6_r8
      elseif(st.gt.-70. .and. st.le.-10.) then
         dcst(i,k) = 5.e-6_r8 * st  + 450.e-6_r8
      elseif(st.gt.-10.) then
         dcst(i,k) = 400.e-6_r8
      end if
   end do
end do

return

end subroutine get_dcst
!!== KZ_DCS

end module module_mp_mg2_acme_v1beta
