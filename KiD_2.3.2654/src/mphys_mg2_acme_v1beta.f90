! *****************************COPYRIGHT*******************************
! based on mphys_morr_two_moment
! *****************************COPYRIGHT*******************************
!
! Module containing interface to Morrison & Gettelman version 2 (mg2)
!  Andrew Gettelman, February 2013
!  Updated May 2013 for MG2 code in MG2_dev_n07_CAM5_3_01

#define ADJUST_SATURATION_BEFORE
!#define ADJUST_SATURATION_AFTER

module mphys_mg2_acme_v1beta

  use parameters, only : num_h_moments, num_h_bins, nspecies, nz, dt, &
       h_names, mom_units, max_char_len, mom_names, nx
  use column_variables
  use physconst, only : p0, r_on_cp, pi

  use diagnostics, only: save_dg, i_dgtime
  use common_physics, only : qsaturation, qisaturation
  
  use wv_sat_methods, only: wv_sat_qsat_water, wv_sat_qsat_ice, wv_sat_methods_init

  use module_mp_mg2_acme_v1beta, only: micro_mg2_acme_v1beta_init, &
       micro_mg2_acme_v1beta_tend

  use mphys_stats
  
  implicit none
  
  ! Logical switches 
  logical :: micro_unset=.True.
  integer:: ih, imom
  character(max_char_len) :: name, units

contains

  subroutine mphys_mg2_acme_v1beta_interface

    real :: t1d(nz), p1d(nz), dz1d(nz),qv1d(nz),qc1d(nz) &
         , qi1d(nz), ni1d(nz) & 
         , qr1d(nz), qs1d(nz) &
         , ns1d(nz), nr1d(nz) &
         , w1d(nz)

    real :: wvar1d(nz)

    real :: qc_tend1d(nz), ni_tend1d(nz)  &
         ,  qi_tend1d(nz), qni_tend1d(nz) & 
         ,  qv_tend1d(nz), t_tend1d(nz)   &
         ,  qr_tend1d(nz), nr_tend1d(nz)  &
         ,  qs_tend1d(nz), ns_tend1d(nz) 
    !         qg_tend1d(nz), ng_tend1d(nz)

    real :: precprt1d, snowrt1d

    real :: effc1d(nz), effi1d(nz), effs1d(nz),        &
         effr1d(nz), effg1d(nz)

    real :: qrcu1d(nz), qscu1d(nz), qicu1d(nz)

    real :: qgsten(nz), qrsten(nz), qisten(nz),  &
         qnisten(nz), qcsten(nz)

    ! KiD_2D diag arrays
    real :: qrsten_2d(nz,nx),precprt2d(nx), snowrt2d(nx) 

    integer :: kts, kte, i, j, k

    !++ag

    ! Init Variables
    character(128) :: errstring
    integer, parameter :: kind = selected_real_kind(12) ! 8 byte real
    integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
    real(r8) :: gravit
    real(r8) :: rair
    real(r8) :: rh2o
    real(r8) :: cpair
    real(r8) :: tmelt_in
    real(r8) :: latvap
    real(r8) :: latice
    real(r8) :: rhmini_in
    real(r8) :: micro_mg_dcs
    logical  :: micro_mg_dcs_tdep
    logical  :: microp_uniform_in
    logical  :: do_cldice_in
    logical  :: use_hetfrz_classnuc_in
    character(len=16) :: micro_mg_precip_frac_method_in  ! type of precip
    real(r8) :: micro_mg_berg_eff_factor_in
    logical  :: allow_sed_supersat_in
    real(r8) :: ice_sed_ai
    real(r8) :: prc_coef1_in
    real(r8) :: prc_exp_in
    real(r8) :: prc_exp1_in
    real(r8) :: cld_sed_in
    logical  :: mg_prc_coeff_fix_in

    real(r8) :: h2otrip_in
    real(r8) :: tboil_in
    real(r8) :: ttrice_in
    real(r8) :: epsilo_in 

    ! From shr_const_mod.F90
    real(R8),parameter :: SHR_CONST_G       = 9.80616_R8      ! acceleration of gravity ~ m/s^2
    real(R8),parameter :: SHR_CONST_BOLTZ   = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
    real(R8),parameter :: SHR_CONST_AVOGAD  = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole
    real(R8),parameter :: SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
    real(R8),parameter :: SHR_CONST_MWDAIR  = 28.966_R8       ! molecular weight dry air ~ kg/kmole
    real(R8),parameter :: SHR_CONST_MWWV    = 18.016_R8       ! molecular weight water vapor
    real(R8),parameter :: SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant     ~ J/K/kg
    real(R8),parameter :: SHR_CONST_RWV     = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg
    real(R8),parameter :: SHR_CONST_TKFRZ   = 273.15_R8       ! freezing T of fresh water          ~ K 
    real(R8),parameter :: SHR_CONST_LATICE  = 3.337e5_R8      ! latent heat of fusion      ~ J/kg
    real(R8),parameter :: SHR_CONST_LATVAP  = 2.501e6_R8      ! latent heat of evaporation ~ J/kg
    real(R8),parameter :: SHR_CONST_CPDAIR  = 1.00464e3_R8    ! specific heat of dry air   ~ J/kg/K

    real(R8),parameter :: qsmall = 1.e-12_R8

    ! DUMMY VARIABLES

    REAL(R8) DUMT,DUMQV,DUMQSS,DUMQC,PCC,DUMS, DNC,DUMNC,ESL,QVL,CPM

    ! Inputs and outputs.

    ! -----------------------------------------------------------------------------
    ! input arguments
    ! -----------------------------------------------------------------------------
    integer :: mgncol            ! number of microphysics columns
    integer :: mgcols(nx)        ! list of microphysics columns
    integer :: nlev              ! number of layers
    integer :: top_lev           ! top level to do microphysics
    real(r8) :: deltatin         ! time step (s)
    real(r8) :: tn(nx,nz)        ! input temperature (K)
    real(r8) :: qn(nx,nz)        ! input h20 vapor mixing ratio (kg/kg)

    ! note: all input cloud variables are grid-averaged
    real(r8) :: qcn(nx,nz)       ! cloud water mixing ratio (kg/kg)
    real(r8) :: qin(nx,nz)       ! cloud ice mixing ratio (kg/kg)
    real(r8) :: ncn(nx,nz)       ! cloud water number conc (1/kg)
    real(r8) :: nin(nx,nz)       ! cloud ice number conc (1/kg)
    real(r8) :: qrn(nx,nz)       ! rain water mixing ratio (kg/kg)
    real(r8) :: qsn(nx,nz)       ! snow ice mixing ratio (kg/kg)
    real(r8) :: nrn(nx,nz)       ! rain water number conc (1/kg)
    real(r8) :: nsn(nx,nz)       ! snow ice number conc (1/kg)
    real(r8) :: pn(nx,nz)        ! air pressure (pa)
    real(r8) :: pdeln(nx,nz)     ! pressure difference across level (pa)
    ! hm add 11-16-11, interface pressure
    real(r8) :: pint(nx,0:nz)    ! level interface pressure (pa)
    real(r8) :: cldn(nx,nz)      ! cloud fraction (no units)
    real(r8) :: liqcldf(nx,nz)   ! liquid cloud fraction (no units)
    real(r8) :: icecldf(nx,nz)   ! ice cloud fraction (no units)
    ! used for scavenging
    ! Inputs for aerosol activation
    real(r8) :: naain(nx,nz)     ! ice nucleation number (from microp_aero_ts) (1/kg)
    real(r8) :: npccnin(nx,nz)   ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)

    ! Note that for these variables, the dust bin is assumed to be the last index.
    ! (For example, in CAM, the last dimension is always size 4.)
    real(r8) :: rndstn(nx,nz,4)  ! radius of each dust bin, for contact freezing (from microp_aero_ts) (m)
    real(r8) :: naconin(nx,nz,4) ! number in each dust bin, for contact freezing  (from microp_aero_ts) (1/m^3)

    ! Used with CARMA cirrus microphysics
    ! (or similar external microphysics model)
    real(r8) :: tnd_qsnown(nx,nz) ! snow mass tendency (kg/kg/s)
    real(r8) :: tnd_nsnown(nx,nz) ! snow number tendency (#/kg/s)
    real(r8) :: re_icen(nx,nz)    ! ice effective radius (m)

    real(r8) :: frzimmn(nx,nz) ! Number tendency due to immersion freezing (1/cm3)
    real(r8) :: frzcntn(nx,nz) ! Number tendency due to contact freezing (1/cm3)
    real(r8) :: frzdepn(nx,nz) ! Number tendency due to deposition nucleation (1/cm3)

    ! -----------------------------------------------------------------------------
    ! output arguments
    ! -----------------------------------------------------------------------------
    real(r8) :: qcsinksum_rate1ordo(nx,nz)    ! 1st order rate for direct cw to precip conversion
    ! real(r8) :: rate1ord_cw2pr_st(nx,nz)    ! 1st order rate for direct cw to precip conversion
    real(r8) :: tlato(nx,nz)         ! latent heating rate       (W/kg)
    real(r8) :: qvlato(nx,nz)        ! microphysical tendency qv (1/s)
    real(r8) :: qctendo(nx,nz)       ! microphysical tendency qc (1/s)
    real(r8) :: qitendo(nx,nz)       ! microphysical tendency qi (1/s)
    real(r8) :: nctendo(nx,nz)       ! microphysical tendency nc (1/(kg*s))
    real(r8) :: nitendo(nx,nz)       ! microphysical tendency ni (1/(kg*s))

    real(r8) :: qrtendo(nx,nz)       ! microphysical tendency qr (1/s)
    real(r8) :: qstendo(nx,nz)       ! microphysical tendency qs (1/s)
    real(r8) :: nrtendo(nx,nz)       ! microphysical tendency nr (1/(kg*s))
    real(r8) :: nstendo(nx,nz)       ! microphysical tendency ns (1/(kg*s))

    real(r8) :: effco(nx,nz)         ! droplet effective radius (micron)
    real(r8) :: effco_fn(nx,nz)      ! droplet effective radius, assuming nc = 1.e8 kg-1
    real(r8) :: effio(nx,nz)         ! cloud ice effective radius (micron)
    real(r8) :: precto(nx)           ! surface precip rate (m/s)
    real(r8) :: precio(nx)           ! cloud ice/snow precip rate (m/s)
    real(r8) :: nevapro(nx,nz)       ! evaporation rate of rain + snow (1/s)
    real(r8) :: evapsnowo(nx,nz)     ! sublimation rate of snow (1/s)
    real(r8) :: praino(nx,nz)        ! production of rain + snow (1/s)
    real(r8) :: prodsnowo(nx,nz)     ! production of snow (1/s)
    real(r8) :: cmeouto(nx,nz)       ! evap/sub of cloud (1/s)
    real(r8) :: deffio(nx,nz)        ! ice effective diameter for optics (radiation) (micron)
    real(r8) :: pgamrado(nx,nz)      ! ice gamma parameter for optics (radiation) (no units)
    real(r8) :: lamcrado(nx,nz)      ! slope of droplet distribution for optics (radiation) (1/m)
    real(r8) :: qsouto(nx,nz)        ! snow mixing ratio (kg/kg)
    real(r8) :: dsouto(nx,nz)        ! snow diameter (m)
    real(r8) :: rflxo(nx,nz+1)       ! grid-box average rain flux (kg m^-2 s^-1) PMC added +1 for consistency w/ module
    real(r8) :: sflxo(nx,nz+1)       ! grid-box average snow flux (kg m^-2 s^-1) PMC added +1 for consistency w/ module
    real(r8) :: qrouto(nx,nz)        ! grid-box average rain mixing ratio (kg/kg)
    real(r8) :: reff_raino(nx,nz)    ! rain effective radius (micron)
    real(r8) :: reff_snowo(nx,nz)    ! snow effective radius (micron)
    real(r8) :: qcsevapo(nx,nz)      ! cloud water evaporation due to sedimentation (1/s)
    real(r8) :: qisevapo(nx,nz)      ! cloud ice sublimation due to sublimation (1/s)
    real(r8) :: qvreso(nx,nz)        ! residual condensation term to ensure RH < 100% (1/s)
    real(r8) :: cmeiout(nx,nz)       ! grid-mean cloud ice sub/dep (1/s)
    real(r8) :: vtrmco(nx,nz)        ! mass-weighted cloud water fallspeed (m/s)
    real(r8) :: vtrmio(nx,nz)        ! mass-weighted cloud ice fallspeed (m/s)
    !++ag
    real(r8) :: umso(nx,nz)   ! mass weighted snow fallspeed (m/s) 
    real(r8) :: umro(nx,nz)   ! mass weighted rain fallspeed (m/s)
    real(r8) :: satadj(nx,nz) ! saturation adjustment tendency (kg kg-1 s-1)
    !--ag

    real(r8) :: qcsedteno(nx,nz)    ! qc sedimentation tendency (1/s)
    real(r8) :: qisedteno(nx,nz)    ! qi sedimentation tendency (1/s)
    real(r8) :: qrsedteno(nx,nz)    ! qc sedimentation tendency (1/s)
    real(r8) :: qssedteno(nx,nz)    ! qi sedimentation tendency (1/s)

    ! microphysical process rates for output (mixing ratio tendencies) (all have units of 1/s)
    real(r8) :: prao(nx,nz)         ! accretion of cloud by rain 
    real(r8) :: prco(nx,nz)         ! autoconversion of cloud to rain
    real(r8) :: mnuccco(nx,nz)      ! mixing ratio tend due to immersion freezing
    real(r8) :: mnuccto(nx,nz)      ! mixing ratio tend due to contact freezing
    real(r8) :: msacwio(nx,nz)      ! mixing ratio tend due to H-M splintering
    real(r8) :: psacwso(nx,nz)      ! collection of cloud water by snow
    real(r8) :: bergso(nx,nz)       ! bergeron process on snow
    real(r8) :: bergo(nx,nz)        ! bergeron process on cloud ice
    real(r8) :: melto(nx,nz)        ! melting of cloud ice
    real(r8) :: homoo(nx,nz)        ! homogeneous freezing cloud water
    real(r8) :: qcreso(nx,nz)       ! residual cloud condensation due to removal of excess supersat
    real(r8) :: prcio(nx,nz)        ! autoconversion of cloud ice to snow
    real(r8) :: praio(nx,nz)        ! accretion of cloud ice by snow
    real(r8) :: qireso(nx,nz)       ! residual ice deposition due to removal of excess supersat
    real(r8) :: mnuccro(nx,nz)      ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
    real(r8) :: pracso(nx,nz)       ! mixing ratio tendency due to accretion of rain by snow (1/s)
    real(r8) :: meltsdto(nx,nz)     ! latent heating rate due to melting of snow  (W/kg)
    real(r8) :: frzrdto(nx,nz)      ! latent heating rate due to homogeneous freezing of rain (W/kg)
    real(r8) :: mnuccdo(nx,nz)      ! mass tendency from ice nucleation
    real(r8) :: nrouto(nx,nz)       ! rain number concentration (1/m3)
    real(r8) :: nsouto(nx,nz)       ! snow number concentration (1/m3)
    real(r8) :: reflo(nx,nz)        ! analytic radar reflectivity        
    real(r8) :: areflo(nx,nz)       ! average reflectivity will zero points outside valid range
    real(r8) :: areflzo(nx,nz)      ! average reflectivity in z.
    real(r8) :: freflo(nx,nz)       ! fractional occurrence of radar reflectivity
    real(r8) :: csrflo(nx,nz)       ! cloudsat reflectivity 
    real(r8) :: acsrflo(nx,nz)      ! cloudsat average
    real(r8) :: fcsrflo(nx,nz)      ! cloudsat fractional occurrence of radar reflectivity
    real(r8) :: rercldo(nx,nz)      ! effective radius calculation for rain + cloud
    real(r8) :: ncaio(nx,nz)        ! output number conc of ice nuclei available (1/m3)
    real(r8) :: ncalo(nx,nz)        ! output number conc of CCN (1/m3)
    real(r8) :: qrouto2(nx,nz)      ! copy of qrout as used to compute drout2
    real(r8) :: qsouto2(nx,nz)      ! copy of qsout as used to compute dsout2
    real(r8) :: nrouto2(nx,nz)      ! copy of nrout as used to compute drout2
    real(r8) :: nsouto2(nx,nz)      ! copy of nsout as used to compute dsout2
    real(r8) :: drouto2(nx,nz)      ! mean rain particle diameter (m)
    real(r8) :: dsouto2(nx,nz)      ! mean snow particle diameter (m)
    real(r8) :: freqso(nx,nz)       ! fractional occurrence of snow
    real(r8) :: freqro(nx,nz)       ! fractional occurrence of rain
    real(r8) :: nficeo(nx,nz)       ! fractional occurrence of ice
    real(r8) :: qcrato(nx,nz)       ! limiter for qc process rates (1=no limit --> 0. no qc)

    real(r8) :: prer_evapo(nx,nz)

    ! -----------------------------------------------------------------------------
    ! Flipped arrays Input:  tn,qn,qcn,ncn,qin,nin,cldn,liqcldf,icecldf
    ! -----------------------------------------------------------------------------
    real(r8) ::flip_pn(nx,nz)
    real(r8) ::flip_pdeln(nx,nz)
    !  real(r8) ::flip_pint(nx,0:nz)
    real(r8) ::flip_tn(nx,nz)
    real(r8) ::flip_qn(nx,nz)
    real(r8) ::flip_qcn(nx,nz)
    real(r8) ::flip_ncn(nx,nz)
    real(r8) ::flip_qin(nx,nz)
    real(r8) ::flip_nin(nx,nz) 
    real(r8) ::flip_qrn(nx,nz)
    real(r8) ::flip_nrn(nx,nz)
    real(r8) ::flip_qsn(nx,nz)
    real(r8) ::flip_nsn(nx,nz) 

    real(r8) ::flip_relvarn(nx,nz)
    real(r8) ::flip_accre_enhann(nx,nz) 

    real(r8) ::flip_cldn(nx,nz)
    real(r8) ::flip_liqcldf(nx,nz)
    real(r8) ::flip_icecldf(nx,nz)  
    real(r8) ::flip_naain(nx,nz) 
    real(r8) ::flip_npccnin(nx,nz)

    ! -----------------------------------------------------------------------------
    ! Flipped arrays output: tlato,qvlato,qctendo,qitendo,nctendo,nitendo
    ! -----------------------------------------------------------------------------

    !tendencies
    real(r8) ::flip_tlato(nx,nz)
    real(r8) ::flip_qvlato(nx,nz)
    real(r8) ::flip_qctendo(nx,nz)
    real(r8) ::flip_qitendo(nx,nz)
    real(r8) ::flip_nitendo(nx,nz)
    real(r8) ::flip_nctendo(nx,nz)
    real(r8) ::flip_qrtendo(nx,nz)
    real(r8) ::flip_qstendo(nx,nz)
    real(r8) ::flip_nrtendo(nx,nz)
    real(r8) ::flip_nstendo(nx,nz)
    !other (for output)
    real(r8) ::flip_qrouto(nx,nz)
    real(r8) ::flip_nrouto(nx,nz) 
    real(r8) ::flip_qrouto2(nx,nz)
    real(r8) ::flip_nrouto2(nx,nz)  
    real(r8) ::flip_prao(nx,nz)
    real(r8) ::flip_prco(nx,nz)
    real(r8) ::flip_vtrmco(nx,nz)
    real(r8) ::flip_vtrmio(nx,nz)
    real(r8) ::flip_umso(nx,nz)
    real(r8) ::flip_umro(nx,nz)
    real(r8) ::flip_cmeiout(nx,nz)

    real(r8) :: flip_qcsedteno(nx,nz)    ! qc sedimentation tendency (1/s)
    real(r8) :: flip_qisedteno(nx,nz)    ! qi sedimentation tendency (1/s)
    real(r8) :: flip_qrsedteno(nx,nz)    ! qc sedimentation tendency (1/s)
    real(r8) :: flip_qssedteno(nx,nz)    ! qi sedimentation tendency (1/s)

    real(r8) :: flip_mnuccco(nx,nz)      ! mixing ratio tend due to immersion freezing
    real(r8) :: flip_mnuccto(nx,nz)      ! mixing ratio tend due to contact freezing
    real(r8) :: flip_msacwio(nx,nz)      ! mixing ratio tend due to H-M splintering
    real(r8) :: flip_psacwso(nx,nz)      ! collection of cloud water by snow
    real(r8) :: flip_bergso(nx,nz)       ! bergeron process on snow
    real(r8) :: flip_bergo(nx,nz)        ! bergeron process on cloud ice
    real(r8) :: flip_melto(nx,nz)        ! melting of cloud ice
    real(r8) :: flip_homoo(nx,nz)        ! homogeneous freezing cloud water
    real(r8) :: flip_qcreso(nx,nz)       ! residual cloud condensation due to removal of excess supersat
    real(r8) :: flip_prcio(nx,nz)        ! autoconversion of cloud ice to snow
    real(r8) :: flip_praio(nx,nz)        ! accretion of cloud ice by snow
    real(r8) :: flip_qireso(nx,nz)       ! residual ice deposition due to removal of excess supersat
    real(r8) :: flip_mnuccro(nx,nz)      ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
    real(r8) :: flip_pracso(nx,nz)       ! mixing ratio tendency due to accretion of rain by snow (1/s)
    real(r8) :: flip_meltsdto(nx,nz)     ! latent heating rate due to melting of snow  (W/kg)
    real(r8) :: flip_frzrdto(nx,nz)      ! latent heating rate due to homogeneous freezing of rain (W/kg)
    real(r8) :: flip_mnuccdo(nx,nz)      ! mass tendency from ice nucleation
    !--ag

    ! -----------------------------------------------------------------------------
    ! Generic Flipped arrays 
    ! -----------------------------------------------------------------------------
    real(r8) :: flip_data(nx,nz) ! temporary fliped data array

#ifdef ADJUST_SATURATION_BEFORE
    ! internal variables for saturation adjustment
    real(r8) :: qvlato_adj(nx,nz)  
    real(r8) :: tlato_adj(nx,nz)   
    real(r8) :: qctendo_adj(nx,nz) 
    real(r8) :: nctendo_adj(nx,nz) 
    real(r8) :: qitendo_adj(nx,nz) 
    real(r8) :: nitendo_adj(nx,nz) 
#endif
    
    !note different order of i,k for MG...
    !also: MG is top down, kid arrays are bottom up.    

    kts = 1
    kte = nz
    j   = 1

    qrouto(:,:) = 0.0_wp
    nrouto(:,:) = 0.0_wp

    flip_qrouto2(:,:) = 0.0_wp
    flip_nrouto2(:,:) = 0.0_wp

    do i=1,nx
       !zero out precipitation
       precto(i) = 0.0_wp
       precio(i) = 0.0_wp

       !interface pressures
       pint(i,:) = 100.0_wp * pmb_half(:,i) 

       do k=1,nz
          ! zero some of these for safety
          tlato(i,k)   = 0.0_wp
          qvlato(i,k)  = 0.0_wp
          qctendo(i,k) = 0.0_wp
          qitendo(i,k) = 0.0_wp
          nctendo(i,k) = 0.0_wp
          nitendo(i,k) = 0.0_wp
          qrtendo(i,k) = 0.0_wp
          qstendo(i,k) = 0.0_wp
          nrtendo(i,k) = 0.0_wp
          nstendo(i,k) = 0.0_wp

          effco(i,k)    = 0.0_wp
          effio(i,k)    = 0.0_wp
          effco_fn(i,k) = 0.0_wp

          !set input arguments

          !temperature
          ! tn(i,k) = (theta(k,i) + (dtheta_adv(k,i)+dtheta_div(k,i))*dt )*exner(k,i)
          tn(i,k) = theta(k,i) * exner(k,i)

          !layer pressure
          pn(i,k) = p0*exner(k,i)**(1.0_wp/r_on_cp)

          !interface pressure calculated using dp = rho * g* dz
          pdeln(i,k)= pn(i,k)/(287.15*tn(i,k)) * 9.81_wp * dz(k)

          ! sequential updates with advection and divergence tendencies (note dt = dtm)

          ! water vapor mass
          ! qn(i,k) = qv(k,i)+ (dqv_adv(k,i)+dqv_div(k,i))*dt
          qn(i,k) = qv(k,i)

          if (qn(i,k).lt.qsmall) qn(i,k) = 0.0_wp

          ! cloud water mass
          ! qcn(i,k) = hydrometeors(k,i,1)%moments(1,1) & 
          !      + (dhydrometeors_adv(k,i,1)%moments(1,1) &
          !      + dhydrometeors_div(k,i,1)%moments(11,))*dt
          qcn(i,k) = hydrometeors(k,i,1)%moments(1,1)

          if (qcn(i,k).lt.qsmall) qcn(i,k) = 0.0_wp

          ! cloud water number
          if (num_h_moments(1) >= 2) &
               ! ncn(i,k) = hydrometeors(k,i,1)%moments(1,2)& 
               ! + (dhydrometeors_adv(k,i,1)%moments(1,2) &
               ! + dhydrometeors_div(k,i,1)%moments(1,2))*dt
               ncn(i,k) = hydrometeors(k,i,1)%moments(1,2)
          
          if (qcn(i,k).lt.qsmall) ncn(i,k) = 0.0_wp

          ! cloud ice mass
          if (num_h_moments(3) >= 1) &
               ! qin(i,k) = hydrometeors(k,i,3)%moments(1,1)& 
               ! + (dhydrometeors_adv(k,i,3)%moments(1,1) &
               ! + dhydrometeors_div(k,i,3)%moments(1,1))*dt
               qin(i,k) = hydrometeors(k,i,3)%moments(1,1)

          if (qin(i,k).lt.qsmall) qin(i,k)=0.0_wp

          ! cloud ice number
          if (num_h_moments(3) >= 2) &
               ! nin(i,k) = hydrometeors(k,i,3)%moments(1,2)& 
               ! + (dhydrometeors_adv(k,i,3)%moments(1,2) &
               ! + dhydrometeors_div(k,i,3)%moments(1,2))*dt
               nin(i,k) = hydrometeors(k,i,3)%moments(1,2)

          if (qin(i,k).lt.qsmall) nin(i,k) = 0.0_wp

          ! rain mass
          if (num_h_moments(2) >= 1) &
               ! qrn(i,k) = hydrometeors(k,i,2)%moments(1,1)& 
               ! + (dhydrometeors_adv(k,i,2)%moments(1,1) &
               ! + dhydrometeors_div(k,i,2)%moments(1,1))*dt
               qrn(i,k) = hydrometeors(k,i,2)%moments(1,1)

          if (qrn(i,k).lt.qsmall) qrn(i,k)=0.0_wp 

          ! rain number
          if (num_h_moments(2) >= 2) &
               ! nrn(i,k) = hydrometeors(k,i,2)%moments(1,2)& 
               ! + (dhydrometeors_adv(k,i,2)%moments(1,2) &
               ! + dhydrometeors_div(k,i,2)%moments(1,2))*dt
               nrn(i,k) = hydrometeors(k,i,2)%moments(1,2)

          if (qrn(i,k).lt.qsmall) nrn(i,k)=0.0_wp

          ! snow mass
          if (num_h_moments(4) >= 1) &
               ! qsn(i,k) = hydrometeors(k,i,4)%moments(1,1)& 
               ! + (dhydrometeors_adv(k,i,4)%moments(1,1) &
               ! + dhydrometeors_div(k,i,4)%moments(1,1))*dt
               qsn(i,k) = hydrometeors(k,i,4)%moments(1,1)

          if (qsn(i,k).lt.qsmall) qsn(i,k)=0.0_wp

          ! snow number
          if (num_h_moments(4) >= 2) &
               ! nsn(i,k) = hydrometeors(k,i,4)%moments(1,2)& 
               ! + (dhydrometeors_adv(k,i,4)%moments(1,2) &
               ! + dhydrometeors_div(k,i,4)%moments(1,2))*dt
               nsn(i,k) = hydrometeors(k,i,4)%moments(1,2)

          if (qsn(i,k).lt.qsmall) nsn(i,k)=0.0_wp

          !set cloud information
          liqcldf(i,k)= 0.0_wp
          icecldf(i,k)= 0.0_wp
          cldn(i,k)=0.0_wp

          if (qcn(i,k) > qsmall) liqcldf(i,k) = 1.0_r8
          if (qin(i,k) > qsmall) icecldf(i,k) = 1.0_r8
          if (qcn(i,k) > qsmall .or. qin(i,k) > qsmall) cldn(i,k) = 1.0_r8

          !set activated aerosol number (uniform for now)
          !          if (tn(i,k) > tmelt_in) &
          npccnin = 10.0e6_wp   ! #/m3   
          !          if (tn(i,k) < tmelt_in - 10.) &
          naain = 0.0e6_wp   ! #/m3             

       end do
    end do


    !print*,'pn(1),pn(nz),pint(0),pint(nz),pdel(1),pdel(nz)=',pn(1,1),pn(1,nz),pint(1,0),pint(1,nz),pdeln(1,1),pdeln(1,nz)

    !set up constants for microphysics here...
    gravit   = SHR_CONST_G      ! Gravity
    rair     = SHR_CONST_RDAIR  ! dry air gas constant: note units (phys_constants are in J/K/kmol)
    rh2o     = SHR_CONST_RWV    ! water vapor gas constant
    cpair    = SHR_CONST_CPDAIR ! specific heat of dry air J/kg/K
    tmelt_in = SHR_CONST_TKFRZ  ! Freezing point of Water (K)
    latvap   = SHR_CONST_LATVAP ! latent heat vaporization J/kg
    latice   = SHR_CONST_LATICE ! latent heat freezing J/kg

    ! Initialise microphysics (see micro_mg2_utils.f90 for old values)
    if (micro_unset) then

       rhmini_in = 0.8_wp       ! Minimum rh for ice cloud fraction > 0.

       micro_mg_dcs      = 195.0e-6_wp ! dcs = 90.e-6_r8 
       micro_mg_dcs_tdep = .true.      ! false to get old ice autoconversion

       microp_uniform_in      = .true.  ! true in old MG2 (false in ACME but causes KiD to fail)
       do_cldice_in           = .true.  ! true in old MG2
       use_hetfrz_classnuc_in = .false. ! not an option in old MG2 (requires more interfacing to bring into KiD)

       micro_mg_precip_frac_method_in = 'in_cloud' ! max_overlap assumed in old MG2
       micro_mg_berg_eff_factor_in    = 0.1_wp     ! 1.0 (scaling not used in old MG2)
       allow_sed_supersat_in          = .false.    ! true => saturation adjstment (old MG2)
                                                   ! false => no saturation adjustment
       ice_sed_ai                     = 500.0_wp   ! ai        =   700 in old MG2
       prc_coef1_in                   = 30500.0_wp ! prc_coef1 =  1350 in old MG2
       prc_exp_in                     = 3.19_wp    ! prc_exp   =  2.47 in old MG2
       prc_exp1_in                    = -1.2_wp    ! prc_exp1  = -1.79 in old MG2
       cld_sed_in                     = 1.0_wp     ! acn value =   1.0 in old MG2
       mg_prc_coeff_fix_in            = .true.     ! false in old MG2

       h2otrip_in = 273.16_wp
       tboil_in   = 373.16_wp
       ttrice_in  = 20.0_wp
       epsilo_in  =  SHR_CONST_MWWV / SHR_CONST_MWDAIR  

       call wv_sat_methods_init(kind, tmelt_in, h2otrip_in, tboil_in, &
            ttrice_in, epsilo_in, errstring)

       call micro_mg2_acme_v1beta_init( &
            kind, gravit, rair, rh2o, cpair,    &
            tmelt_in, latvap, latice,           &
            rhmini_in, micro_mg_dcs, micro_mg_dcs_tdep, &
            microp_uniform_in, do_cldice_in, use_hetfrz_classnuc_in, &
            micro_mg_precip_frac_method_in, micro_mg_berg_eff_factor_in, &
            allow_sed_supersat_in, ice_sed_ai, prc_coef1_in,prc_exp_in,  &
            prc_exp1_in, cld_sed_in, mg_prc_coeff_fix_in, errstring)
       !--ag
       micro_unset=.False.
    end if

    !microphysics routine for each timestep goes here...

    mgncol    = nx
    mgcols(:) = 1
    nlev      = nz
    top_lev   = 1  
    deltatin  = dt

    !set dust information (for contact freezing):
    !size
    rndstn(:,:,1)  = 0.5e-6_wp      ! radius (m)
    rndstn(:,:,2)  = 1.0e-6_wp      ! radius (m)
    rndstn(:,:,3)  = 2.0e-6_wp      ! radius (m)
    rndstn(:,:,4)  = 10.0e-6_wp     ! radius (m)
    !number
    naconin(:,:,:) = 1000.0_wp  ! 100 m-3

#ifdef ADJUST_SATURATION_BEFORE
    ! zero saturation adjustment tendencies
    qvlato_adj(i,k)  = 0.0_wp 
    tlato_adj(i,k)   = 0.0_wp 
    qctendo_adj(i,k) = 0.0_wp 
    nctendo_adj(i,k) = 0.0_wp 
    qitendo_adj(i,k) = 0.0_wp 
    nitendo_adj(i,k) = 0.0_wp 

    ! Add saturation adjustment...based on m2005
    do i=1,nx
       do k=1,nz
          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
          ! WATER SATURATION  

          DUMT  = tn(i,k) !+ dt*tlato(i,k) 
          DUMQV = qn(i,k) !+ dt*qvlato(i,k)

          call wv_sat_qsat_water(DUMT, pn(i,k), ESL, QVL, 1)

          ! DUMQSS = qsaturation(tn(k,i),pmb(k,i))
          DUMQSS = QVL

          DUMQC = qcn(i,k) !+ dt*qctendo(i,k)
          DUMQC = MAX(DUMQC,0.0_wp)
          DUMNC = ncn(i,k) !+ dt*nctendo(i,k)

          ! SATURATION ADJUSTMENT FOR LIQUID
          PCC= 0.0_wp

          DUMS = DUMQV-DUMQSS

          CPM = cpair*(1.0_wp+0.887_wp*DUMQC) 

          ! Do condensation only (doesnt seem to work without a positive limit warm2: 
          ! why doesn't evap work?)
          !             if (DUMS.GT.qsmall) &
          PCC = DUMS / (1.+ latvap**2 * DUMQSS/(CPM * rh2o * DUMT**2)) / dt

          ! limit evap to qc
          IF (PCC*dt + DUMQC.LT.qsmall) THEN
             PCC = -DUMQC/dt
          END IF

          !++ag
          satadj(i,k)=PCC
          !--ag

          ! Figure out additional number concentration (assume size of 8 microns...)
          DNC =  3.0_r8 * PCC / (4.0_r8*3.14_r8* 8.0e-6_r8**3*997.0_r8)

          ! make sure it does not get less than zero
          IF (DNC*dt + DUMNC.LT.qsmall) THEN
             DNC = -DUMNC/dt
          END IF

          ! if (i.eq.1.and.k.eq. 25) &
          !      print*, 'z,pmb,DUMT,DUMqv,DUMQS,DUMQC,pcc,dnc,tlato=',&
          !      z(k),pmb(k,i),DUMT,DUMQV,QVL,DUMQC,PCC,DNC,tlato(k,i)
          
          ! apply tendencies
          qvlato_adj(i,k)  = -PCC
          tlato_adj(i,k)   = PCC * latvap/CPM
          qctendo_adj(i,k) = PCC
          nctendo_adj(i,k) = DNC 

          ! limters to make sure if 
          ! (a) non negative mass and number and 
          ! (b) no mass, then no number
          if (qn(i,k) + dt * qvlato(i,k).lt.qsmall) then !+++djg
             qvlato_adj(i,k) = -qn(i,k)/dt
          end if

          if (qcn(i,k) + dt * qctendo_adj(i,k).lt.qsmall) then
             qctendo_adj(i,k)=-qcn(i,k)/dt
             nctendo_adj(i,k)=-ncn(i,k)/dt
          end if

          if (qin(i,k) + dt * qitendo_adj(i,k).lt.qsmall) then
             qitendo_adj(i,k)=-qin(i,k)/dt
             nitendo_adj(i,k)=-nin(i,k)/dt
          end if

          ! perform sequential update
          qn(i,k)  = qn(i,k)  + dt * qvlato_adj(i,k)
          tn(i,k)  = tn(i,k)  + dt * tlato_adj(i,k)
          qcn(i,k) = qcn(i,k) + dt * qctendo_adj(i,k)
          ncn(i,k) = ncn(i,k) + dt * nctendo_adj(i,k)
          qin(i,k) = qin(i,k) + dt * qitendo_adj(i,k)
          nin(i,k) = nin(i,k) + dt * nitendo_adj(i,k)

          if (tn(i,k) < 0.0d0) then
             write(*,*) "Error: negative temperature"
          end if
         
          ! set small values are zero
          if (qn(i,k).lt.qsmall)  qn(i,k)  = 0.0_wp
          if (qcn(i,k).lt.qsmall) qcn(i,k) = 0.0_wp
          if (qcn(i,k).lt.qsmall) ncn(i,k) = 0.0_wp
          if (qin(i,k).lt.qsmall) qin(i,k) = 0.0_wp
          if (qin(i,k).lt.qsmall) nin(i,k) = 0.0_wp
          
       end do
    end do
#endif

    !=================================================
    ! 2. NOW FLIP FROM SURF-FIRST TO TOP-OF-ATMOS FIRST
    !=================================================

    !now: need to flip all inputs in the vertical: MG assumes 1=top , nz=surf
    ! Flipped arrays Input:  tn,qn,qcn,ncn,qin,nin,cldn,liqcldf,icecldf

    flip_pn(:,:)    = pn(:,nz:1:-1)
    flip_pdeln(:,:) = pdeln(:,nz:1:-1)
    ! flip_pint(:,:)  = pint(:,nz:0:-1)
    flip_tn(:,:)    = tn(:,nz:1:-1)
    flip_qn(:,:)    = qn(:,nz:1:-1)
    flip_qcn(:,:)   = qcn(:,nz:1:-1)
    flip_qin(:,:)   = qin(:,nz:1:-1)
    flip_ncn(:,:)   = ncn(:,nz:1:-1)
    flip_nin(:,:)   = nin(:,nz:1:-1)
    flip_qrn(:,:)   = qrn(:,nz:1:-1)
    flip_qsn(:,:)   = qsn(:,nz:1:-1)
    flip_nrn(:,:)   = nrn(:,nz:1:-1)
    flip_nsn(:,:)   = nsn(:,nz:1:-1)
    flip_cldn(:,:)  = cldn(:,nz:1:-1)
    flip_liqcldf(:,:) = liqcldf(:,nz:1:-1)
    flip_icecldf(:,:) = icecldf(:,nz:1:-1)

    !Flipped activated aerosol number concentration
    flip_naain(:,:)=naain(:,nz:1:-1)
    flip_npccnin(:,:)=npccnin(:,nz:1:-1)

    !Flipped relative variance (qcvar) and accretion enhancement (new for MG2 after CAM5_3_01)
    flip_relvarn(:,:)      = 0.0_wp
    flip_accre_enhann(:,:) = 1.0_wp

    !=================================================
    ! 6. NOW COMPUTE MICRO TENDS:
    !=================================================
    !Authors: Hugh Morrison, Andrew Gettelman, NCAR, Peter Caldwell, LLNL
    ! e-mail: morrison@ucar.edu, andrew@ucar.edu, caldwell19@llnl.gov

    call micro_mg2_acme_v1beta_tend ( &
         ! Input 
         mgncol, mgcols,               nlev, top_lev,                    &
         deltatin,                                                       &
         flip_tn,                      flip_qn,                          &
         flip_qcn,                     flip_qin,                         &
         flip_ncn,                     flip_nin,                         &
         flip_qrn,                     flip_qsn,                         &
         flip_nrn,                     flip_nsn,                         &
         flip_relvarn,                 flip_accre_enhann,                &
         flip_pn,                      flip_pdeln,                       &
         flip_cldn,                    flip_liqcldf,       flip_icecldf, &
         ! Output
         qcsinksum_rate1ordo,                                            &
         flip_naain,                   flip_npccnin,                     &
         rndstn,                       naconin,                          &
         flip_tlato,                   flip_qvlato,                      &
         flip_qctendo,                 flip_qitendo,                     &
         flip_nctendo,                 flip_nitendo,                     &
         flip_qrtendo,                 flip_qstendo,                     &
         flip_nrtendo,                 flip_nstendo,                     &
         effco,                        effco_fn,           effio,        &
         precto,                       precio,                           &
         nevapro,                      evapsnowo,                        &
         praino,                       prodsnowo,                        &
         cmeouto,                      deffio,                           &
         pgamrado,                     lamcrado,                         &
         qsouto,                       dsouto,                           &
         rflxo,                        sflxo,              flip_qrouto,  &
         reff_raino,                   reff_snowo,                       &
         qcsevapo,                     qisevapo,           qvreso,       &
         flip_cmeiout,                 flip_vtrmco,        flip_vtrmio,  &
         flip_umso,                    flip_umro,                        &
         flip_qcsedteno,               flip_qisedteno,                   &
         flip_qrsedteno,               flip_qssedteno,                   &
         flip_prao,                    flip_prco,                        &
         flip_mnuccco,                 flip_mnuccto,       flip_msacwio, &
         flip_psacwso,                 flip_bergso,        flip_bergo,   &
         flip_melto,                   flip_homoo,                       &
         flip_qcreso,                  flip_prcio,         flip_praio,   &
         flip_qireso,                  flip_mnuccro,       flip_pracso,  &
         flip_meltsdto,                flip_frzrdto,       flip_mnuccdo, &
         flip_nrouto,                  nsouto,                           &
         reflo,                        areflo,             areflzo,      &
         freflo,                       csrflo,             acsrflo,      &
         fcsrflo,                      rercldo,                          &
         ncaio,                        ncalo,                            &
         flip_qrouto2,                 qsouto2,                          &
         flip_nrouto2,                 nsouto2,                          &
         drouto2,                      dsouto2,                          &
         freqso,                       freqro,                           &
         nficeo,                       qcrato,                           &
         tnd_qsnown,                   tnd_nsnown,         re_icen,      &
         errstring,                                                      &
         prer_evapo,                                                     &
         frzimmn,                      frzcntn,            frzdepn)

    !=================================================
    ! 7. NOW FLIP BACK TO SURF-FIRST (KiD ORDERING):
    !=================================================
    ! Key tendencies to hook up prognostically:
    !tlato,qvlato,qctendo,qitendo,nctendo,nitendo

    !note: flip_tlato is dry static energy: W/kg divide by cpair
    tlato(:,:)   = flip_tlato(:,nz:1:-1) / cpair
    qvlato(:,:)  = flip_qvlato(:,nz:1:-1)
    qctendo(:,:) = flip_qctendo(:,nz:1:-1)
    qitendo(:,:) = flip_qitendo(:,nz:1:-1)
    nctendo(:,:) = flip_nctendo(:,nz:1:-1)
    nitendo(:,:) = flip_nitendo(:,nz:1:-1)
    qrtendo(:,:) = flip_qrtendo(:,nz:1:-1)
    qstendo(:,:) = flip_qstendo(:,nz:1:-1)
    nrtendo(:,:) = flip_nrtendo(:,nz:1:-1)
    nstendo(:,:) = flip_nstendo(:,nz:1:-1)

    nrouto(:,:)  = flip_nrouto2(:,nz:1:-1)
    qrouto(:,:)  = flip_qrouto2(:,nz:1:-1)
    prco(:,:)    = flip_prco(:,nz:1:-1)
    prao(:,:)    = flip_prao(:,nz:1:-1)
    cmeiout(:,:) = flip_cmeiout(:,nz:1:-1)
    umso(:,:)    = flip_umso(:,nz:1:-1)
    umro(:,:)    = flip_umro(:,nz:1:-1)
    vtrmco(:,:)  = flip_vtrmco(:,nz:1:-1)
    vtrmio(:,:)  = flip_vtrmio(:,nz:1:-1)

    qcsedteno(:,:) = flip_qcsedteno(:,nz:1:-1)
    qisedteno(:,:) = flip_qisedteno(:,nz:1:-1)
    mnuccco(:,:)   = flip_mnuccco(:,nz:1:-1)  
    mnuccto(:,:)   = flip_mnuccto(:,nz:1:-1) 
    msacwio(:,:)   = flip_msacwio(:,nz:1:-1)  
    psacwso(:,:)   = flip_psacwso(:,nz:1:-1)  
    bergso(:,:)    = flip_bergso(:,nz:1:-1)   
    bergo(:,:)     = flip_bergo(:,nz:1:-1)    
    melto(:,:)     = flip_melto(:,nz:1:-1)    
    homoo(:,:)     = flip_homoo(:,nz:1:-1)  
    qcreso(:,:)    = flip_qcreso(:,nz:1:-1)             
    prcio(:,:)     = flip_prcio(:,nz:1:-1)    
    praio(:,:)     = flip_praio(:,nz:1:-1)    
    qireso(:,:)    = flip_qireso(:,nz:1:-1)             
    mnuccro(:,:)   = flip_mnuccro(:,nz:1:-1)  
    pracso(:,:)    = flip_pracso(:,nz:1:-1)   
    meltsdto(:,:)  = flip_meltsdto(:,nz:1:-1) 
    frzrdto(:,:)   = flip_frzrdto(:,nz:1:-1) 
    mnuccdo(:,:)   = flip_mnuccdo(:,nz:1:-1)       

#ifdef ADJUST_SATURATION_AFTER
    !Add saturation adjustment...based on m2005
    do i=1,nx
       do k=1,nz
          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
          ! WATER SATURATION  

          DUMT = tn(i,k) + dt*tlato(i,k)
          DUMQV = qn(i,k) + dt*qvlato(i,k)

          call wv_sat_qsat_water(DUMT, pn(i,k), ESL, QVL, 1)

          !               DUMQSS = qsaturation(tn(k,i),pmb(k,i))
          DUMQSS = QVL

          DUMQC = qcn(i,k) + dt*qctendo(i,k)
          DUMQC = MAX(DUMQC,0.0_wp)
          DUMNC = ncn(i,k) + dt*nctendo(i,k)

          ! SATURATION ADJUSTMENT FOR LIQUID
          PCC= 0.0_wp

          DUMS = DUMQV-DUMQSS

          CPM = cpair*(1.0_wp+0.887_wp*DUMQC) 

          !Do condensation only (doesnt seem to work without a positive limit warm2: why doesn't evap work?)
          !             if (DUMS.GT.qsmall) &
          PCC = DUMS / (1.+ latvap**2 * DUMQSS/(CPM * rh2o * DUMT**2)) / dt

          !limit evap to qc
          IF (PCC*dt + DUMQC.LT.qsmall) THEN
             PCC = -DUMQC/dt
          END IF

          !++ag
          satadj(i,k)=PCC
          !--ag

          !Figure out additional number concentration (assume size of 8 microns...)
          DNC =  3.0_r8 * PCC / (4.0_r8*3.14_r8* 8.0e-6_r8**3*997.0_r8)

          !make sure it does not get less than zero
          IF (DNC*dt + DUMNC.LT.qsmall) THEN
             DNC = -DUMNC/dt
          END IF


          !             if (i.eq.1.and.k.eq. 25) &
          !                  print*, 'z,pmb,DUMT,DUMqv,DUMQS,DUMQC,pcc,dnc,tlato=',&
          !                  z(k),pmb(k,i),DUMT,DUMQV,QVL,DUMQC,PCC,DNC,tlato(k,i)

          !apply tendencies
          qvlato(i,k) = qvlato(i,k) - PCC
          tlato(i,k) = tlato(i,k) + PCC * latvap/CPM
          qctendo(i,k) = qctendo(i,k) + PCC
          nctendo(i,k) = nctendo(i,k) + DNC 

          !limters to make sure if (a) non negative mass and number and (b) no mass, then no number
          if (qcn(i,k) + dt * qctendo(i,k).lt.qsmall) then
             qctendo(i,k)=-qcn(i,k)/dt
             nctendo(i,k)=-ncn(i,k)/dt
          end if

          if (qin(i,k) + dt * qitendo(i,k).lt.qsmall) then
             qitendo(i,k)=-qin(i,k)/dt
             nitendo(i,k)=-nin(i,k)/dt
          end if
       end do
    end do
#endif

#ifdef ADJUST_SATURATION_BEFORE
    do i=1,nx
       do k=1,nz
          ! add in saturation from before microphysics to output tendencies
          qvlato(i,k)  = qvlato(i,k)  + qvlato_adj(i,k) 
          tlato(i,k)   = tlato(i,k)   + tlato_adj(i,k)  
          qctendo(i,k) = qctendo(i,k) + qctendo_adj(i,k)
          nctendo(i,k) = nctendo(i,k) + nctendo_adj(i,k)
          qitendo(i,k) = qitendo(i,k) + qitendo_adj(i,k)
          nitendo(i,k) = nitendo(i,k) + nitendo_adj(i,k)

          ! limiters?
       end do
    end do
#endif   

    !=================================================
    ! 9. WRITE OUTPUT:
    !=================================================

    precprt2d = precto
    snowrt2d  = precio

    ! save tendencies
    do i=1,nx
       do k=1,nz
          dtheta_mphys(k,i) = tlato(i,k)/exner(k,i)
          dqv_mphys(k,i)    = qvlato(i,k)

          dhydrometeors_mphys(k,i,1)%moments(1,1) = qctendo(i,k)
          dhydrometeors_mphys(k,i,1)%moments(1,2) = nctendo(i,k)
          dhydrometeors_mphys(k,i,3)%moments(1,1) = qitendo(i,k)
          dhydrometeors_mphys(k,i,3)%moments(1,2) = nitendo(i,k)
          dhydrometeors_mphys(k,i,2)%moments(1,1) = qrtendo(i,k)
          dhydrometeors_mphys(k,i,2)%moments(1,2) = nrtendo(i,k)
          dhydrometeors_mphys(k,i,4)%moments(1,1) = qstendo(i,k)
          dhydrometeors_mphys(k,i,4)%moments(1,2) = nstendo(i,k)
       end do
    end do

    ! Save some diagnostics
    !Rain/Snow (what are units?)

    name='total_surface_ppt'
    units='m s-1'

    if (nx == 1) then
       call save_dg(precprt2d(1), name, i_dgtime,  units, dim='time')
    else
       call save_dg(precprt2d(1:nx), name, i_dgtime,  units, dim='time')
    endif

    name='surface_ppt_for_snow'
    units='m s-1'
    if (nx == 1) then
       call save_dg(snowrt2d(1), name, i_dgtime,  units, dim='time')
    else
       call save_dg(snowrt2d(1:nx), name, i_dgtime,  units, dim='time')
    endif

    !    name='rain_mass'
    !    units='kg kg-1'
    !    if (nx == 1) then
    !       call save_dg(qrouto(1,:), name, i_dgtime,  units, dim='z')
    !    else
    !       call save_dg(qrouto(1:nx,:), name, i_dgtime,  units, dim='z')
    !    endif   

    !    name='rain_number'
    !    units='kg-1'
    !    if (nx == 1) then
    !       call save_dg(nrouto(1,:), name, i_dgtime,  units, dim='z')
    !    else
    !       call save_dg(nrouto(1:nx,:), name, i_dgtime,  units, dim='z')
    !    endif   


    name='autoconversion_rate'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(prco(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(prco(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='accretion_rate'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(prao(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(prao(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='ice_condensation_rate'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(cmeiout(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(cmeiout(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='mass_weighted_liquid_fallspeed'
    units='m s-1'
    if (nx == 1) then
       call save_dg(vtrmco(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(vtrmco(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='mass_weighted_ice_fallspeed'
    units='m s-1'
    if (nx == 1) then
       call save_dg(vtrmio(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(vtrmio(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='mass_weighted_rain_fallspeed'
    units='m s-1'
    if (nx == 1) then
       call save_dg(umro(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(umro(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='mass_weighted_snow_fallspeed'
    units='m s-1'
    if (nx == 1) then
       call save_dg(umso(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(umso(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='saturation_adjustment_tendency'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(satadj(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(satadj(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! MICROPHYSICS BUDGET TERMS....use same names as CAM interface (helps with budget code)

    name='QCSEDTENO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(qcsedteno(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(qcsedteno(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='QISEDTENO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(qisedteno(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(qisedteno(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='PRCO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(prco(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(prco(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='PRAO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(prao(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(prao(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='MNUCCCO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(mnuccco(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(mnuccco(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='MNUCCTO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(mnuccto(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(mnuccto(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='MSACWIO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(msacwio(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(msacwio(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='PSACWSO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(psacwso(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(psacwso(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='BERGSO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(bergso(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(bergso(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='BERGO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(bergo(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(bergo(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='MELTO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(melto(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(melto(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='HOMOO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(homoo(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(homoo(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='QCRESO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(qcreso(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(qcreso(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='PRCIO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(prcio(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(prcio(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='PRAIO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(praio(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(praio(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='QIRESO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(qireso(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(qireso(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='MNUCCRO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(mnuccro(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(mnuccro(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='PRACSO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(pracso(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(pracso(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='MELTSDTO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(meltsdto(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(meltsdto(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='FRZRDTO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(frzrdto(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(frzrdto(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='MNUCCDO'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(mnuccdo(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(mnuccdo(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    name='MPDLIQ'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(qctendo(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(qctendo(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------
    ! output limiter information
    ! ---------------------------------------------------------------------------
     
    flip_data(:,:) = qric_limiter(:,nz:1:-1)
    name='qric_lim'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = qric_limiter_mag(:,nz:1:-1)
    name='qric_lim_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = qric_qsmall(:,nz:1:-1)
    name='qric_qsmall'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = qric_qsmall_mag(:,nz:1:-1)
    name='qric_qsmall_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = nric_qsmall(:,nz:1:-1)
    name='nric_qsmall'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = nric_qsmall_mag(:,nz:1:-1)
    name='nric_qsmall_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = nric_neg(:,nz:1:-1)
    name='nric_neg'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = nric_neg_mag(:,nz:1:-1)
    name='nric_neg_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = qsic_limiter(:,nz:1:-1)
    name='qsic_lim'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = qsic_limiter_mag(:,nz:1:-1)
    name='qsic_lim_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = qsic_qsmall(:,nz:1:-1)
    name='qsic_qsmall'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = qsic_qsmall_mag(:,nz:1:-1)
    name='qsic_qsmall_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = nric_qsmall(:,nz:1:-1)
    name='nsic_qsmall'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = nric_qsmall_mag(:,nz:1:-1)
    name='nsic_qsmall_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = nsic_neg(:,nz:1:-1)
    name='nsic_neg'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = nsic_neg_mag(:,nz:1:-1)
    name='nsic_neg_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = qc_conservation(:,nz:1:-1)
    name='qc_conserv'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = qc_conservation_mag(:,nz:1:-1)
    name='qc_conserv_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = ice_nucleation_limiter(:,nz:1:-1)
    name='ice_nuc_lim'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = ice_nucleation_limiter_mag(:,nz:1:-1)
    name='ice_nuc_lim_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = ice_deposition_limiter(:,nz:1:-1)
    name='ice_dep_lim'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = ice_deposition_limiter_mag(:,nz:1:-1)
    name='ice_dep_lim_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = nc_conservation(:,nz:1:-1)
    name='nc_conserv'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = nc_conservation_mag(:,nz:1:-1)
    name='nc_conserv_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = qr_conservation(:,nz:1:-1)
    name='qr_conserv'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = qr_conservation_mag(:,nz:1:-1)
    name='qr_conserv_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = nr_conservation(:,nz:1:-1)
    name='nr_conserv'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = nr_conservation_mag(:,nz:1:-1)
    name='nr_conserv_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = qi_conservation(:,nz:1:-1)
    name='qi_conserv'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = qi_conservation_mag(:,nz:1:-1)
    name='qi_conserv_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = ni_conservation(:,nz:1:-1)
    name='ni_conserv'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = ni_conservation_mag(:,nz:1:-1)
    name='ni_conserv_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = qs_conservation(:,nz:1:-1)
    name='qs_conserv'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = qs_conservation_mag(:,nz:1:-1)
    name='qs_conserv_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    ! flip_data(:,:) = ns_conservationq(:,nz:1:-1)
    ! name='ns_conserv'
    ! units=''
    ! if (nx == 1) then
    !    call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    ! else
    !    call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    ! endif

    ! flip_data(:,:) = ns_conservation_mag(:,nz:1:-1)
    ! name='ns_conserv_mag'
    ! units=''
    ! if (nx == 1) then
    !    call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    ! else
    !    call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    ! endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = qiqs_sublimation_qr_evaporation_limiter(:,nz:1:-1)
    name='qiqs_sub_qr_evap_lim'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = rain_evaporation_limiter_mag(:,nz:1:-1)
    name='rain_evap_lim_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = snow_sublimation_limiter_mag(:,nz:1:-1)
    name='snow_sub_lim_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = ice_sublimation_limiter_mag(:,nz:1:-1)
    name='ice_sub_lim_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    flip_data(:,:) = ni_tendency_limiter(:,nz:1:-1)
    name='ni_tend_lim'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    flip_data(:,:) = ni_tendency_limiter_mag(:,nz:1:-1)
    name='ni_tend_lim_mag'
    units=''
    if (nx == 1) then
       call save_dg(flip_data(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(flip_data(1:nx,:), name, i_dgtime,  units, dim='z')
    endif

    ! ---------------------------------------------------------------------------

    name='nsteps_qi'
    units='#'
    if (nx == 1) then
       call save_dg(1.0*nsteps_qi, name, i_dgtime,  units, dim='time')
    endif

    name='nsteps_qc'
    units='#'
    if (nx == 1) then
       call save_dg(1.0*nsteps_qc, name, i_dgtime,  units, dim='time')
    endif

    name='nsteps_qr'
    units='#'
    if (nx == 1) then
       call save_dg(1.0*nsteps_qr, name, i_dgtime,  units, dim='time')
    endif

    name='nsteps_qs'
    units='#'
    if (nx == 1) then
       call save_dg(1.0*nsteps_qs, name, i_dgtime,  units, dim='time')
    endif

    ! ---------------------------------------------------------------------------  

  end Subroutine mphys_mg2_acme_v1beta_interface

end module mphys_mg2_acme_v1beta
