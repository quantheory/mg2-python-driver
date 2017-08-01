module Sedimentation
  use micro_mg2_acme_v1beta_utils, only: r8
  implicit none
  private
  save
  public :: sed_CalcFallVelocity, sed_AdvanceOneStep
  integer, parameter, public :: MG_LIQUID = 1
  integer, parameter, public :: MG_ICE = 2
  integer, parameter, public :: MG_RAIN = 3
  integer, parameter, public :: MG_SNOW = 4

contains

  subroutine sed_CalcFallVelocity(q,qtend,n,ntend,cloud_frac,rho,nlev,i, &
    mg_type,g,an,rhof,fq,fn,ncons,nnst,gamma_b_plus1,gamma_b_plus4)
    use micro_mg2_acme_v1beta_utils, only: size_dist_param_liq, &
                                           size_dist_param_basic, &
                                           mg_ice_props, mg_liq_props, &
                                           mg_snow_props, mg_rain_props, &
                                           qsmall, bi, bc, br, bs
    implicit none
    real(r8), intent(in)            :: q(:,:), qtend(:,:), n(:,:), ntend(:,:)
    real(r8), intent(in)            :: rho(:,:), an(:,:), rhof(:,:)
    real(r8), intent(in)            :: cloud_frac(:,:), g
    integer,  intent(in)            :: nlev, i, mg_type
    real(r8), intent(in), optional  :: nnst, gamma_b_plus1, gamma_b_plus4
    logical, intent(in), optional   :: ncons
    real(r8), intent(out)           :: fq(:), fn(:)

    real(r8) :: qtemp(nlev), ntemp(nlev), lam(nlev), pgam(nlev)
    integer :: k, nstep

    do k=1,nlev

      ! calculate sedimentation for cloud water and ice
      !================================================================================

      ! update in-cloud cloud mixing ratio and number concentration
      ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
      ! note: these are in-cloud values***, hence we divide by cloud fraction

      qtemp(k) = q(i,k)/cloud_frac(i,k)
      ntemp(k) = max(n(i,k)/cloud_frac(i,k),0._r8)

      ! switch for specification of droplet and crystal number or ice cloud number
      if ((mg_type == MG_LIQUID  .or. mg_type == MG_ICE) .and. ncons) then
        ntemp(k) = nnst/rho(i,k)
      end if

      if(mg_type == MG_LIQUID) then
        ! obtain new slope parameter to avoid possible singularity
        ! NOTE: original variables were pgam(i,k), lamc(i,k)
        call size_dist_param_liq(mg_liq_props, qtemp(k), ntemp(k), rho(i,k), &
             pgam(k), lam(k))
        ! calculate number and mass weighted fall velocity for droplets and cloud ice
        !-------------------------------------------------------------------
        if (qtemp(k) .ge. qsmall) then
          ! NOTE: original variable was vtrmc(i,k) = everything after g*rho
          fq(k) = g*rho(i,k)*an(i,k)*gamma(4._r8+bc+pgam(k))/ &
               (lam(k)**bc*gamma(pgam(k)+4._r8))

          fn(k) = g*rho(i,k)* &
               an(i,k)*gamma(1._r8+bc+pgam(k))/ &
               (lam(k)**bc*gamma(pgam(k)+1._r8))
        else
          fq(k) = 0._r8
          fn(k)= 0._r8
        end if

      else if (mg_type == MG_ICE) then
        ! obtain new slope parameter to avoid possible singularity
        ! NOTE: original variables were ntemp(i,k),lami(i,k)
        call size_dist_param_basic(mg_ice_props, qtemp(k), ntemp(k), lam(k))
        ! calculate number and mass weighted fall velocity for cloud ice
        if (qtemp(k) .ge. qsmall) then
          ! NOTE: original variable was vtrmi(i,k) = everything after g*rho
          fq(k) = g*rho(i,k)*min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bi), &
                1.2_r8*rhof(i,k))
          fn(k) = g*rho(i,k)* &
                min(an(i,k)*gamma_b_plus1/lam(k)**bi,1.2_r8*rhof(i,k))
        else
          fq(k) = 0._r8
          fn(k)= 0._r8
        end if

      else if (mg_type == MG_RAIN) then
        ! fallspeed for rain
        ! NOTE: original variables were ntemp(i,k),lamr(i,k)
        call size_dist_param_basic(mg_rain_props, qtemp(k), ntemp(k), lam(k))
        if (lam(k) .ge. qsmall) then
           ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)
           ! NOTE: original variables were unr(i,k) and umr(i,k) = everything after g*rho
           fq(k) = g*rho(i,k)* &
              min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**br),9.1_r8*rhof(i,k))
           fn(k) = g*rho(i,k)* &
              min(an(i,k)*gamma_b_plus1/lam(k)**br,9.1_r8*rhof(i,k))
        else
           fq(k)=0._r8
           fn(k)=0._r8
        end if

      else if (mg_type == MG_SNOW) then
        ! fallspeed for snow
        ! NOTE: original variables were ntemp(i,k),lams(i,k)
        call size_dist_param_basic(mg_snow_props, qtemp(k), ntemp(k), lam(k))
        if (lam(k) .ge. qsmall) then
           ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
           ! NOTE: original variables were uns(i,k) and ums(i,k) = everything after g*rho
           fq(k) = g*rho(i,k)* &
              min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bs),1.2_r8*rhof(i,k))
           fn(k) = g*rho(i,k)* &
              min(an(i,k)*gamma_b_plus1/lam(k)**bs,1.2_r8*rhof(i,k))
        else
           fq(k)=0._r8
           fn(k)=0._r8
        end if
      end if

    end do

  end subroutine sed_CalcFallVelocity

  subroutine sed_AdvanceOneStep(q,fq,n,fn,pdel,deltat,deltat_sed,nlev,i,mg_type,g, &
    qtend,ntend,prect,qsedtend,cloud_frac,qvlat,tlat,xxl,preci,qsevap)
    implicit none
    real(r8), intent(in)              :: fq(:), fn(:), pdel(:,:)
    real(r8), intent(in)              :: deltat, deltat_sed, g
    integer, intent(in)               :: nlev, i, mg_type
    real(r8), intent(inout)           :: q(:,:), qtend(:,:), n(:,:), ntend(:,:)
    real(r8), intent(inout)           :: prect(:)
    real(r8), intent(in), optional    :: cloud_frac(:,:), xxl
    real(r8), intent(inout), optional :: qvlat(:,:), tlat(:,:), preci(:)
    real(r8), intent(out), optional   :: qsevap(:,:)
    real(r8), intent(out)             :: qsedtend(:,:)

    real(r8) :: fluxQ(0:nlev), fluxN(0:nlev), ratio(nlev), deltaFluxQ, deltaFluxN
    real(r8) :: deltafluxQ_evap, dfluxp, dfluxm
    integer :: k

    ! set physical flux at top boundary (bottom boundary is outflow)
    fluxQ(0) = 0._r8
    fluxN(0) = 0._r8

    ! compute cell interface fluxes using 2nd order ENO scheme
    fluxQ(1) = fq(1)*q(i,1)
    fluxN(1) = fn(1)*n(i,1)
    fluxQ(nlev) = fq(nlev)*q(i,nlev)
    fluxN(nlev) = fn(nlev)*n(i,nlev)

    do k = 2,nlev-1
      dfluxp = (fq(k+1)*q(i,k+1)-fq(k)*q(i,k))/(0.5_r8*(pdel(i,k+1)+pdel(i,k)))
      dfluxm = (fq(k)*q(i,k)-fq(k-1)*q(i,k-1))/(0.5_r8*(pdel(i,k)+pdel(i,k-1)))
      if (abs(dfluxm) < abs(dfluxp)) then
        fluxQ(k) = fq(k)*q(i,k) + 0.5_r8*pdel(i,k)*dfluxm
      else
        fluxQ(k) = fq(k)*q(i,k) + 0.5_r8*pdel(i,k)*dfluxp
      end if
      dfluxp = (fn(k+1)*n(i,k+1)-fn(k)*n(i,k))/(0.5_r8*(pdel(i,k+1)+pdel(i,k)))
      dfluxm = (fn(k)*n(i,k)-fn(k-1)*n(i,k-1))/(0.5_r8*(pdel(i,k)+pdel(i,k-1)))
      if (abs(dfluxm) < abs(dfluxp)) then
        fluxN(k) = fn(k)*n(i,k) + 0.5_r8*pdel(i,k)*dfluxm
      else
        fluxN(k) = fn(k)*n(i,k) + 0.5_r8*pdel(i,k)*dfluxp
      end if
    end do

    ! for cloud liquid and ice, if cloud fraction increases with height
    ! then add flux from above to both vapor and cloud water of current level
    ! this means that flux entering clear portion of cell from above evaporates
    ! instantly
    ! note: this is not an issue with precip, since we assume max overlap
    if (mg_type == MG_ICE .or. mg_type == MG_LIQUID) then
      ratio(1) = 1._r8
      ratio(2:nlev) = cloud_frac(i,2:nlev)/cloud_frac(i,1:(nlev-1))
    else
      ratio = 1._r8
    end if

    do k = 1,nlev
      ratio(k)=min(ratio(k),1._r8)
      deltafluxQ = (fluxQ(k)-ratio(k)*fluxQ(k-1))/pdel(i,k)
      deltafluxN = (fluxN(k)-ratio(k)*fluxN(k-1))/pdel(i,k)
      ! add fallout terms to eulerian tendencies
      qtend(i,k) = qtend(i,k) - deltafluxQ*deltat_sed/deltat
      ntend(i,k) = ntend(i,k) - deltafluxN*deltat_sed/deltat
      ! sedimentation tendency for output
      qsedtend(i,k) = qsedtend(i,k) - deltafluxQ*deltat_sed/deltat
      if (mg_type == MG_ICE .or. mg_type == MG_LIQUID) then
        ! add terms to to evap/sub of cloud water
        deltafluxQ_evap = (ratio(k)-1._r8)*fluxQ(k-1)/pdel(i,k)
        qvlat(i,k) = qvlat(i,k) - deltafluxQ_evap*deltat_sed/deltat
        qsevap(i,k) = qsevap(i,k) - deltafluxQ_evap*deltat_sed/deltat
        tlat(i,k) = tlat(i,k) + deltafluxQ_evap*xxl*deltat_sed/deltat
      end if

      q(i,k) = q(i,k) - deltat_sed*deltafluxQ
      n(i,k) = n(i,k) - deltat_sed*deltafluxN

      if (q(i,k) < 0._r8) then
        stop "q negative"
      else if (n(i,k) < 0._r8) then
        stop "n negative"
      end if
    end do

    ! units below are m/s
    ! sedimentation flux at surface is added to precip flux at surface
    ! to get total precip (cloud + precip water) rate

    prect(i) = prect(i) + deltat_sed/deltat*fluxQ(nlev)/g/1000._r8
    if (mg_type == MG_ICE .or. mg_type == MG_SNOW) then
      preci(i) = preci(i) + deltat_sed/deltat*fluxQ(nlev)/g/1000._r8
    end if

  end subroutine

end module Sedimentation
