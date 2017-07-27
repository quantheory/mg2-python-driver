module Sedimentation
  use micro_mg2_acme_v1beta_utils, only: r8
  implicit none
  private
  save
  public :: sed_CalcFallVelocity
  integer, parameter, public :: MG_LIQUID = 1
  integer, parameter, public :: MG_ICE = 2
  integer, parameter, public :: MG_RAIN = 3
  integer, parameter, public :: MG_SNOW = 4

contains

  subroutine sed_CalcFallVelocity(q,qtend,n,ntend,cloud_frac,rho,deltat,nlev,i, &
    mg_type,g,an,rhof,fq,fn,ncons,nnst,gamma_b_plus1,gamma_b_plus4)
    use micro_mg2_acme_v1beta_utils, only: size_dist_param_liq, &
                                           size_dist_param_basic, &
                                           mg_ice_props, mg_liq_props, &
                                           mg_snow_props, mg_rain_props, &
                                           qsmall, bi, bc, br, bs
    implicit none
    real(r8), intent(in)            :: q(:,:), qtend(:,:), ntend(:,:)
    real(r8), intent(in)            :: rho(:,:), an(:,:), rhof(:,:)
    real(r8), intent(in)            :: cloud_frac(:,:)
    real(r8), intent(in)            :: deltat, g
    integer,  intent(in)            :: nlev, i, mg_type
    real(r8), intent(inout)         :: n(:,:) ! limiter in size_dist_param_*
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

      qtemp(k) = (q(i,k)+qtend(i,k)*deltat)/cloud_frac(i,k)
      ntemp(k) = max((n(i,k)+ntend(i,k)*deltat)/cloud_frac(i,k),0._r8)

      ! switch for specification of droplet and crystal number or ice cloud number
      if ((mg_type == MG_LIQUID  .or. mg_type == MG_ICE) .and. ncons) then
        ntemp(k) = nnst/rho(i,k)
      end if

      if(mg_type == MG_LIQUID) then
        ! obtain new slope parameter to avoid possible singularity
        call size_dist_param_liq(mg_liq_props, qtemp(k), ntemp(k), rho(i,k), &
             pgam(k), lam(k))
        ! calculate number and mass weighted fall velocity for droplets and cloud ice
        !-------------------------------------------------------------------
        if (qtemp(k) .ge. qsmall) then
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
        call size_dist_param_basic(mg_ice_props, qtemp(k), ntemp(k), lam(k))
        ! calculate number and mass weighted fall velocity for cloud ice
        if (qtemp(k) .ge. qsmall) then
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
        call size_dist_param_basic(mg_rain_props, qtemp(k), ntemp(k), lam(k))
        if (lam(k) .ge. qsmall) then
           ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)
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
        call size_dist_param_basic(mg_snow_props, qtemp(k), ntemp(k), lam(k))
        if (lam(k) .ge. qsmall) then
           ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
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

  end subroutine

end module
