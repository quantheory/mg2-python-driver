module Sedimentation
  use micro_mg2_acme_v1beta_utils, only: r8
  implicit none
  private
  save
  public :: sed_CalcFallVelocity, sed_AdvanceOneStep
  integer, parameter, public  :: MG_LIQUID = 1
  integer, parameter, public  :: MG_ICE = 2
  integer, parameter, public  :: MG_RAIN = 3
  integer, parameter, public  :: MG_SNOW = 4
  real(r8), parameter, public :: CFL = 0.99_r8

contains

  subroutine sed_CalcFallVelocity(q,qtend,n,ntend,cloud_frac,rho,pdel,nlev,i, &
    mg_type,deltat,g,an,rhof,fq,fn,cfl,ncons,nnst,gamma_b_plus1,gamma_b_plus4)
    use micro_mg2_acme_v1beta_utils, only: size_dist_param_liq, &
                                           size_dist_param_basic, &
                                           mg_ice_props, mg_liq_props, &
                                           mg_snow_props, mg_rain_props, &
                                           qsmall, bi, bc, br, bs
    implicit none
    real(r8), intent(in)            :: q(:,:), qtend(:,:), n(:,:), ntend(:,:)
    real(r8), intent(in)            :: rho(:,:), pdel(:,:), an(:,:), rhof(:,:)
    real(r8), intent(in)            :: cloud_frac(:,:), deltat, g
    integer,  intent(in)            :: nlev, i, mg_type
    real(r8), intent(in), optional  :: nnst, gamma_b_plus1, gamma_b_plus4
    logical, intent(in), optional   :: ncons
    real(r8), intent(out)           :: cfl, fq(:), fn(:)

    real(r8) :: qtemp(nlev), ntemp(nlev), lam(nlev), pgam(nlev)
    real(r8) :: qic_dlam(2,nlev), nic_dlam(2,nlev)
    real(r8) :: qdfq(2), ndfn(2), s1, s2, tr, det
    integer :: k, nstep

    ! initialize values to zero
    fq = 0._r8
    fn = 0._r8
    cfl = 0._r8

    ! use quantity in cloud
    qtemp = q(i,:)/cloud_frac(i,:)
    ntemp = n(i,:)/cloud_frac(i,:)
    where(ntemp < 0._r8) ntemp = 0._r8
    ! switch for specification of droplet and crystal number or ice cloud number
    if ((mg_type == MG_LIQUID  .or. mg_type == MG_ICE) .and. ncons) then
      ntemp = nnst/rho(i,:)
    end if

    ! compute lam and pgam using elemental function call
    if (mg_type == MG_LIQUID) then
      call size_dist_param_liq(mg_liq_props, qtemp(:), ntemp(:), rho(i,:), &
                               pgam(:), lam(:))

    else if (mg_type == MG_ICE) then
      call size_dist_param_basic(mg_ice_props, qtemp(:), ntemp(:), lam(:))

    else if (mg_type == MG_RAIN) then
      call size_dist_param_basic(mg_rain_props, qtemp(:), ntemp(:), lam(:), &
                                 qic_dlam_dqic=qic_dlam(1,:), &
                                 nic_dlam_dqic=nic_dlam(1,:), &
                                 qic_dlam_dnic=qic_dlam(2,:), &
                                 nic_dlam_dnic=nic_dlam(2,:))

    else if (mg_type == MG_SNOW) then
       call size_dist_param_basic(mg_snow_props, qtemp(:), ntemp(:), lam(:))
    else
       print *, "Invalid mg_type to mg2_sedimentation"
       stop
    end if

    ! Loop through levels to compute fq and fn while updating the CFL number
    do k=1,nlev

      ! initialize derivative quantities to zero
      qdfq = 0.d0
      ndfn = 0.d0

      if(mg_type == MG_LIQUID) then
        if (qtemp(k) .ge. qsmall) then
          fq(k) = g*rho(i,k)*an(i,k)*gamma(4._r8+bc+pgam(k))/ &
               (lam(k)**bc*gamma(pgam(k)+4._r8))

          fn(k) = g*rho(i,k)* &
               an(i,k)*gamma(1._r8+bc+pgam(k))/ &
               (lam(k)**bc*gamma(pgam(k)+1._r8))
        end if

      else if (mg_type == MG_ICE) then
        if (qtemp(k) .ge. qsmall) then
          fq(k) = g*rho(i,k)*min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bi), &
                1.2_r8*rhof(i,k))
          fn(k) = g*rho(i,k)* &
                min(an(i,k)*gamma_b_plus1/lam(k)**bi,1.2_r8*rhof(i,k))
        end if

      else if (mg_type == MG_RAIN) then
        if (lam(k) .ge. qsmall) then
          fq(k) = g*rho(i,k)* &
            min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**br),9.1_r8*rhof(i,k))
          fn(k) = g*rho(i,k)* &
            min(an(i,k)*gamma_b_plus1/lam(k)**br,9.1_r8*rhof(i,k))
          ! if fq wasn't capped, compute q*dfq/dq and q*dfq/dn
          if (fq(k) < g*rho(i,k)*9.1_r8*rhof(i,k)) then
            qdfq(:) = g*rho(i,k)* &
              an(i,k)*(-br)*gamma_b_plus4/(6._r8*lam(k)**(br+1._r8))*qic_dlam(:,k)
          end if
          ! if fn wasn't capped, compute n*dfn/dq and n*dfn/dn
          if (fn(k) < g*rho(i,k)*9.1_r8*rhof(i,k)) then
            ndfn(:) = g*rho(i,k)* &
              an(i,k)*(-br)*gamma_b_plus1/lam(k)**(br+1._r8)*nic_dlam(:,k)
          end if
        end if

      else if (mg_type == MG_SNOW) then
        if (lam(k) .ge. qsmall) then
           fq(k) = g*rho(i,k)* &
              min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bs),1.2_r8*rhof(i,k))
           fn(k) = g*rho(i,k)* &
              min(an(i,k)*gamma_b_plus1/lam(k)**bs,1.2_r8*rhof(i,k))
        end if
      end if

      ! original approach assumes Jacobian(F)(q,n,p) = [fq, 0; 0, fn], giving eigenvalues
      ! of fq and fn for propagation speeds... however this is incorrect.
      ! Noting that Jacobian(F)(q,n,p) = [fq + q*dfq/dq, q*dfq/dn; n*dfn/dq, fn + n*dfn/dn],
      ! the eigenvalues are computed numerically.  The dependence of F explicitly
      ! on p can be ignored.

      ! obtain eigenvalues of F_q for propagation speeds
      tr = fq(k) + qdfq(1) + fn(k) + ndfn(2)
      det = (fq(k) + qdfq(1))*(fn(k) + ndfn(2)) - qdfq(2)*ndfn(1)
      s1 = 0.5_r8*(tr + sqrt(tr*tr - 4._r8*det))
      s2 = 0.5_r8*(tr - sqrt(tr*tr - 4._r8*det))

      ! Update CFL number
      cfl = max(cfl,s1*deltat/pdel(i,k),s2*deltat/pdel(i,k))
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
        fluxQ(k) = fq(k)*q(i,k)! + 0.5_r8*pdel(i,k)*dfluxm
      else
        fluxQ(k) = fq(k)*q(i,k)! + 0.5_r8*pdel(i,k)*dfluxp
      end if
      dfluxp = (fn(k+1)*n(i,k+1)-fn(k)*n(i,k))/(0.5_r8*(pdel(i,k+1)+pdel(i,k)))
      dfluxm = (fn(k)*n(i,k)-fn(k-1)*n(i,k-1))/(0.5_r8*(pdel(i,k)+pdel(i,k-1)))
      if (abs(dfluxm) < abs(dfluxp)) then
        fluxN(k) = fn(k)*n(i,k)! + 0.5_r8*pdel(i,k)*dfluxm
      else
        fluxN(k) = fn(k)*n(i,k)! + 0.5_r8*pdel(i,k)*dfluxp
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
