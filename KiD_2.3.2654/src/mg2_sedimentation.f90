module Sedimentation
  use micro_mg2_acme_v1beta_utils, only: r8
  implicit none
  private
  save
  public :: sed_CalcCFL, sed_CalcFlux, sed_AdvanceOneStep
  integer, parameter, public  :: MG_LIQUID = 1
  integer, parameter, public  :: MG_ICE = 2
  integer, parameter, public  :: MG_RAIN = 3
  integer, parameter, public  :: MG_SNOW = 4
  real(r8), parameter, public :: CFL = 0.9_r8

contains

  subroutine sed_CalcCFL(q,n,cloud_frac,rho,pdel,nlev,i, &
    mg_type,deltat,g,an,rhof,computed,alphaq,alphan,s1,s2,w1,w2,cfl,ncons,nnst,gamma_b_plus1,gamma_b_plus4)
    use gptl,                        only: gptlstart_handle, gptlstop_handle
    use gptl_kid
    use micro_mg2_acme_v1beta_utils, only: size_dist_param_liq, &
                                           size_dist_param_basic, &
                                           mg_ice_props, mg_liq_props, &
                                           mg_snow_props, mg_rain_props, &
                                           qsmall, bi, bc, br, bs, &
                                           limiter_is_on
    implicit none
    real(r8), intent(in)            :: q(:,:), n(:,:)
    real(r8), intent(in)            :: rho(:,:), pdel(:,:), an(:,:), rhof(:,:)
    real(r8), intent(in)            :: cloud_frac(:,:), deltat, g
    integer, intent(in)             :: nlev, i, mg_type
    logical, intent(inout)          :: computed(:)
    real(r8), intent(in), optional  :: nnst, gamma_b_plus1, gamma_b_plus4
    logical, intent(in), optional   :: ncons
    real(r8), intent(inout)         :: alphaq(:), alphan(:)
    real(r8), intent(inout)         :: s1(:), s2(:), w1(:,:), w2(:,:)
    real(r8), intent(out)           :: cfl

    real(r8) :: qic(nlev), nic(nlev), lam(nlev), pgam(nlev), lambr_bounds(2)
    real(r8) :: cq, cn, clam, tr, det, sqrtdisc, lambr(nlev)
    integer :: k, ngptl

    ! use quantity in cloud
    qic = q(i,:)/cloud_frac(i,:)
    nic = n(i,:)/cloud_frac(i,:)

    ! switch for specification of droplet and crystal number or ice cloud number
    if ((mg_type == MG_LIQUID  .or. mg_type == MG_ICE) .and. ncons) then
      nic = nnst/rho(i,:)
    end if

    ! calculate lambda
    select case (mg_type)

      case (MG_ICE)
        gptl_ret = gptlstart_handle('Lambda Calculation (ice)', gptl_lambda_ice)
        do ngptl=1,gptl_loop
          do k=1,nlev
            if (.not.computed(k)) call size_dist_param_basic(mg_ice_props, qic(k), nic(k), lam(k))
          end do
        end do
        gptl_ret = gptlstop_handle('Lambda Calculation (ice)', gptl_lambda_ice)

      case (MG_LIQUID)
        gptl_ret = gptlstart_handle('Lambda Calculation (cloud)', gptl_lambda_cloud)
        do ngptl=1,gptl_loop
          do k=1,nlev
            if (.not.computed(k)) call size_dist_param_liq(mg_liq_props, qic(k), nic(k), rho(i,k), pgam(k), lam(k))
          end do
        end do
        gptl_ret = gptlstop_handle('Lambda Calculation (cloud)', gptl_lambda_cloud)

      case (MG_RAIN)
        gptl_ret = gptlstart_handle('Lambda Calculation (rain)', gptl_lambda_rain)
        do ngptl=1,gptl_loop
#ifdef SED_COMBINELAMBDA
          lambr_bounds(1) = mg_rain_props%lambda_bounds(1)**br
          lambr_bounds(2) = mg_rain_props%lambda_bounds(2)**br
          do k=1,nlev
            if (qic(k) > qsmall .and. (.not.computed(k))) then
              ! add upper limit to in-cloud number concentration to prevent
              ! numerical error
              if (limiter_is_on(mg_rain_props%min_mean_mass)) then
                nic(k) = min(nic(k), qic(k) / mg_rain_props%min_mean_mass)
              end if
              ! lambda^b = (c nic/qic)^(b/d)
              lambr(k) = (mg_rain_props%shape_coef * nic(k)/qic(k))**(br/mg_rain_props%eff_dim)
              ! check for slope
              ! adjust vars
              if (lambr(k) < lambr_bounds(1)) then
                lambr(k) = lambr_bounds(1)
                nic(k) = lambr(k)**(mg_rain_props%eff_dim/br) * qic(k)/mg_rain_props%shape_coef
              else if (lambr(k) > lambr_bounds(2)) then
                lambr(k) = lambr_bounds(2)
                nic(k) = lambr(k)**(mg_rain_props%eff_dim/br) * qic(k)/mg_rain_props%shape_coef
              end if
            end if
          end do
#else
        do k=1,nlev
          if (.not.computed(k)) call size_dist_param_basic(mg_rain_props, qic(k), nic(k), lam(k))
        end do
#endif
        end do
        gptl_ret = gptlstop_handle('Lambda Calculation (rain)', gptl_lambda_rain)

      case (MG_SNOW)
        gptl_ret = gptlstart_handle('Lambda Calculation (snow)', gptl_lambda_snow)
        do ngptl=1,gptl_loop
          do k=1,nlev
            if (.not.computed(k)) call size_dist_param_basic(mg_snow_props, qic(k), nic(k), lam(k))
          end do
        end do
        gptl_ret = gptlstop_handle('Lambda Calculation (snow)', gptl_lambda_snow)

      case default
        print *, "Invalid mg_type to mg2_sedimentation"
        stop

    end select


#ifndef SED_COMPFLAG
    alphaq(:) = 0._r8
    alphan(:) = 0._r8
    s1(:) = 0._r8
    s2(:) = 0._r8
#endif

    ! Loop through levels to compute alphaq, alphan, and possibly eigenspaces
    select case (mg_type)

      case (MG_ICE)
        gptl_ret = gptlstart_handle('Fall Speed Calculation (ice)', gptl_fall_ice)
        do ngptl=1,gptl_loop
          do k=1,nlev
            if (qic(k) .ge. qsmall) then
              alphaq(k) = g*rho(i,k)*min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bi), &
                    1.2_r8*rhof(i,k))
              alphan(k) = g*rho(i,k)* &
                    min(an(i,k)*gamma_b_plus1/lam(k)**bi,1.2_r8*rhof(i,k))
              s1(k) = alphaq(k)
              s2(k) = alphan(k)
#ifdef SED_COMPFLAG
              computed(k) = .true.
#endif
#ifdef SED_USEWPA
              ! obtain eigenvectors
              w1(1,k) = 1._r8
              w1(2,k) = 0._r8
              w2(1,k) = 0._r8
              w2(2,k) = 1._r8
#endif
            end if
          end do
        end do
        gptl_ret = gptlstop_handle('Fall Speed Calculation (ice)', gptl_fall_ice)

      case (MG_LIQUID)
        gptl_ret = gptlstart_handle('Fall Speed Calculation (cloud)', gptl_fall_cloud)
        do ngptl=1,gptl_loop
          do k=1,nlev
            if (qic(k) .ge. qsmall) then
              alphaq(k) = g*rho(i,k)*an(i,k)*gamma(4._r8+bc+pgam(k))/ &
                          (lam(k)**bc*gamma(pgam(k)+4._r8))

              alphan(k) = g*rho(i,k)* &
                          an(i,k)*gamma(1._r8+bc+pgam(k))/ &
                          (lam(k)**bc*gamma(pgam(k)+1._r8))
              s1(k) = alphaq(k)
              s2(k) = alphan(k)
#ifdef SED_COMPFLAG
              computed(k) = .true.
#endif
#ifdef SED_USEWPA
              ! obtain eigenvectors
              w1(1,k) = 1._r8
              w1(2,k) = 0._r8
              w2(1,k) = 0._r8
              w2(2,k) = 1._r8
#endif
            end if
          end do
        end do
        gptl_ret = gptlstop_handle('Fall Speed Calculation (cloud)', gptl_fall_cloud)

      case (MG_RAIN)
        gptl_ret = gptlstart_handle('Fall Speed Calculation (rain)', gptl_fall_rain)
        do ngptl=1,gptl_loop
          do k=1,nlev
            if (qic(k) > qsmall .and. (.not.computed(k))) then
              cq = g*rho(i,k)*an(i,k)*gamma_b_plus4/6._r8
              cn = g*rho(i,k)*an(i,k)*gamma_b_plus1
#ifdef SED_COMBINELAMBDA
              alphaq(k) = min(cq/lambr(k), 9.1_r8*g*rho(i,k)*rhof(i,k))
              alphan(k) = min(cn/lambr(k), 9.1_r8*g*rho(i,k)*rhof(i,k))
#else
              alphaq(k) = min(cq/lam(k)**br, 9.1_r8*g*rho(i,k)*rhof(i,k))
              alphan(k) = min(cn/lam(k)**br, 9.1_r8*g*rho(i,k)*rhof(i,k))
#endif

#ifdef SED_NONLINEAR
              ! obtain eigenvalues considering
              ! Jacobian([alphaq(q,n)*q, alphan(q,n)*n]^T) =
              ! [alphaq + q*dalphaq/dq,          q*dalphaq/dn
              !           n*dalphan/dq, alphan + n*dalphan/dn]
              tr = cq*(1._r8 + br/mg_rain_props%eff_dim) + &
                   cn*(1._r8 - br/mg_rain_props%eff_dim)
              det = cq*cn
              sqrtdisc = sqrt(tr*tr - 4._r8*det)
              s1(k) = 0.5_r8/lambr(k) * (tr - sqrtdisc)
              s2(k) = 0.5_r8/lambr(k) * (tr + sqrtdisc)
#else
              ! obtain eigenvalues considering
              ! Jacobian([alphaq*q, alphan*n]^T) =
              ! [alphaq, 0.0
              !     0.0, alphan]
              s1(k) = alphaq(k)
              s2(k) = alphan(k)
#endif
#ifdef SED_COMPFLAG
              computed(k) = .true.
#endif
#ifdef SED_USEWPA
              ! obtain eigenvectors
#ifdef SED_NONLINEAR
              w1(1,k) = mg_rain_props%shape_coef*(br/mg_rain_props%eff_dim) / &
                         lambr(k)**(mg_rain_props%eff_dim/br)
              w1(2,k) = 0.5_r8/cq * &
                        ( cq*(1._r8 + br/mg_rain_props%eff_dim) - &
                        cn*(1._r8 - br/mg_rain_props%eff_dim) + sqrtdisc )
              w2(1,k) = mg_rain_props%shape_coef*(br/mg_rain_props%eff_dim) / &
                         lambr(k)**(mg_rain_props%eff_dim/br)
              w2(2,k) = 0.5_r8/cq * &
                        ( cq*(1._r8 + br/mg_rain_props%eff_dim) - &
                        cn*(1._r8 - br/mg_rain_props%eff_dim) - sqrtdisc )
#else
              w1(1,k) = 1._r8
              w1(2,k) = 0._r8
              w2(1,k) = 0._r8
              w2(2,k) = 1._r8
#endif
#endif
            end if
          end do
        end do
        gptl_ret = gptlstop_handle('Fall Speed Calculation (rain)', gptl_fall_rain)

      case (MG_SNOW)
        gptl_ret = gptlstart_handle('Fall Speed Calculation (snow)', gptl_fall_snow)
        do ngptl=1,gptl_loop
          do k=1,nlev
            if (lam(k) .ge. qsmall) then
              alphaq(k) = g*rho(i,k)* &
                min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bs),1.2_r8*rhof(i,k))
              alphan(k) = g*rho(i,k)* &
                min(an(i,k)*gamma_b_plus1/lam(k)**bs,1.2_r8*rhof(i,k))
              s1(k) = alphaq(k)
              s2(k) = alphan(k)
#ifdef SED_COMPFLAG
              computed(k) = .true.
#endif
#ifdef SED_USEWPA
              ! obtain eigenvectors
              w1(1,k) = 1._r8
              w1(2,k) = 0._r8
              w2(1,k) = 0._r8
              w2(2,k) = 1._r8
#endif
            end if
          end do
        end do
        gptl_ret = gptlstop_handle('Fall Speed Calculation (snow)', gptl_fall_snow)

      case default
        print *, "Invalid mg_type to mg2_sedimentation"
        stop

    end select

    ! compute CFL number
    cfl = max(maxval(s1(:)*deltat/pdel(i,:)),maxval(s2(:)*deltat/pdel(i,:)))

  end subroutine sed_CalcCFL

  subroutine sed_CalcFlux(q,n,alphaq,alphan,pdel,deltat,s1,s2,w1,w2, &
                          nlev,i,mg_type,fq,fn)
    use micro_mg2_acme_v1beta_utils, only: qsmall
    implicit none
    real(r8), intent(in)            :: q(:,:), n(:,:)
    real(r8), intent(in)            :: pdel(:,:), deltat
    integer,  intent(in)            :: nlev, i, mg_type
    real(r8), intent(out)           :: fq(:,:), fn(:,:)
    real(r8), intent(inout)         :: s1(:), s2(:), w1(:,:), w2(:,:)
    real(r8), intent(inout)         :: alphaq(:), alphan(:)
    ! note that "inout" variables here are really just "in" variables that are
    ! modified
    real(r8) :: qtemp(0:nlev), ntemp(0:nlev)
    real(r8) :: a, b, c, d, dq, dn, coeff1(nlev), coeff2(nlev), sig1(2), sig2(2)
    real(r8) :: mult1, mult2, phi1, phi2, theta
    integer :: k

    ! initialize values to zero
    fq = 0._r8
    fn = 0._r8

    ! initialize qtemp and ntemp with ghost cell
    qtemp(0) = 0._r8
    qtemp(1:nlev) = q(i,:)
    ntemp(0) = 0._r8
    ntemp(1:nlev) = n(i,:)

    ! Compute flux
#ifdef SED_USEWPA
    do k=1,nlev
      if (s1(k) < qsmall .and. k > 1 .and. s1(k-1) > qsmall) then
        ! to handle waves going from a wet cell into a dry cell,
        ! copy eigenvectors and eigenvalues of Jacobian(F(q,n))
        ! from wet x_{k-1} cell into dry x_k cell
        s1(k) = s1(k-1)
        w1(:,k) = w1(:,k-1)
        s2(k) = s2(k-1)
        w2(:,k) = w2(:,k-1)
        alphaq(k) = alphaq(k-1)
        alphan(k) = alphan(k-1)
      end if

      if (s1(k) > qsmall) then
      ! scaled jump in x_{k-1} and x_k values will be added to flux at
      ! x_{k-1/2} using eigenvectors and eigenvalues of Jacobian(F(q,b))
      ! in x_k cell
        dq = alphaq(k-1)/alphaq(k)*qtemp(k-1) - qtemp(k)
        dn = alphan(k-1)/alphan(k)*ntemp(k-1) - ntemp(k)
        ! solve for coeff1 and coeff2 in [dq, dn]^T = coeff1*w1 + coeff2*w2
        a = w1(1,k)
        c = w1(2,k)
        b = w2(1,k)
        d = w2(2,k)
        if (abs(a*d - b*c) < 1.d-12) then
          print *, a, b, c, d
          stop "matrix nearly singular"
        end if
        coeff1(k) = (d*dq - b*dn)/(a*d - b*c)
        coeff2(k) = (a*dn - c*dq)/(a*d - b*c)
        ! Add first order terms to flux
        fq(k-1,2) = s1(k)*coeff1(k)*w1(1,k) + s2(k)*coeff2(k)*w2(1,k)
        fn(k-1,2) = s1(k)*coeff1(k)*w1(2,k) + s2(k)*coeff2(k)*w2(2,k)
        ! Need flux limiter value for second order term
        phi1 = 1._r8 ! phi ~ 1 for smooth data
        phi2 = 1._r8
        if (k > 1) then
          theta = coeff1(k-1)/coeff1(k)
          if (abs(theta) < 1._r8) phi1 = theta ! minmod
          !phi1 = max(0.r_8,min(1._r8,2._r8*theta),min(2._r8,theta)) ! superbee
          !phi1 = max(0._r8,min(0.5_r8*(1._r8+theta),2._r8,2._r8*theta)) ! MC
          !phi1 = (theta + abs(theta))/(1._r8 + abs(theta)) ! van Leer
          theta = coeff2(k-1)/coeff2(k)
          if (abs(theta) < 1._r8) phi2 = theta ! minmod
          !phi2 = max(0.r_8,min(1._r8,2._r8*theta),min(2._r8,theta)) ! superbee
          !phi2 = max(0._r8,min(0.5_r8*(1._r8+theta),2._r8,2._r8*theta)) ! MC
          !phi2 = (theta + abs(theta))/(1._r8 + abs(theta)) ! van Leer
        end if
        ! Add second order terms to flux
        mult1 = 0.5_r8*s1(k)*(1._r8-s1(k)*deltat/pdel(i,k-1))
        mult2 = 0.5_r8*s2(k)*(1._r8-s2(k)*deltat/pdel(i,k-1))
        fq(k-1,:) = fq(k-1,:) + mult1*phi1*coeff1(k)*w1(1,k) &
                              + mult2*phi2*coeff2(k)*w2(1,k)
        fn(k-1,:) = fn(k-1,:) + mult1*phi1*coeff1(k)*w1(2,k) &
                              + mult2*phi2*coeff2(k)*w2(2,k)
      end if

      if (k == nlev) then
        ! store flux out for surface precipitation calculation
        fq(nlev,2) = alphaq(nlev)*qtemp(nlev)
      end if
    end do
#else
    fq(:,1) = alphaq(:)*qtemp(:)
    fq(:,2) = fq(:,1)
    fn(:,1) = alphan(:)*ntemp(:)
    fn(:,2) = fn(:,1)
#endif

  end subroutine sed_CalcFlux


  subroutine sed_AdvanceOneStep(q,fq,n,fn,pdel,deltat,deltat_sed,nlev,i,mg_type,g, &
    qtend,ntend,prect,qsedtend,cloud_frac,qvlat,tlat,xxl,preci,qsevap)
    use micro_mg2_acme_v1beta_utils, only: qsmall
    implicit none
    real(r8), intent(in)           :: pdel(:,:)
    real(r8), intent(in)              :: deltat, deltat_sed, g
    integer, intent(in)               :: nlev, i, mg_type
    real(r8), intent(inout)           :: q(:,:), qtend(:,:), n(:,:), ntend(:,:)
    real(r8), intent(inout)           :: fq(:,:), fn(:,:), prect(:)
    real(r8), intent(inout)           :: qsedtend(:,:)
    real(r8), intent(in), optional    :: cloud_frac(:,:), xxl
    real(r8), intent(inout), optional :: qvlat(:,:), tlat(:,:), preci(:)
    real(r8), intent(inout), optional :: qsevap(:,:)

    real(r8) :: ratio(nlev), deltaFluxQ, deltaFluxN, deltafluxQ_evap
    integer :: k

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
      deltafluxQ = (fq(k,1)-ratio(k)*fq(k-1,2))/pdel(i,k)
      deltafluxN = (fn(k,1)-ratio(k)*fn(k-1,2))/pdel(i,k)
      ! add fallout terms to eulerian tendencies
      qtend(i,k) = qtend(i,k) - deltafluxQ*deltat_sed/deltat
      ntend(i,k) = ntend(i,k) - deltafluxN*deltat_sed/deltat
      ! sedimentation tendency for output
      qsedtend(i,k) = qsedtend(i,k) - deltafluxQ*deltat_sed/deltat
      if (mg_type == MG_ICE .or. mg_type == MG_LIQUID) then
        ! add terms to to evap/sub of cloud water
        deltafluxQ_evap = (ratio(k)-1._r8)*fq(k-1,2)/pdel(i,k)
        qvlat(i,k) = qvlat(i,k) - deltafluxQ_evap*deltat_sed/deltat
        qsevap(i,k) = qsevap(i,k) - deltafluxQ_evap*deltat_sed/deltat
        tlat(i,k) = tlat(i,k) + deltafluxQ_evap*xxl*deltat_sed/deltat
      end if

      q(i,k) = q(i,k) - deltat_sed*deltafluxQ
      n(i,k) = n(i,k) - deltat_sed*deltafluxN

      if (q(i,k) < -1.d-10) then
        print *, q(i,k)
        stop "q negative"
      else if (n(i,k) < -1.d-10) then
        print *, n(i,k)
        stop "n negative"
      end if
    end do

    ! units below are m/s
    ! sedimentation flux at surface is added to precip flux at surface
    ! to get total precip (cloud + precip water) rate

    prect(i) = prect(i) + deltat_sed/deltat*fq(nlev,2)/g/1000._r8
    if (mg_type == MG_ICE .or. mg_type == MG_SNOW) then
      preci(i) = preci(i) + deltat_sed/deltat*fq(nlev,2)/g/1000._r8
    end if

  end subroutine

end module Sedimentation
