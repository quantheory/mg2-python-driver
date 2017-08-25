module Sedimentation
  use micro_mg2_acme_v1beta_utils, only: r8
  implicit none
  private
  save
  public :: sed_CalcFluxAndCFL, sed_AdvanceOneStep
  integer, parameter, public  :: MG_LIQUID = 1
  integer, parameter, public  :: MG_ICE = 2
  integer, parameter, public  :: MG_RAIN = 3
  integer, parameter, public  :: MG_SNOW = 4
  real(r8), parameter, public :: CFL = 0.9_r8

contains

  subroutine sed_CalcFluxAndCFL(q,qtend,n,ntend,cloud_frac,rho,pdel,nlev,i, &
    mg_type,deltat,g,an,rhof,fq,fn,cfl,ncons,nnst,gamma_b_plus1,gamma_b_plus4)
    use micro_mg2_acme_v1beta_utils, only: size_dist_param_liq, &
                                           size_dist_param_basic, &
                                           mg_ice_props, mg_liq_props, &
                                           mg_snow_props, mg_rain_props, &
                                           qsmall, bi, bc, br, bs, &
                                           limiter_is_on
    implicit none
    real(r8), intent(in)            :: q(:,:), qtend(:,:), n(:,:), ntend(:,:)
    real(r8), intent(in)            :: rho(:,:), pdel(:,:), an(:,:), rhof(:,:)
    real(r8), intent(in)            :: cloud_frac(:,:), deltat, g
    integer,  intent(in)            :: nlev, i, mg_type
    real(r8), intent(in), optional  :: nnst, gamma_b_plus1, gamma_b_plus4
    logical, intent(in), optional   :: ncons
    real(r8), intent(out)           :: cfl, fq(:,:), fn(:,:)

    real(r8) :: qic(nlev), nic(nlev), qtemp(0:nlev), ntemp(0:nlev)
    real(r8) :: alphaq(0:nlev), alphan(0:nlev), lambda_bounds(2)
    real(r8) :: lam(nlev), pgam(nlev), cq, cn, clam, tr, sqrtdisc, lambr(nlev)
    real(r8) :: s1(nlev), s2(nlev), a, b, c, d, dq, dn, coeff(2)
    real(r8) :: w1(2,nlev), w2(2,nlev)
    integer :: k
    logical :: limited(nlev)

    ! initialize values to zero
    cfl = 0._r8
    fq = 0._r8
    fn = 0._r8
    alphaq = 0._r8
    alphan = 0._r8
    s1 = 0._r8
    s2 = 0._r8

    ! initialize qtemp and ntemp with ghost cell
    qtemp(0) = 0._r8
    qtemp(1:nlev) = q(i,:)
    ntemp(0) = 0._r8
    ntemp(1:nlev) = n(i,:)

    ! use quantity in cloud
    qic = q(i,:)/cloud_frac(i,:)
    nic = n(i,:)/cloud_frac(i,:)

    ! switch for specification of droplet and crystal number or ice cloud number
    if ((mg_type == MG_LIQUID  .or. mg_type == MG_ICE) .and. ncons) then
      nic = nnst/rho(i,:)
    end if

    ! compute lam and pgam using elemental function call
    if (mg_type == MG_LIQUID) then
      call size_dist_param_liq(mg_liq_props, qic, nic, rho(i,:), &
                               pgam, lam)

    else if (mg_type == MG_ICE) then
      call size_dist_param_basic(mg_ice_props, qic, nic, lam)

    else if (mg_type == MG_RAIN) then
      lambda_bounds(1) = mg_rain_props%lambda_bounds(1)**br
      lambda_bounds(2) = mg_rain_props%lambda_bounds(2)**br
      do k=1,nlev
        if (qic(k) > qsmall) then
          ! add upper limit to in-cloud number concentration to prevent
          ! numerical error
          if (limiter_is_on(mg_rain_props%min_mean_mass)) then
            nic(k) = min(nic(k), qic(k) / mg_rain_props%min_mean_mass)
          end if
          ! lambda^b = (c nic/qic)^(b/d)
          lambr(k) = (mg_rain_props%shape_coef * nic(k)/qic(k))**(br/mg_rain_props%eff_dim)
          ! check for slope
          ! adjust vars
          limited(k) = .false.
          if (lambr(k) < lambda_bounds(1)) then
            lambr(k) = lambda_bounds(1)
            nic(k) = lambr(k)**(mg_rain_props%eff_dim/br) * qic(k)/mg_rain_props%shape_coef
            limited(k) = .true.
          else if (lambr(k) > lambda_bounds(2)) then
            lambr(k) = lambda_bounds(2)
            nic(k) = lambr(k)**(mg_rain_props%eff_dim/br) * qic(k)/mg_rain_props%shape_coef
            limited(k) = .true.
          end if
        end if
      end do

    else if (mg_type == MG_SNOW) then
       call size_dist_param_basic(mg_snow_props, qic, nic, lam)
    else
       print *, "Invalid mg_type to mg2_sedimentation"
       stop
    end if

    ! Loop through levels to compute alphaq, alphan, and the CFL number
    do k=1,nlev

      if(mg_type == MG_LIQUID) then
        if (qic(k) .ge. qsmall) then
          alphaq(k) = g*rho(i,k)*an(i,k)*gamma(4._r8+bc+pgam(k))/ &
               (lam(k)**bc*gamma(pgam(k)+4._r8))

          alphan(k) = g*rho(i,k)* &
               an(i,k)*gamma(1._r8+bc+pgam(k))/ &
               (lam(k)**bc*gamma(pgam(k)+1._r8))
          s1(k) = alphaq(k)
          s2(k) = alphan(k)
          ! update fluxes
          fq(k,:) = s1(k)*q(i,k)
          fn(k,:) = s2(k)*n(i,k)
        end if

      else if (mg_type == MG_ICE) then
        if (qic(k) .ge. qsmall) then
          alphaq(k) = g*rho(i,k)*min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bi), &
                1.2_r8*rhof(i,k))
          alphan(k) = g*rho(i,k)* &
                min(an(i,k)*gamma_b_plus1/lam(k)**bi,1.2_r8*rhof(i,k))
          s1(k) = alphaq(k)
          s2(k) = alphan(k)
          ! update fluxes
          fq(k,:) = s1(k)*q(i,k)
          fn(k,:) = s2(k)*n(i,k)
        end if

      else if (mg_type == MG_RAIN) then

        ! use wave propagation algorithm to determine flux modifiers
        if (qic(k) > qsmall) then
        ! scaled jump in x_{k-1} and x_k values will be added to flux at
        ! x_{k-1/2} using eigenvectors and eigenvalues of Jacobian(F(q,b))
        ! in x_k cell
          cq = g*rho(i,k)*an(i,k)*gamma_b_plus4/6._r8
          cn = g*rho(i,k)*an(i,k)*gamma_b_plus1
          alphaq(k) = min(cq/lambr(k), 9.1_r8*g*rho(i,k)*rhof(i,k))
          alphan(k) = min(cn/lambr(k), 9.1_r8*g*rho(i,k)*rhof(i,k))
          ! compute scaled jump
          dq = alphaq(k-1)/alphaq(k)*qtemp(k-1) - qtemp(k)
          dn = alphan(k-1)/alphan(k)*ntemp(k-1) - ntemp(k)
          ! obtain eigenvectors and eigenvalues
          if (limited(k)) then
          ! Jacobian([alphaq*q, alphan*n]^T) =
          ! [alphaq,      0
          !       0, alphan]
            s1(k) = alphaq(k)
            w1(1,k) = 1._r8
            w1(2,k) = 0._r8
            s2(k) = alphan(k)
            w2(1,k) = 0._r8
            w2(2,k) = 1._r8
          else
          ! Jacobian([alphaq(q,n)*q, alphan(q,n)*n]^T) =
          ! [alphaq + q*dalphaq/dq,          q*dalphaq/dn
          !           n*dalphan/dq, alphan + n*dalphan/dn]
            tr = cq*(1._r8 - br/mg_rain_props%eff_dim) + &
                 cn*(1._r8 + br/mg_rain_props%eff_dim)
            sqrtdisc = sqrt( (1._r8 + br/mg_rain_props%eff_dim)**2 - &
                             4._r8*br/mg_rain_props%eff_dim * cq/(cq-cn) )
            s1(k) = 0.5_r8/lambr(k) * (tr - (cq-cn)*sqrtdisc)
            w1(1,k) = qic(k)/nic(k)
            w1(2,k) = -0.5_r8*mg_rain_props%eff_dim/br/cq * &
                      ( cq*(1._r8 - br/mg_rain_props%eff_dim) - &
                      cn*(1._r8 + br/mg_rain_props%eff_dim) + (cq-cn)*sqrtdisc )
            s2(k) = 0.5_r8/lambr(k) * (tr + (cq-cn)*sqrtdisc)
            w2(1,k) = qic(k)/nic(k)
            w2(2,k) = -0.5_r8*mg_rain_props%eff_dim/br/cq * &
                      ( cq*(1._r8 - br/mg_rain_props%eff_dim) - &
                      cn*(1._r8 + br/mg_rain_props%eff_dim) - (cq-cn)*sqrtdisc )
          end if

        else if (k > 1 .and. qic(k-1) > qsmall) then
        ! to handle waves going from a non-zero cell into a zero cell,
        ! jump in x_{k-1} and x_k values will be added to flux at
        ! x_{k-1/2} using eigenvectors and eigenvalues of Jacobian(F(q,b))
        ! in x_{k-1} cell
          dq = qtemp(k-1) - qtemp(k)
          dn = ntemp(k-1) - ntemp(k)
          s1(k) = s1(k-1)
          w1(:,k) = w1(:,k-1)
          s2(k) = s2(k-1)
          w2(:,k) = w2(:,k-1)
        end if

        if (s1(k) > 0._r8) then
        ! Update fluxes at x_{k-1/2} to include the waves coeff1*w1 and coeff2*w2
          ! solve for coeff1, coeff2 in [dq, dn]^T = coeff1*w1 + coeff2*w2
          a = w1(1,k)
          c = w1(2,k)
          b = w2(1,k)
          d = w2(2,k)
          if (abs(a*d - b*c) < 1.d-12) then
            print *, lambr
            print *, q(i,k), n(i,k)
            print *, a, b, c, d
            stop "matrix nearly singular"
          end if
          coeff(1) = (d*dq - b*dn)/(a*d - b*c)
          coeff(2) = (a*dn - c*dq)/(a*d - b*c)

          fq(k-1,2) = s1(k)*coeff(1)*w1(1,k) + s2(k)*coeff(2)*w2(1,k)
          fn(k-1,2) = s1(k)*coeff(1)*w1(2,k) + s2(k)*coeff(2)*w2(2,k)
        end if

        if (k == nlev) then
          ! store flux out for surface precipitation calculation
          fq(nlev,2) = alphaq(nlev)*q(i,nlev)
        end if

      else if (mg_type == MG_SNOW) then
        if (lam(k) .ge. qsmall) then
           alphaq(k) = g*rho(i,k)* &
              min(an(i,k)*gamma_b_plus4/(6._r8*lam(k)**bs),1.2_r8*rhof(i,k))
           alphan(k) = g*rho(i,k)* &
              min(an(i,k)*gamma_b_plus1/lam(k)**bs,1.2_r8*rhof(i,k))
           s1(k) = alphaq(k)
           s2(k) = alphan(k)
           ! update fluxes
           fq(k,:) = s1(k)*q(i,k)
           fn(k,:) = s2(k)*n(i,k)
        end if
      end if

      ! Update CFL number
      cfl = max(cfl,s1(k)*deltat/pdel(i,k),s2(k)*deltat/pdel(i,k))

      ! TODO: detect shock formation

    end do

  end subroutine sed_CalcFluxAndCFL

  subroutine sed_AdvanceOneStep(q,fq,n,fn,pdel,deltat,deltat_sed,nlev,i,mg_type,g, &
    qtend,ntend,prect,qsedtend,cloud_frac,qvlat,tlat,xxl,preci,qsevap)
    implicit none
    real(r8), intent(in)           :: pdel(:,:)
    real(r8), intent(in)              :: deltat, deltat_sed, g
    integer, intent(in)               :: nlev, i, mg_type
    real(r8), intent(inout)           :: q(:,:), qtend(:,:), n(:,:), ntend(:,:)
    real(r8), intent(inout)           :: fq(:,:), fn(:,:), prect(:)
    real(r8), intent(in), optional    :: cloud_frac(:,:), xxl
    real(r8), intent(inout), optional :: qvlat(:,:), tlat(:,:), preci(:)
    real(r8), intent(out), optional   :: qsevap(:,:)
    real(r8), intent(out)             :: qsedtend(:,:)

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

      if (q(i,k) < -1.d-16) then
        print *, q(i,k)
        stop "q negative"
      else if (n(i,k) < -1.d-16) then
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
