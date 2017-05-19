module mphys_stats

  ! This module contains counters for run accumlating run time stats
  !
  ! Note: Array values are set in the microphysics and so are idexed as top down
  ! while kid assumes bottom up ordering in arrays. So need to flip arrays before
  ! outputting with netcdf
  !
  ! All limtier counters are preceded by a comment with: LIMITER
  
  use parameters, only : nx, nz
  use namelists, only  : KiD_outdir, KiD_outfile, fileNameOut
  
  implicit none
  private
  save

  public :: init_stats, finalize_stats
  public :: zero_stats
  public :: write_stats, write_stats_summary
  
  integer :: ncalls = 0
  
  ! limit on max rain mixing ratio (relative change)
  integer,          public, dimension(nx,nz) :: qric_limiter
  double precision, public, dimension(nx,nz) :: qric_limiter_mag
  double precision :: qric_limiter_mag_max = 0.0d0

  ! setting small rain values to zero (absolute change)
  integer,          public, dimension(nx,nz) :: qric_qsmall
  integer,          public, dimension(nx,nz) :: nric_qsmall
  double precision, public, dimension(nx,nz) :: qric_qsmall_mag
  double precision, public, dimension(nx,nz) :: nric_qsmall_mag
  double precision :: qric_qsmall_mag_max = 0.0d0
  double precision :: nric_qsmall_mag_max = 0.0d0

  ! setting negative rain number values to zero (absolute change)
  integer,          public, dimension(nx,nz) :: nric_neg
  double precision, public, dimension(nx,nz) :: nric_neg_mag
  double precision :: nric_neg_mag_max = 0.0d0
  
  ! limit on max snow mixing ratio (relative change)
  integer,          public, dimension(nx,nz) :: qsic_limiter
  double precision, public, dimension(nx,nz) :: qsic_limiter_mag
  double precision :: qsic_limiter_mag_max = 0.0d0

  ! setting small snow values to zero (absolute change)
  integer,          public, dimension(nx,nz) :: qsic_qsmall
  integer,          public, dimension(nx,nz) :: nsic_qsmall
  double precision, public, dimension(nx,nz) :: qsic_qsmall_mag
  double precision, public, dimension(nx,nz) :: nsic_qsmall_mag
  double precision :: qsic_qsmall_mag_max = 0.0d0
  double precision :: nsic_qsmall_mag_max = 0.0d0

  ! setting negative snow number values to zero (absolute change)
  integer,          public, dimension(nx,nz) :: nsic_neg
  double precision, public, dimension(nx,nz) :: nsic_neg_mag
  double precision :: nsic_neg_mag_max = 0.0d0

  ! -----------------------------------------------------------------------------
  ! Begin: Big "adminstration" loop for conservation
  ! -----------------------------------------------------------------------------
  
  ! cloud water mass (tendency scaling factor) 
  integer,          public, dimension(nx,nz) :: qc_conservation
  double precision, public, dimension(nx,nz) :: qc_conservation_mag  
  double precision :: qc_conservation_mag_max = 0.0d0

  ! mnuccd tendency (absolute change)
  integer,          public, dimension(nx,nz) :: ice_nucleation_limiter
  double precision, public, dimension(nx,nz) :: ice_nucleation_limiter_mag
  double precision :: ice_nucleation_limiter_mag_max = 0.0d0

  ! vap_dep tendency (absolute change)
  integer,          public, dimension(nx,nz) :: ice_deposition_limiter
  double precision, public, dimension(nx,nz) :: ice_deposition_limiter_mag
  double precision :: ice_deposition_limiter_mag_max = 0.0d0

  ! cloud water number (tendency scaling factor)
  integer,          public, dimension(nx,nz) :: nc_conservation
  double precision, public, dimension(nx,nz) :: nc_conservation_mag  
  double precision :: nc_conservation_mag_max = 0.0d0

  ! rain mass (tendency scaling factor)
  integer,          public, dimension(nx,nz) :: qr_conservation
  double precision, public, dimension(nx,nz) :: qr_conservation_mag  
  double precision :: qr_conservation_mag_max = 0.0d0

  ! rain number (tendency scaling factor)
  integer,          public, dimension(nx,nz) :: nr_conservation
  double precision, public, dimension(nx,nz) :: nr_conservation_mag  
  double precision :: nr_conservation_mag_max = 0.0d0

  ! cloud ice mass (tendency scaling factor)
  integer,          public, dimension(nx,nz) :: qi_conservation
  double precision, public, dimension(nx,nz) :: qi_conservation_mag  
  double precision :: qi_conservation_mag_max = 0.0d0

  ! cloud ice number (tendency scaling factor)
  integer,          public, dimension(nx,nz) :: ni_conservation
  double precision, public, dimension(nx,nz) :: ni_conservation_mag  
  double precision :: ni_conservation_mag_max = 0.0d0

  ! snow mass (tendency scaling factor)
  integer,          public, dimension(nx,nz) :: qs_conservation
  double precision, public, dimension(nx,nz) :: qs_conservation_mag  
  double precision :: qs_conservation_mag_max = 0.0d0

  ! snow number (tendency scaling factor)
  integer,          public, dimension(nx,nz) :: ns_conservation
  double precision, public, dimension(nx,nz) :: ns_conservation_mag  
  double precision :: ns_conservation_mag_max = 0.0d0

  ! the next three (rain, snow, and ice) limiters appeared before qc conservation 
  ! in an earlier version of MG2
  integer, public, dimension(nx,nz) :: qiqs_sublimation_qr_evaporation_limiter

  ! rain evaporation (absolute change)
  double precision, public, dimension(nx,nz) :: rain_evaporation_limiter_mag 
  double precision :: rain_evaporation_limiter_mag_max = 0.0d0

  ! snow sublimation (absolute change)
  double precision, public, dimension(nx,nz) :: snow_sublimation_limiter_mag 
  double precision :: snow_sublimation_limiter_mag_max = 0.0d0
  
  ! ice sublimation (absolute change)
  double precision, public, dimension(nx,nz) :: ice_sublimation_limiter_mag  
  double precision :: ice_sublimation_limiter_mag_max = 0.0d0

  ! -----------------------------------------------------------------------------
  ! 1st End: Big "adminstration" loop for conservation
  ! -----------------------------------------------------------------------------
 
  ! cloud ice number (relative change)
  integer,          public, dimension(nx,nz) :: ni_tendency_limiter
  double precision, public, dimension(nx,nz) :: ni_tendency_limiter_mag
  double precision :: ni_tendency_limiter_mag_max = 0.0d0

  ! integer,          public, dimension(nx,nz) :: nr_evap_limiter
  ! double precision, public, dimension(nx,nz) :: nr_evap_limiter_mag
    
  ! -----------------------------------------------------------------------------
  ! 2nd End: Big "adminstration" loop for conservation
  ! -----------------------------------------------------------------------------

  ! saturation adjustmet at end of sedimentaiton has limiters but 
  ! allow_sed_supersat flag is false, so these are not applied

  ! -----------------------------------------------------------------------------
  ! Sedimentation substeps
  ! -----------------------------------------------------------------------------
  integer, public :: nsteps_qi
  integer, public :: nsteps_qc
  integer, public :: nsteps_qr
  integer, public :: nsteps_qs

contains
  
  subroutine init_stats()
    
    character(100) :: fname ! output file name
    integer        :: fid   ! output file id

    fid = 101
    fname = trim(fileNameOut)//'.limiters.txt'
    open(unit = fid, file = fname)
    write(fid,'(A)') "================================================================================"
    write(fid,'(A)') "Limiter Counts"
    close(fid)

    call zero_stats()

  end subroutine init_stats



  subroutine zero_stats()

    ! limiter stats values

    qric_limiter      = 0
    qric_limiter_mag  = 0.0d0

    qric_qsmall     = 0
    nric_qsmall     = 0    
    qric_qsmall_mag = 0.0d0
    nric_qsmall_mag = 0.0d0
    
    nric_neg     = 0
    nric_neg_mag = 0.0d0
    
    qsic_limiter     = 0
    qsic_limiter_mag = 0.0d0
    
    qsic_qsmall     = 0
    nsic_qsmall     = 0
    qsic_qsmall_mag = 0.0d0
    nsic_qsmall_mag = 0.0d0
    
    nsic_neg     = 0
    nsic_neg_mag = 0.0d0
    
    qc_conservation     = 0
    qc_conservation_mag = 0.0d0
    
    ice_nucleation_limiter     = 0
    ice_nucleation_limiter_mag = 0.0d0
    
    ice_deposition_limiter     = 0
    ice_deposition_limiter_mag = 0.0d0
    
    nc_conservation     = 0
    nc_conservation_mag = 0.0d0
    
    qr_conservation     = 0
    qr_conservation_mag = 0.0d0  
    
    nr_conservation     = 0
    nr_conservation_mag = 0.0d0  
    
    qi_conservation     = 0
    qi_conservation_mag = 0.0d0  
    
    ni_conservation     = 0
    ni_conservation_mag = 0.0d0  
    
    qs_conservation     = 0
    qs_conservation_mag = 0.0d0  
    
    ns_conservation     = 0
    ns_conservation_mag = 0.0d0  
    
    qiqs_sublimation_qr_evaporation_limiter = 0
    rain_evaporation_limiter_mag = 0.0d0 
    snow_sublimation_limiter_mag = 0.0d0 
    ice_sublimation_limiter_mag  = 0.0d0  
    
    ni_tendency_limiter     = 0
    ni_tendency_limiter_mag = 0.0d0

    ! sedimentation substeps
    nsteps_qi = 0 
    nsteps_qc = 0 
    nsteps_qr = 0 
    nsteps_qs = 0 

  end subroutine zero_stats



  subroutine write_stats(step)

    integer, intent(in) :: step
    
    character(100) :: fname ! output file name
    integer        :: fid   ! output file id

    integer :: i,j
    character(len=25) :: fmt

    fmt = '(A,I0,A,I0,A,I0,E23.15)'

    fid = 101
    fname = trim(fileNameOut)//'.limiters.txt'
    open(unit = fid, file = fname, status="old", position="append", action="write")

    write(fid,'(A)')      "================================================================================"
    write(fid,'(A,I0,A)') "Step",step," (nx,nz)" 
    write(fid,'(A)')      "================================================================================"

    ! Write output in reverse order to match KiD NetCDF file
    ! MG2 1 = top,  nz = surf 
    ! KiD 1 = surf, nz = top
    do j=nz,1,-1
       do i=1,nx

          write(fid,fmt) "qric_limiter(",i,",",nz-j+1,")                                   ",&
               qric_limiter(i,j), qric_limiter_mag(i,j)

          write(fid,fmt) "qric_qsmall(",i,",",nz-j+1,")                                    ",& 
               qric_qsmall(i,j), qric_qsmall_mag(i,j)

          write(fid,fmt) "nric_qsmall(",i,",",nz-j+1,")                                    ",& 
               nric_qsmall(i,j), nric_qsmall_mag(i,j)

          write(fid,fmt) "nric_neg(",i,",",nz-j+1,")                                       ",& 
               nric_neg(i,j), nric_neg_mag(i,j)

          write(fid,fmt) "qsic_limiter(",i,",",nz-j+1,")                                   ",& 
               qsic_limiter(i,j), qsic_limiter_mag(i,j)

          write(fid,fmt) "qsic_qsmall(",i,",",nz-j+1,")                                    ",& 
               qsic_qsmall(i,j), qsic_qsmall_mag(i,j)

          write(fid,fmt) "nsic_qsmall(",i,",",nz-j+1,")                                    ",& 
               nsic_qsmall(i,j), nsic_qsmall_mag(i,j)

          write(fid,fmt) "nsic_neg(",i,",",nz-j+1,")                                       ",& 
               nsic_neg(i,j), nsic_neg_mag(i,j)

          write(fid,fmt) "qc_conservation(",i,",",nz-j+1,")                                ",& 
               qc_conservation(i,j), qc_conservation_mag(i,j)

          write(fid,fmt) "ice_nucleation_limiter(",i,",",nz-j+1,")                         ",& 
               ice_nucleation_limiter(i,j), ice_nucleation_limiter_mag(i,j)

          write(fid,fmt) "ice_deposition_limiter(",i,",",nz-j+1,")                         ",& 
               ice_deposition_limiter(i,j), ice_deposition_limiter_mag(i,j)

          write(fid,fmt) "nc_conservation(",i,",",nz-j+1,")                                ",& 
               nc_conservation(i,j), nc_conservation_mag(i,j)

          write(fid,fmt) "qr_conservation(",i,",",nz-j+1,")                                ",& 
               qr_conservation(i,j), qr_conservation_mag(i,j)

          write(fid,fmt) "nr_conservation(",i,",",nz-j+1,")                                ",& 
               nr_conservation(i,j), nr_conservation_mag(i,j)

          write(fid,fmt) "qi_conservation(",i,",",nz-j+1,")                                ",& 
               qi_conservation(i,j), qi_conservation_mag(i,j)

          write(fid,fmt) "ni_conservation(",i,",",nz-j+1,")                                ",& 
               ni_conservation(i,j), ni_conservation_mag(i,j)

          write(fid,fmt) "qs_conservation(",i,",",nz-j+1,")                                ",& 
               qs_conservation(i,j), qs_conservation_mag(i,j)

          write(fid,fmt) "ns_conservation(",i,",",nz-j+1,")                                ",& 
               ns_conservation(i,j), ns_conservation_mag(i,j)

          write(fid,fmt) "rain_evaporation_limiter(",i,",",nz-j+1,")                       ",& 
               qiqs_sublimation_qr_evaporation_limiter(i,j), &
               rain_evaporation_limiter_mag(i,j)

          write(fid,fmt) "snow_sublimation_limiter(",i,",",nz-j+1,")                       ",& 
               qiqs_sublimation_qr_evaporation_limiter(i,j), &
               snow_sublimation_limiter_mag(i,j)

          write(fid,fmt) "ice_sublimation_limiter(",i,",",nz-j+1,")                        ",& 
               qiqs_sublimation_qr_evaporation_limiter(i,j), &
               ice_sublimation_limiter_mag(i,j)

          write(fid,fmt) "ni_tendency_limiter(",i,",",nz-j+1,")                            ",& 
               ni_tendency_limiter(i,j), ni_tendency_limiter_mag(i,j)

          ! write(fid,fmt) "nr_evaporation_limiter(",i,",",nz-j+1,")                         ",& 
          !      nr_evaporation_limiter(i,j), nr_evaporation_limiter_mag(i,j)

       end do
       write(fid,'(A)') "--------------------------------------------------------------------------------"
    end do

    close(fid)
    
  end subroutine write_stats



  subroutine get_stats(array1, array2, ntrue, max, min, avg)
    
    integer, intent(in)          :: array1(nx,nz)
    double precision, intent(in) :: array2(nx,nz)

    integer, intent(out)          :: ntrue
    double precision, intent(out) :: max, min, avg
    
    if (any(array1 == 1)) then
       ntrue = sum(array1)
       max   = maxval(array2, MASK = array1==1)
       min   = minval(array2, MASK = array1==1)
       avg   = sum(array2, MASK = array1==1) / ntrue
    else
       ntrue = 0
       max   = 0.0d0
       min   = 0.0d0
       avg   = 0.0d0
    end if

  end subroutine get_stats



  subroutine write_stats_summary(step)

    integer, intent(in) :: step
    
    character(len=100) :: fname  ! output file name
    integer            :: fid    ! output file id
    
    integer          :: nlim
    double precision :: maxmag, minmag, avgmag

    ! output format
    character(len=100) :: fmt='(A30,4X,I4,4X,E23.15,4X,E23.15,4X,E23.15)'

    ! update number of calls
    ncalls = ncalls + 1
       
    fid = 101
    fname = trim(fileNameOut)//'.limiters.txt'
    open(unit = fid, file = fname, status="old", position="append", action="write")
    
    write(fid,'(A)')    "================================================================================"
    write(fid,'(A,I0)') "Step ",step
    write(fid,'(A)')    "================================================================================"

    call get_stats(qric_limiter, qric_limiter_mag, nlim, maxmag, minmag, avgmag)       
    write(fid,fmt) "qric_limiter", nlim, avgmag, maxmag, minmag
    qric_limiter_mag_max = max(qric_limiter_mag_max, maxmag)

    call get_stats(qric_qsmall, qric_qsmall_mag, nlim, maxmag, minmag, avgmag)       
    write(fid,fmt) "qric_qsmall", nlim, avgmag, maxmag, minmag
    qric_qsmall_mag_max = max(qric_qsmall_mag_max, maxmag)

    call get_stats(nric_qsmall, nric_qsmall_mag, nlim, maxmag, minmag, avgmag)       
    write(fid,fmt) "nric_qsmall", nlim, avgmag, maxmag, minmag
    nric_qsmall_mag_max = max(nric_qsmall_mag_max, maxmag)

    call get_stats(nric_neg, nric_neg_mag, nlim, maxmag, minmag, avgmag)       
    write(fid,fmt) "nric_neg", nlim, avgmag, maxmag, minmag
    nric_neg_mag_max = max(nric_neg_mag_max, maxmag)

    call get_stats(qsic_limiter, qsic_limiter_mag, nlim, maxmag, minmag, avgmag)       
    write(fid,fmt) "qsic_limiter", nlim, avgmag, maxmag, minmag
    qsic_limiter_mag_max = max(qsic_limiter_mag_max, maxmag)

    call get_stats(qsic_qsmall, qsic_qsmall_mag, nlim, maxmag, minmag, avgmag)       
    write(fid,fmt) "qsic_qsmall", nlim, avgmag, maxmag, minmag
    qsic_qsmall_mag_max = max(qsic_qsmall_mag_max, maxmag)

    call get_stats(nsic_qsmall, nsic_qsmall_mag, nlim, maxmag, minmag, avgmag)       
    write(fid,fmt) "nsic_qsmall", nlim, avgmag, maxmag, minmag
    nsic_qsmall_mag_max = max(nsic_qsmall_mag_max, maxmag)

    call get_stats(nsic_neg, nsic_neg_mag, nlim, maxmag, minmag, avgmag)       
    write(fid,fmt) "nsic_neg", nlim, avgmag, maxmag, minmag
    nsic_neg_mag_max = max(nsic_neg_mag_max, maxmag)

    call get_stats(qc_conservation, qc_conservation_mag, nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "qc_conservation", nlim, avgmag, maxmag, minmag
    qc_conservation_mag_max = max(qc_conservation_mag_max, maxmag)

    call get_stats(ice_nucleation_limiter, ice_nucleation_limiter_mag, &
         nlim, maxmag, minmag, avgmag)
    write(fid,fmt) "ice_nucleation_limiter", nlim, avgmag, maxmag, minmag
    ice_nucleation_limiter_mag_max = max(ice_nucleation_limiter_mag_max, maxmag)

    call get_stats(ice_deposition_limiter, ice_deposition_limiter_mag, &
         nlim, maxmag, minmag, avgmag)
    write(fid,fmt) "ice_deposition_limiter", nlim, avgmag, maxmag, minmag
    ice_deposition_limiter_mag_max = max(ice_deposition_limiter_mag_max, maxmag)
    
    call get_stats(nc_conservation, nc_conservation_mag, nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "nc_conservation", nlim, avgmag, maxmag, minmag
    nc_conservation_mag_max = max(nc_conservation_mag_max, maxmag)

    call get_stats(qr_conservation, qr_conservation_mag, nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "qr_conservation", nlim, avgmag, maxmag, minmag
    qr_conservation_mag_max = max(qr_conservation_mag_max, maxmag)

    call get_stats(nr_conservation, nr_conservation_mag, nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "nr_conservation", nlim, avgmag, maxmag, minmag
    nr_conservation_mag_max = max(nr_conservation_mag_max, maxmag)

    call get_stats(qi_conservation, qi_conservation_mag, nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "qi_conservation", nlim, avgmag, maxmag, minmag
    qi_conservation_mag_max = max(qi_conservation_mag_max, maxmag)

    call get_stats(ni_conservation, ni_conservation_mag, nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "ni_conservation", nlim, avgmag, maxmag, minmag
    ni_conservation_mag_max = max(ni_conservation_mag_max, maxmag)

    call get_stats(qs_conservation, qs_conservation_mag, nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "qs_conservation", nlim, avgmag, maxmag, minmag
    qs_conservation_mag_max = max(qs_conservation_mag_max, maxmag)

    call get_stats(ns_conservation, ns_conservation_mag, nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "ns_conservation", nlim, avgmag, maxmag, minmag
    ns_conservation_mag_max = max(ns_conservation_mag_max, maxmag)

    call get_stats(qiqs_sublimation_qr_evaporation_limiter, rain_evaporation_limiter_mag, &
         nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "rain_evaporation_limiter_mag", nlim, avgmag, maxmag, minmag
    rain_evaporation_limiter_mag_max = max(rain_evaporation_limiter_mag_max, maxmag)

    call get_stats(qiqs_sublimation_qr_evaporation_limiter, snow_sublimation_limiter_mag, &
         nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "snow_sublimation_limiter_mag", nlim, avgmag, maxmag, minmag
    snow_sublimation_limiter_mag_max = max(snow_sublimation_limiter_mag_max, maxmag) 

    call get_stats(qiqs_sublimation_qr_evaporation_limiter, ice_sublimation_limiter_mag, &
         nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "ice_sublimation_limiter_mag", nlim, avgmag, maxmag, minmag
    ice_sublimation_limiter_mag_max = max(ice_sublimation_limiter_mag_max, maxmag)

    call get_stats(ni_tendency_limiter, ni_tendency_limiter_mag, &
         nlim, maxmag, minmag, avgmag) 
    write(fid,fmt) "ni_tendency_limiter", nlim, avgmag, maxmag, minmag
    ni_tendency_limiter_mag_max = max(ni_tendency_limiter_mag_max, maxmag)

    close(fid)
    
  end subroutine write_stats_summary


  
  subroutine finalize_stats

    character(100) :: fname ! output file name
    integer        :: fid   ! output file id

    character(len=100) :: fmt='(A30,4X,E23.15)'

    ! add sum over columns as a summary to end of file

    fid = 101
    fname = trim(fileNameOut)//'.limiters.txt'
    open(unit = fid, file = fname, status="old", position="append", action="write")
    
    write(fid,'(A)') "================================================================================"
    write(fid,'(A)') "End of Run"
    write(fid,'(A)') "================================================================================"
    write(fid,'(A)') "Max Values Over Time"
    write(fid,'(A)') "================================================================================"

    write(fid,fmt) "qric_limiter", qric_limiter_mag_max
    write(fid,fmt) "qric_qsmall", qric_qsmall_mag_max
    write(fid,fmt) "nric_qsmall", nric_qsmall_mag_max
    write(fid,fmt) "nric_neg", nric_neg_mag_max
    write(fid,fmt) "qsic_limiter", qsic_limiter_mag_max
    write(fid,fmt) "qsic_qsmall", qsic_qsmall_mag_max
    write(fid,fmt) "nsic_qsmall", nsic_qsmall_mag_max
    write(fid,fmt) "nsic_neg", nsic_neg_mag_max
    write(fid,fmt) "qc_conservation", qc_conservation_mag_max
    write(fid,fmt) "ice_nucleation_limiter", ice_nucleation_limiter_mag_max 
    write(fid,fmt) "ice_deposition_limiter", ice_deposition_limiter_mag_max
    write(fid,fmt) "nc_conservation", nc_conservation_mag_max
    write(fid,fmt) "qr_conservation", qr_conservation_mag_max
    write(fid,fmt) "nr_conservation", nr_conservation_mag_max
    write(fid,fmt) "qi_conservation", qi_conservation_mag_max
    write(fid,fmt) "ni_conservation", ni_conservation_mag_max
    write(fid,fmt) "qs_conservation", qs_conservation_mag_max
    write(fid,fmt) "ns_conservation", ns_conservation_mag_max
    write(fid,fmt) "rain_evaporation_limiter_mag", rain_evaporation_limiter_mag_max
    write(fid,fmt) "snow_sublimation_limiter_mag", snow_sublimation_limiter_mag_max
    write(fid,fmt) "ice_sublimation_limiter_mag", ice_sublimation_limiter_mag_max
    write(fid,fmt) "ni_tendency_limiter", ni_tendency_limiter_mag_max

    write(fid,'(A)') "================================================================================"
    
    close(fid)
    
  end subroutine finalize_stats
  
end module mphys_stats
