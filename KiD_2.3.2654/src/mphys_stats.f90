module mphys_stats

  ! This module contains counters for run accumlating run time stats
  !
  ! Note: Array values are set in the microphysics and so are idexed as top down
  ! while kid assumes bottom up ordering in arrays. So need to flip arrays before
  ! outputting with netcdf
  
  use parameters, only : nx, nz
  use switches, only   : fname_stats
  
  implicit none
  private
  save

  public :: init_stats, finalize_stats
  public :: write_stats
  
  ! file id
  integer :: fid
  character(100) :: fname ! output file

  ! limit on max rain mixing ratio
  integer,          public, dimension(nx,nz) :: qric_limiter
  double precision, public, dimension(nx,nz) :: qric_limiter_mag

  ! setting small rain values to zero
  integer,          public, dimension(nx,nz) :: qric_qsmall
  integer,          public, dimension(nx,nz) :: nric_qsmall
  double precision, public, dimension(nx,nz) :: qric_qsmall_mag
  double precision, public, dimension(nx,nz) :: nric_qsmall_mag

  ! setting negative rain number values to zero
  integer,          public, dimension(nx,nz) :: nric_neg
  double precision, public, dimension(nx,nz) :: nric_neg_mag
  
  ! limit on max snow mixing ratio
  integer,          public, dimension(nx,nz) :: qsic_limiter
  double precision, public, dimension(nx,nz) :: qsic_limiter_mag

  ! setting small snow values to zero
  integer,          public, dimension(nx,nz) :: qsic_qsmall
  integer,          public, dimension(nx,nz) :: nsic_qsmall
  double precision, public, dimension(nx,nz) :: qsic_qsmall_mag
  double precision, public, dimension(nx,nz) :: nsic_qsmall_mag

  ! setting negative snow number values to zero
  integer,          public, dimension(nx,nz) :: nsic_neg
  double precision, public, dimension(nx,nz) :: nsic_neg_mag

  ! conservation limiters
  integer,          public, dimension(nx,nz) :: qc_conservation
  double precision, public, dimension(nx,nz) :: qc_conservation_mag  

  integer,          public, dimension(nx,nz) :: ice_nucleation_limiter
  double precision, public, dimension(nx,nz) :: ice_nucleation_limiter_mag

  integer,          public, dimension(nx,nz) :: ice_deposition_limiter
  double precision, public, dimension(nx,nz) :: ice_deposition_limiter_mag

  integer,          public, dimension(nx,nz) :: nc_conservation
  double precision, public, dimension(nx,nz) :: nc_conservation_mag  

  integer,          public, dimension(nx,nz) :: qr_conservation
  double precision, public, dimension(nx,nz) :: qr_conservation_mag  

  integer,          public, dimension(nx,nz) :: nr_conservation
  double precision, public, dimension(nx,nz) :: nr_conservation_mag  

  integer,          public, dimension(nx,nz) :: qi_conservation
  double precision, public, dimension(nx,nz) :: qi_conservation_mag  

  integer,          public, dimension(nx,nz) :: ni_conservation
  double precision, public, dimension(nx,nz) :: ni_conservation_mag  

  integer,          public, dimension(nx,nz) :: qs_conservation
  double precision, public, dimension(nx,nz) :: qs_conservation_mag  

  integer,          public, dimension(nx,nz) :: ns_conservation
  double precision, public, dimension(nx,nz) :: ns_conservation_mag  

  ! these conservation limiters appeared earlier in an older MG2
  integer,          public, dimension(nx,nz) :: qiqs_sublimation_qr_evaporation_limiter
  double precision, public, dimension(nx,nz) :: rain_evaporation_limiter_mag 
  double precision, public, dimension(nx,nz) :: snow_sublimation_limiter_mag 
  double precision, public, dimension(nx,nz) :: ice_sublimation_limiter_mag  
  
  integer,          public, dimension(nx,nz) :: ni_tendency_limiter
  double precision, public, dimension(nx,nz) :: ni_tendency_limiter_mag

  ! integer,          public, dimension(nx,nz) :: nr_evap_limiter
  ! double precision, public, dimension(nx,nz) :: nr_evap_limiter_mag
    
contains
  
  subroutine init_stats()
    
    fid = 101
    fname = trim(fname_stats) // '_limiters.txt'
    open(unit = fid, file = fname)
    write(fid,'(A)') "================================================================================"
    write(fid,'(A)') "Limiter Counts"

    qric_limiter = 0.0d0
    qric_limiter_mag  = 0.0d0
    
    qric_qsmall = 0.0d0
    nric_qsmall = 0.0d0
    qric_qsmall_mag = 0.0d0
    nric_qsmall_mag = 0.0d0
    
    nric_neg = 0.0d0
    nric_neg_mag = 0.0d0
    
    qsic_limiter = 0.0d0
    qsic_limiter_mag = 0.0d0
    
    qsic_qsmall = 0.0d0
    nsic_qsmall = 0.0d0
    qsic_qsmall_mag = 0.0d0
    nsic_qsmall_mag = 0.0d0
    
    nsic_neg = 0.0d0
    nsic_neg_mag = 0.0d0
    
    qc_conservation = 0.0d0
    qc_conservation_mag = 0.0d0
    
    ice_nucleation_limiter = 0.0d0
    ice_nucleation_limiter_mag = 0.0d0
    
    ice_deposition_limiter = 0.0d0
    ice_deposition_limiter_mag = 0.0d0
    
    nc_conservation = 0.0d0
    nc_conservation_mag = 0.0d0
    
    qr_conservation = 0.0d0
    qr_conservation_mag = 0.0d0  
    
    nr_conservation = 0.0d0
    nr_conservation_mag = 0.0d0  
    
    qi_conservation = 0.0d0
    qi_conservation_mag = 0.0d0  
    
    ni_conservation = 0.0d0
    ni_conservation_mag = 0.0d0  
    
    qs_conservation = 0.0d0
    qs_conservation_mag = 0.0d0  
    
    ns_conservation = 0.0d0
    ns_conservation_mag = 0.0d0  
    
    qiqs_sublimation_qr_evaporation_limiter = 0.0d0
    rain_evaporation_limiter_mag = 0.0d0 
    snow_sublimation_limiter_mag = 0.0d0 
    ice_sublimation_limiter_mag = 0.0d0  
    
    ni_tendency_limiter = 0.0d0
    ni_tendency_limiter_mag = 0.0d0

  end subroutine init_stats

  subroutine write_stats(step)

    integer, intent(in) :: step
    
    integer :: i,j
    character(len=25) :: fmt

    fmt = '(A,I0,A,I0,A,I0,E23.15)'

    write(fid,'(A)')      "================================================================================"
    write(fid,'(A,I0,A)') "Step",step," (nx,nz)" 
    write(fid,'(A)')      "================================================================================"
    do j=1,nz
       do i=1,nx

          write(fid,fmt) "qric_limiter(",i,",",j,")                                   ",&
               qric_limiter(i,j), qric_limiter_mag(i,j)

          write(fid,fmt) "qric_qsmall(",i,",",j,")                                    ",& 
               qric_qsmall(i,j), qric_qsmall_mag(i,j)

          write(fid,fmt) "nric_qsmall(",i,",",j,")                                    ",& 
               nric_qsmall(i,j), nric_qsmall_mag(i,j)

          write(fid,fmt) "nric_neg(",i,",",j,")                                       ",& 
               nric_neg(i,j), nric_neg_mag(i,j)

          write(fid,fmt) "qsic_limiter(",i,",",j,")                                   ",& 
               qsic_limiter(i,j), qsic_limiter_mag(i,j)

          write(fid,fmt) "qsic_qsmall(",i,",",j,")                                    ",& 
               qsic_qsmall(i,j), qsic_qsmall_mag(i,j)

          write(fid,fmt) "nsic_qsmall(",i,",",j,")                                    ",& 
               nsic_qsmall(i,j), nsic_qsmall_mag(i,j)

          write(fid,fmt) "nsic_neg(",i,",",j,")                                       ",& 
               nsic_neg(i,j), nsic_neg_mag(i,j)

          write(fid,fmt) "qc_conservation(",i,",",j,")                                ",& 
               qc_conservation(i,j), qc_conservation_mag(i,j)

          write(fid,fmt) "ice_nucleation_limiter(",i,",",j,")                         ",& 
               ice_nucleation_limiter(i,j), ice_nucleation_limiter_mag(i,j)

          write(fid,fmt) "ice_deposition_limiter(",i,",",j,")                         ",& 
               ice_deposition_limiter(i,j), ice_deposition_limiter_mag(i,j)

          write(fid,fmt) "nc_conservation(",i,",",j,")                                ",& 
               nc_conservation(i,j), nc_conservation_mag(i,j)

          write(fid,fmt) "qr_conservation(",i,",",j,")                                ",& 
               qr_conservation(i,j), qr_conservation_mag(i,j)

          write(fid,fmt) "nr_conservation(",i,",",j,")                                ",& 
               nr_conservation(i,j), nr_conservation_mag(i,j)

          write(fid,fmt) "qi_conservation(",i,",",j,")                                ",& 
               qi_conservation(i,j), qi_conservation_mag(i,j)

          write(fid,fmt) "ni_conservation(",i,",",j,")                                ",& 
               ni_conservation(i,j), ni_conservation_mag(i,j)

          write(fid,fmt) "qs_conservation(",i,",",j,")                                ",& 
               qs_conservation(i,j), qs_conservation_mag(i,j)

          write(fid,fmt) "ns_conservation(",i,",",j,")                                ",& 
               ns_conservation(i,j), ns_conservation_mag(i,j)

          write(fid,fmt) "rain_evaporation_limiter(",i,",",j,")                       ",& 
               qiqs_sublimation_qr_evaporation_limiter(i,j), &
               rain_evaporation_limiter_mag(i,j)

          write(fid,fmt) "snow_sublimation_limiter(",i,",",j,")                       ",& 
               qiqs_sublimation_qr_evaporation_limiter(i,j), &
               snow_sublimation_limiter_mag(i,j)

          write(fid,fmt) "ice_sublimation_limiter(",i,",",j,")                        ",& 
               qiqs_sublimation_qr_evaporation_limiter(i,j), &
               ice_sublimation_limiter_mag(i,j)

          write(fid,fmt) "ni_tendency_limiter(",i,",",j,")                            ",& 
               ni_tendency_limiter(i,j), ni_tendency_limiter_mag(i,j)

          ! write(fid,fmt) "nr_evaporation_limiter(",i,",",j,")                         ",& 
          !      nr_evaporation_limiter(i,j), nr_evaporation_limiter_mag(i,j)

       end do
       write(fid,'(A)') "--------------------------------------------------------------------------------"
    end do
    
  end subroutine write_stats
  
  subroutine finalize_stats

    ! add sum over columns as a summary to end of file

    write(fid,'(A)') "================================================================================"
    write(fid,'(A)') "EOF"
    write(fid,'(A)') "================================================================================"
    close(fid)
    
  end subroutine finalize_stats
  
end module mphys_stats
