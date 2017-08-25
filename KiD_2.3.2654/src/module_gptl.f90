module gptl_kid
  implicit none
  public
  save
  integer :: gptl_ret
  integer :: gptl_sed_ice = 0
  integer :: gptl_calccfl_total_ice = 0
  integer :: gptl_lambda_ice = 0
  integer :: gptl_cfl_ice = 0
  integer :: gptl_flux_ice = 0
  integer :: gptl_advance_sol_ice = 0

  integer :: gptl_sed_cloud = 0
  integer :: gptl_calccfl_total_cloud = 0
  integer :: gptl_lambda_cloud = 0
  integer :: gptl_cfl_cloud = 0
  integer :: gptl_flux_cloud = 0
  integer :: gptl_advance_sol_cloud = 0

  integer :: gptl_sed_rain = 0
  integer :: gptl_calccfl_total_rain = 0
  integer :: gptl_lambda_rain = 0
  integer :: gptl_cfl_rain = 0
  integer :: gptl_flux_rain = 0
  integer :: gptl_advance_sol_rain = 0

  integer :: gptl_sed_snow = 0
  integer :: gptl_calccfl_total_snow = 0
  integer :: gptl_lambda_snow = 0
  integer :: gptl_cfl_snow = 0
  integer :: gptl_flux_snow = 0
  integer :: gptl_advance_sol_snow = 0


end module gptl_kid
