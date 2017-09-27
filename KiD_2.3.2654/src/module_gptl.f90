module gptl_kid
  implicit none
  public
  save
  integer :: gptl_ret
  integer :: gptl_loop = 10000

  integer :: gptl_lambda_ice = 0
  integer :: gptl_fall_ice = 0
  integer :: gptl_advance_sol_ice = 0

  integer :: gptl_lambda_cloud = 0
  integer :: gptl_fall_cloud = 0
  integer :: gptl_advance_sol_cloud = 0

  integer :: gptl_lambda_rain = 0
  integer :: gptl_fall_rain = 0
  integer :: gptl_advance_sol_rain = 0

  integer :: gptl_lambda_snow = 0
  integer :: gptl_fall_snow = 0
  integer :: gptl_advance_sol_snow = 0


end module gptl_kid
