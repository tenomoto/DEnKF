program test_mg
  use mg_module, only: mg_prolong_test, mg_vcycle_test
  implicit none

  call mg_prolong_test
  call mg_vcycle_test

end program test_mg
