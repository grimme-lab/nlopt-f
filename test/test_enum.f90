! This file is part of nlopt-f.
! SPDX-Identifier: Apache-2.0 OR MIT
!
! Licensed under either of Apache License, Version 2.0 or MIT license
! at your option; you may not use this file except in compliance with
! the License.
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module test_enum
  use testdrive, only : unittest_type, new_unittest, error_type, check
  use nlopt_enum
  implicit none
  private

  public :: collect_enum

  integer, parameter :: algorithms(*) = [&
    NLOPT_GN_DIRECT, &
    NLOPT_GN_DIRECT_L, &
    NLOPT_GN_DIRECT_L_RAND, &
    NLOPT_GN_DIRECT_NOSCAL, &
    NLOPT_GN_DIRECT_L_NOSCAL, &
    NLOPT_GN_DIRECT_L_RAND_NOSCAL, &
    NLOPT_GN_ORIG_DIRECT, &
    NLOPT_GN_ORIG_DIRECT_L, &
    NLOPT_GD_STOGO, &
    NLOPT_GD_STOGO_RAND, &
    NLOPT_LD_LBFGS_NOCEDAL, &
    NLOPT_LD_LBFGS, &
    NLOPT_LN_PRAXIS, &
    NLOPT_LD_VAR1, &
    NLOPT_LD_VAR2, &
    NLOPT_LD_TNEWTON, &
    NLOPT_LD_TNEWTON_RESTART, &
    NLOPT_LD_TNEWTON_PRECOND, &
    NLOPT_LD_TNEWTON_PRECOND_RESTART, &
    NLOPT_GN_CRS2_LM, &
    NLOPT_GN_MLSL, &
    NLOPT_GD_MLSL, &
    NLOPT_GN_MLSL_LDS, &
    NLOPT_GD_MLSL_LDS, &
    NLOPT_LD_MMA, &
    NLOPT_LN_COBYLA, &
    NLOPT_LN_NEWUOA, &
    NLOPT_LN_NEWUOA_BOUND, &
    NLOPT_LN_NELDERMEAD, &
    NLOPT_LN_SBPLX, &
    NLOPT_LN_AUGLAG, &
    NLOPT_LD_AUGLAG, &
    NLOPT_LN_AUGLAG_EQ, &
    NLOPT_LD_AUGLAG_EQ, &
    NLOPT_LN_BOBYQA, &
    NLOPT_GN_ISRES, &
    NLOPT_AUGLAG, &
    NLOPT_AUGLAG_EQ, &
    NLOPT_G_MLSL, &
    NLOPT_G_MLSL_LDS, &
    NLOPT_LD_SLSQP, &
    NLOPT_LD_CCSAQ, &
    NLOPT_GN_ESCH, &
    NLOPT_GN_AGS]

  integer(nlopt_result), parameter :: results(*) = [&
    NLOPT_FAILURE, &
    NLOPT_INVALID_ARGS, &
    NLOPT_OUT_OF_MEMORY, &
    NLOPT_ROUNDOFF_LIMITED, &
    NLOPT_FORCED_STOP, &
    NLOPT_SUCCESS, &
    NLOPT_STOPVAL_REACHED, &
    NLOPT_FTOL_REACHED, &
    NLOPT_XTOL_REACHED, &
    NLOPT_MAXEVAL_REACHED, &
    NLOPT_MAXTIME_REACHED]

contains

!> Collect all exported unit tests
subroutine collect_enum(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
     new_unittest("result-string", test_result_string), &
     new_unittest("algorithm-name", test_algorithm_name), &
     new_unittest("algorithm-invalid", test_algorithm_invalid), &
     new_unittest("algorithm-string", test_algorithm_string) &
     ]
end subroutine collect_enum

subroutine test_algorithm_name(error)
  !> Error handling
  type(error_type), allocatable, intent(out) :: error

  character(len=:), allocatable :: str
  integer :: i

  do i = 1, size(algorithms)
    str = algorithm_name(algorithms(i))
    call check(error, len(str) > 0)
    if (allocated(error)) exit
  end do
  if (allocated(error)) return

end subroutine test_algorithm_name

subroutine test_algorithm_string(error)
  !> Error handling
  type(error_type), allocatable, intent(out) :: error

  character(len=:), allocatable :: str
  integer :: i

  do i = 1, size(algorithms)
    str = algorithm_to_string(algorithms(i))
    call check(error, algorithm_from_string(str), algorithms(i))
    if (allocated(error)) exit
  end do
  if (allocated(error)) return

end subroutine test_algorithm_string

subroutine test_algorithm_invalid(error)
  !> Error handling
  type(error_type), allocatable, intent(out) :: error

  character(len=:), allocatable :: str
  integer :: i

  str = algorithm_to_string(huge(i))
  call check(error, allocated(str), "Return value for unknown algorithm should be allocated")
  if (allocated(error)) return
  call check(error, len(str) == 0, "String for unknown algorithm should be empty")
  if (allocated(error)) return

  call check(error, algorithm_from_string(""), -1, &
    & "Invalid algorithm identifier should yield -1 as identifier")
  if (allocated(error)) return

end subroutine test_algorithm_invalid

subroutine test_result_string(error)
  !> Error handling
  type(error_type), allocatable, intent(out) :: error

  character(len=:), allocatable :: str
  integer :: i

  do i = 1, size(results)
    str = result_to_string(results(i))
    call check(error, result_from_string(str), results(i))
    if (allocated(error)) exit
  end do
  if (allocated(error)) return

end subroutine test_result_string

end module test_enum
