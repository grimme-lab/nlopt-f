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

module test_opt
  use nlopt_wrap
  use nlopt_enum
  use testdrive, only : unittest_type, new_unittest, error_type, check
  implicit none
  private

  public :: collect_opt

  integer, parameter :: wp = kind(0.0d0)
  type :: constraint_data
    real(wp) :: d(2)
  end type
contains

!> Collect all exported unit tests
subroutine collect_opt(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
     new_unittest("example", test_example) &
     ]
end subroutine collect_opt

subroutine test_example(error)
  !> Error handling
  type(error_type), allocatable, intent(out) :: error

  type(nlopt_opt) :: opt

  real(wp) :: lb(2), x(2), minf
  integer :: stat
  type(constraint_data), target :: d1, d2
  real(wp), parameter :: xtol = 1.0e-4_wp

  call create(opt, algorithm_from_string("LD_MMA"), 2)

  call opt%get_lower_bounds(lb, stat)
  call check(error, stat, NLOPT_SUCCESS)
  if (allocated(error)) return

  lb(2) = 0.0_wp
  call opt%set_lower_bounds(lb, stat)
  call check(error, stat, NLOPT_SUCCESS)
  if (allocated(error)) return

  d1%d = [+2.0_wp, +0.0_wp]
  d2%d = [-1.0_wp, +1.0_wp]
  associate(&
      & f => nlopt_func(myfunc), &
      & fc1 => nlopt_func(myconstraint, d1), &
      & fc2 => nlopt_func(myconstraint, d2))
    call opt%set_min_objective(f, stat)
    call check(error, stat, NLOPT_SUCCESS)
    if (allocated(error)) return

    call opt%add_inequality_constraint(fc1, 1.0e-8_wp, stat)
    call check(error, stat, NLOPT_SUCCESS)
    if (allocated(error)) return

    call opt%add_inequality_constraint(fc2, 1.0e-8_wp, stat)
    call check(error, stat, NLOPT_SUCCESS)
    if (allocated(error)) return

    call opt%set_xtol_rel(xtol, stat)
    call check(error, stat, NLOPT_SUCCESS)
    if (allocated(error)) return

    x = [1.234_wp, 5.678_wp]
    call opt%optimize(x, minf, stat)
    call check(error, stat, NLOPT_XTOL_REACHED)
    if (allocated(error)) return
  end associate

  call check(error, x(1), 0.33333333_wp, thr=xtol, rel=.true.)
  call check(error, x(2), 0.29629629_wp, thr=xtol, rel=.true.)
  call check(error, minf, 0.54433105_wp, thr=xtol, rel=.true.)

  call destroy(opt)
end subroutine test_example

function myfunc(x, gradient, func_data) result(f)
  real(wp), intent(in) :: x(:)
  real(wp), intent(inout), optional :: gradient(:)
  class(*), intent(in), optional :: func_data
  real(wp) :: f

  if (present(gradient)) then
    gradient(1) = 0.0
    gradient(2) = 0.5 / dsqrt(x(2))
  endif
  f = dsqrt(x(2))
end function myfunc

function myconstraint(x, gradient, func_data) result(f)
  real(wp), intent(in) :: x(:)
  real(wp), intent(inout), optional :: gradient(:)
  class(*), intent(in), optional :: func_data
  real(wp) :: f

  select type(func_data)
  type is(constraint_data)
    associate(a => func_data%d(1), b => func_data%d(2))
      if (present(gradient)) then
        gradient(1) = 3. * a * (a*x(1) + b)**2
        gradient(2) = -1.0
      endif
      f = (a*x(1) + b)**3 - x(2)
    end associate
  end select
end function myconstraint

end module test_opt
