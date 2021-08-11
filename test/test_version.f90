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

module test_version
  use testdrive, only : unittest_type, new_unittest, error_type, check
  use nlopt_interface, only : nlopt_version
  implicit none
  private

  public :: collect_version

contains

!> Collect all exported unit tests
subroutine collect_version(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
     new_unittest("nlopt-version", test_nlopt_version) &
     ]
end subroutine collect_version

subroutine test_nlopt_version(error)
  !> Error handling
  type(error_type), allocatable, intent(out) :: error

  integer :: major, minor, bugfix

  call nlopt_version(major, minor, bugfix)
  call check(error, major, 2)
end subroutine test_nlopt_version

end module test_version
