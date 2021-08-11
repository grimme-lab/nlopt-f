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

!> Implementation of the callback functionality for the NLopt library.
!>
!> This module provides ISO-C compatible callback functions which are
!> passed to the NLopt library and redirect the computation to a regular
!> Fortran procedure, which in turn can be more generic since it does
!> have to adhere to the constraints placed on ISO-C compatible functions.
!>
!> This is archived by providing small types, for the actual procedure reference
!> and an optional polymorphic data container, which are carried along in
!> NLopt as opaque data pointer and resolved again in the callback wrapper
!> before invoking the actual procedure.
!>
!> While the wrapper types reference the data, they do not own the data,
!> therefore the caller side is responsible to keep the object alive as long
!> it is referenced from the [[nlopt_opt]] instance. One option is to create
!> the wrappers by association:
!>
!>```fortran
!>associate(f => nlopt_func(objective_function, data_container)
!>  call create(opt, algoritm, n)
!>  ! ...
!>  call opt%set_min_objective(f)
!>  ! ...
!>  call opt%optimize(stat)
!>  ! ...
!>  call destroy(opt)
!>end associate
!>```
module nlopt_callback
  use, intrinsic :: iso_c_binding
  implicit none

  abstract interface
    function nlopt_func_interface(x, gradient, func_data) result(f)
      import :: c_int, c_double, c_ptr
      implicit none
      real(c_double), intent(in) :: x(:)
      real(c_double), intent(inout), optional :: gradient(:)
      class(*), intent(in), optional :: func_data
      real(c_double) :: f
    end function

    subroutine nlopt_mfunc_interface(result, x, gradient, func_data)
      import :: c_int, c_double, c_ptr
      implicit none
      real(c_double), intent(inout) :: result(:)
      real(c_double), intent(in) :: x(:)
      real(c_double), intent(inout), optional :: gradient(:)
      class(*), intent(in), optional :: func_data
    end subroutine

    subroutine nlopt_precond_interface(x, v, vpre, func_data)
      import :: c_int, c_double
      implicit none
      real(c_double), intent(in) :: x(:)
      real(c_double), intent(in) :: v(:)
      real(c_double), intent(inout) :: vpre(:)
      class(*), intent(in), optional :: func_data
    end subroutine
  end interface

  type :: nlopt_func
    private
    procedure(nlopt_func_interface), nopass, pointer :: f => null()
    class(*), pointer :: f_data => null()
  end type nlopt_func

  interface nlopt_func
    module procedure :: create_nlopt_func
  end interface nlopt_func

  type :: nlopt_mfunc
    private
    procedure(nlopt_mfunc_interface), nopass, pointer :: f => null()
    class(*), pointer :: f_data => null()
  end type nlopt_mfunc

  interface nlopt_mfunc
    module procedure :: create_nlopt_mfunc
  end interface nlopt_mfunc

  type :: nlopt_precond
    private
    procedure(nlopt_func_interface), nopass, pointer :: f => null()
    procedure(nlopt_precond_interface), nopass, pointer :: pre => null()
    class(*), pointer :: f_data => null()
  end type nlopt_precond

  interface nlopt_precond
    module procedure :: create_nlopt_precond
  end interface nlopt_precond

contains


  function create_nlopt_func(f, f_data) result(self)
    procedure(nlopt_func_interface) :: f
    class(*), intent(in), target, optional :: f_data
    type(nlopt_func) :: self

    self%f => f
    nullify(self%f_data)
    if (present(f_data)) self%f_data => f_data
  end function create_nlopt_func


  function wrap_func(n, x, gradient, func_data) result(f) &
      bind(c, name="nlopt_wrap_func")
    integer(c_int), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(inout), optional :: gradient(n)
    type(c_ptr), value :: func_data
    real(c_double) :: f

    type(nlopt_func), pointer :: func
    real(c_double), allocatable :: stub(:)

    call c_f_pointer(func_data, func)
    f = func%f(x, gradient, func%f_data)
  end function


  function create_nlopt_mfunc(f, f_data) result(self)
    procedure(nlopt_mfunc_interface) :: f
    class(*), intent(in), target, optional :: f_data
    type(nlopt_mfunc) :: self

    self%f => f
    nullify(self%f_data)
    if (present(f_data)) self%f_data => f_data
  end function create_nlopt_mfunc

  subroutine wrap_mfunc(m, result, n, x, gradient, func_data) &
      bind(c, name="nlopt_wrap_mfunc")
    integer(c_int), value :: m
    real(c_double), intent(inout) :: result(m)
    integer(c_int), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(inout), optional :: gradient(n)
    type(c_ptr), value :: func_data

    type(nlopt_mfunc), pointer :: mfunc

    call c_f_pointer(func_data, mfunc)
    call mfunc%f(result, x, gradient, mfunc%f_data)
  end subroutine


  function create_nlopt_precond(f, pre, f_data) result(self)
     procedure(nlopt_func_interface) :: f
     procedure(nlopt_precond_interface) :: pre
     class(*), intent(in), target, optional :: f_data
     type(nlopt_precond) :: self

     self%f => f
     self%pre => pre
    nullify(self%f_data)
     if (present(f_data)) self%f_data => f_data
  end function create_nlopt_precond

  subroutine wrap_precond(n, x, v, vpre, func_data) &
      bind(c, name="nlopt_wrap_precond")
    integer(c_int), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(in) :: v(n)
    real(c_double), intent(inout) :: vpre(n)
    type(c_ptr), value :: func_data

    type(nlopt_precond), pointer :: precond

    call c_f_pointer(func_data, precond)
    call precond%pre(x, v, vpre, precond%f_data)
  end subroutine wrap_precond

end module nlopt_callback
