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

module nlopt_wrap
  use, intrinsic :: iso_c_binding, only : c_null_ptr, c_int, c_char, c_double, c_null_char, &
     & c_ptr, c_loc, c_funloc, c_associated, c_f_pointer
  use nlopt_enum, only : &
    NLOPT_GN_DIRECT, NLOPT_GN_DIRECT_L, NLOPT_GN_DIRECT_L_RAND, NLOPT_GN_DIRECT_NOSCAL, &
    NLOPT_GN_DIRECT_L_NOSCAL, NLOPT_GN_DIRECT_L_RAND_NOSCAL, NLOPT_GN_ORIG_DIRECT, &
    NLOPT_GN_ORIG_DIRECT_L, NLOPT_GD_STOGO, NLOPT_GD_STOGO_RAND, NLOPT_LD_LBFGS_NOCEDAL, &
    NLOPT_LD_LBFGS, NLOPT_LN_PRAXIS, NLOPT_LD_VAR1, NLOPT_LD_VAR2, NLOPT_LD_TNEWTON, &
    NLOPT_LD_TNEWTON_RESTART, NLOPT_LD_TNEWTON_PRECOND, NLOPT_LD_TNEWTON_PRECOND_RESTART, &
    NLOPT_GN_CRS2_LM, NLOPT_GN_MLSL, NLOPT_GD_MLSL, NLOPT_GN_MLSL_LDS, NLOPT_GD_MLSL_LDS, &
    NLOPT_LD_MMA, NLOPT_LN_COBYLA, NLOPT_LN_NEWUOA, NLOPT_LN_NEWUOA_BOUND, &
    NLOPT_LN_NELDERMEAD, NLOPT_LN_SBPLX, NLOPT_LN_AUGLAG, NLOPT_LD_AUGLAG, &
    NLOPT_LN_AUGLAG_EQ, NLOPT_LD_AUGLAG_EQ, NLOPT_LN_BOBYQA, NLOPT_GN_ISRES, &
    NLOPT_AUGLAG, NLOPT_AUGLAG_EQ, NLOPT_G_MLSL, NLOPT_G_MLSL_LDS, NLOPT_LD_SLSQP, &
    NLOPT_LD_CCSAQ, NLOPT_GN_ESCH, NLOPT_GN_AGS, NLOPT_NUM_ALGORITHMS, nlopt_algorithm, &
    NLOPT_FAILURE, NLOPT_INVALID_ARGS, NLOPT_OUT_OF_MEMORY, NLOPT_ROUNDOFF_LIMITED, &
    NLOPT_FORCED_STOP, NLOPT_SUCCESS, NLOPT_STOPVAL_REACHED, NLOPT_FTOL_REACHED, &
    NLOPT_XTOL_REACHED, NLOPT_MAXEVAL_REACHED, NLOPT_MAXTIME_REACHED, NLOPT_NUM_RESULTS, &
    nlopt_result
  use nlopt_interface, only : &
    nlopt_srand, nlopt_srand_time, &
    nlopt_version, nlopt_create, nlopt_destroy, nlopt_copy, nlopt_optimize, &
    nlopt_set_min_objective, nlopt_set_max_objective, nlopt_set_precond_min_objective, &
    nlopt_set_precond_max_objective, nlopt_get_algorithm, nlopt_get_dimension, &
    nlopt_get_errmsg, nlopt_set_param, nlopt_get_param, nlopt_has_param, &
    nlopt_num_params, nlopt_nth_param, nlopt_set_lower_bounds, nlopt_set_lower_bounds1, &
    nlopt_set_lower_bound, nlopt_get_lower_bounds, nlopt_set_upper_bounds, &
    nlopt_set_upper_bounds1, nlopt_set_upper_bound, nlopt_get_upper_bounds, &
    nlopt_remove_inequality_constraints, nlopt_add_inequality_constraint, &
    nlopt_add_precond_inequality_constraint, nlopt_add_inequality_mconstraint, &
    nlopt_remove_equality_constraints, nlopt_add_equality_constraint, &
    nlopt_add_precond_equality_constraint, nlopt_add_equality_mconstraint, &
    nlopt_set_stopval, nlopt_get_stopval, nlopt_set_ftol_rel, nlopt_get_ftol_rel, &
    nlopt_set_ftol_abs, nlopt_get_ftol_abs, nlopt_set_xtol_rel, nlopt_get_xtol_rel, &
    nlopt_set_xtol_abs1, nlopt_set_xtol_abs, nlopt_get_xtol_abs, nlopt_set_x_weights1, &
    nlopt_set_x_weights, nlopt_get_x_weights, nlopt_set_maxeval, nlopt_get_maxeval, &
    nlopt_get_numevals, nlopt_set_maxtime, nlopt_get_maxtime, nlopt_force_stop, &
    nlopt_set_force_stop, nlopt_get_force_stop, nlopt_set_local_optimizer, &
    nlopt_set_population, nlopt_get_population, nlopt_set_vector_storage, &
    nlopt_get_vector_storage, nlopt_set_default_initial_step, nlopt_set_initial_step, &
    nlopt_set_initial_step1, nlopt_get_initial_step
  use nlopt_callback, only : nlopt_func, nlopt_mfunc, nlopt_precond, &
    wrap_func, wrap_mfunc, wrap_precond
  implicit none
  private

  public :: nlopt_opt, create, destroy
  public :: nlopt_result_enum, nlopt_algorithm_enum
  public :: nlopt_func, nlopt_mfunc, nlopt_precond

  type :: enum_result
    integer :: &
      FAILURE          = NLOPT_FAILURE, &
      INVALID_ARGS     = NLOPT_INVALID_ARGS, &
      OUT_OF_MEMORY    = NLOPT_OUT_OF_MEMORY, &
      ROUNDOFF_LIMITED = NLOPT_ROUNDOFF_LIMITED, &
      FORCED_STOP      = NLOPT_FORCED_STOP, &
      SUCCESS          = NLOPT_SUCCESS, &
      STOPVAL_REACHED  = NLOPT_STOPVAL_REACHED, &
      FTOL_REACHED     = NLOPT_FTOL_REACHED, &
      XTOL_REACHED     = NLOPT_XTOL_REACHED, &
      MAXEVAL_REACHED  = NLOPT_MAXEVAL_REACHED, &
      MAXTIME_REACHED  = NLOPT_MAXTIME_REACHED
  end type
  type(enum_result) :: nlopt_result_enum = enum_result()

  type :: enum_algorithm
    integer :: &
      GN_DIRECT                  = NLOPT_GN_DIRECT, &
      GN_DIRECT_L                = NLOPT_GN_DIRECT_L, &
      GN_DIRECT_L_RAND           = NLOPT_GN_DIRECT_L_RAND, &
      GN_DIRECT_NOSCAL           = NLOPT_GN_DIRECT_NOSCAL, &
      GN_DIRECT_L_NOSCAL         = NLOPT_GN_DIRECT_L_NOSCAL, &
      GN_DIRECT_L_RAND_NOSCAL    = NLOPT_GN_DIRECT_L_RAND_NOSCAL, &
      GN_ORIG_DIRECT             = NLOPT_GN_ORIG_DIRECT, &
      GN_ORIG_DIRECT_L           = NLOPT_GN_ORIG_DIRECT_L, &
      GD_STOGO                   = NLOPT_GD_STOGO, &
      GD_STOGO_RAND              = NLOPT_GD_STOGO_RAND, &
      LD_LBFGS_NOCEDAL           = NLOPT_LD_LBFGS_NOCEDAL, &
      LD_LBFGS                   = NLOPT_LD_LBFGS, &
      LN_PRAXIS                  = NLOPT_LN_PRAXIS, &
      LD_VAR1                    = NLOPT_LD_VAR1, &
      LD_VAR2                    = NLOPT_LD_VAR2, &
      LD_TNEWTON                 = NLOPT_LD_TNEWTON, &
      LD_TNEWTON_RESTART         = NLOPT_LD_TNEWTON_RESTART, &
      LD_TNEWTON_PRECOND         = NLOPT_LD_TNEWTON_PRECOND, &
      LD_TNEWTON_PRECOND_RESTART = NLOPT_LD_TNEWTON_PRECOND_RESTART, &
      GN_CRS2_LM                 = NLOPT_GN_CRS2_LM, &
      GN_MLSL                    = NLOPT_GN_MLSL, &
      GD_MLSL                    = NLOPT_GD_MLSL, &
      GN_MLSL_LDS                = NLOPT_GN_MLSL_LDS, &
      GD_MLSL_LDS                = NLOPT_GD_MLSL_LDS, &
      LD_MMA                     = NLOPT_LD_MMA, &
      LN_COBYLA                  = NLOPT_LN_COBYLA, &
      LN_NEWUOA                  = NLOPT_LN_NEWUOA, &
      LN_NEWUOA_BOUND            = NLOPT_LN_NEWUOA_BOUND, &
      LN_NELDERMEAD              = NLOPT_LN_NELDERMEAD, &
      LN_SBPLX                   = NLOPT_LN_SBPLX, &
      LN_AUGLAG                  = NLOPT_LN_AUGLAG, &
      LD_AUGLAG                  = NLOPT_LD_AUGLAG, &
      LN_AUGLAG_EQ               = NLOPT_LN_AUGLAG_EQ, &
      LD_AUGLAG_EQ               = NLOPT_LD_AUGLAG_EQ, &
      LN_BOBYQA                  = NLOPT_LN_BOBYQA, &
      GN_ISRES                   = NLOPT_GN_ISRES, &
      AUGLAG                     = NLOPT_AUGLAG, &
      AUGLAG_EQ                  = NLOPT_AUGLAG_EQ, &
      G_MLSL                     = NLOPT_G_MLSL, &
      G_MLSL_LDS                 = NLOPT_G_MLSL_LDS, &
      LD_SLSQP                   = NLOPT_LD_SLSQP, &
      LD_CCSAQ                   = NLOPT_LD_CCSAQ, &
      GN_ESCH                    = NLOPT_GN_ESCH, &
      GN_AGS                     = NLOPT_GN_AGS
  end type
  type(enum_algorithm) :: nlopt_algorithm_enum = enum_algorithm()

  integer, parameter :: wp = kind(0.0d0)
  integer, parameter :: ik = kind(0)

  interface create
    module procedure :: create
  end interface

  interface destroy
    module procedure :: destroy
  end interface

  type :: nlopt_opt
    private
    type(c_ptr) :: handle = c_null_ptr
  contains
    procedure, private :: copy
    generic :: assignment(=) => copy
    final :: destroy

    procedure :: optimize
    procedure :: set_min_objective
    procedure :: set_max_objective
    procedure :: set_precond_min_objective
    procedure :: set_precond_max_objective
    procedure :: get_algorithm
    procedure :: get_dimension
    procedure :: get_errmsg
    procedure :: set_param
    procedure :: get_param
    procedure :: has_param
    procedure :: num_params
    procedure :: nth_param
    procedure :: set_lower_bounds
    procedure :: set_lower_bounds1
    procedure :: set_lower_bound
    procedure :: get_lower_bounds
    procedure :: set_upper_bounds
    procedure :: set_upper_bounds1
    procedure :: set_upper_bound
    procedure :: get_upper_bounds
    procedure :: remove_inequality_constraints
    procedure :: add_inequality_constraint
    procedure :: add_precond_inequality_constraint
    procedure :: add_inequality_mconstraint
    procedure :: remove_equality_constraints
    procedure :: add_equality_constraint
    procedure :: add_precond_equality_constraint
    procedure :: add_equality_mconstraint
    procedure :: set_stopval
    procedure :: get_stopval
    procedure :: set_ftol_rel
    procedure :: get_ftol_rel
    procedure :: set_ftol_abs
    procedure :: get_ftol_abs
    procedure :: set_xtol_rel
    procedure :: get_xtol_rel
    procedure :: set_xtol_abs1
    procedure :: set_xtol_abs
    procedure :: get_xtol_abs
    procedure :: set_x_weights1
    procedure :: set_x_weights
    procedure :: get_x_weights
    procedure :: set_maxeval
    procedure :: get_maxeval
    procedure :: get_numevals
    procedure :: set_maxtime
    procedure :: get_maxtime
    procedure :: force_stop
    procedure :: set_force_stop
    procedure :: get_force_stop
    procedure :: set_local_optimizer
    procedure :: set_population
    procedure :: get_population
    procedure :: set_vector_storage
    procedure :: get_vector_storage
    procedure :: set_default_initial_step
    procedure :: set_initial_step
    procedure :: set_initial_step1
    procedure :: get_initial_step
  end type nlopt_opt

  interface nlopt_opt
    module procedure :: create_opt
  end interface nlopt_opt

  interface
    function strlen(str) result(len) bind(c)
      import :: c_ptr, c_int
      type(c_ptr), value :: str
      integer(c_int) :: len
    end function strlen
  end interface

contains

  function create_opt(algorithm, n) result(self)
    integer(nlopt_algorithm), intent(in) :: algorithm
    integer(ik), intent(in) :: n
    type(nlopt_opt) :: self

    call create(self, algorithm, n)
  end function create_opt

  subroutine create(self, algorithm, n)
    type(nlopt_opt), intent(out) :: self
    integer(nlopt_algorithm), intent(in) :: algorithm
    integer(ik), intent(in) :: n
    self%handle = nlopt_create(algorithm, int(n, c_int))
  end subroutine create

  subroutine destroy(self)
    type(nlopt_opt), intent(inout) :: self
    if (c_associated(self%handle)) call nlopt_destroy(self%handle)
    self%handle = c_null_ptr
  end subroutine destroy

  subroutine copy(lhs, rhs)
    class(nlopt_opt), intent(inout) :: lhs
    type(nlopt_opt), intent(in) :: rhs

    call destroy(lhs)
    if (c_associated(rhs%handle)) lhs%handle = nlopt_copy(rhs%handle)
  end subroutine copy


  subroutine handle_result(stat, istat)
    integer(ik), intent(out), optional :: stat
    integer(nlopt_result), intent(in) :: istat

    if (present(stat)) then
      stat = istat
    else
      if (istat < NLOPT_SUCCESS) then
        error stop
      end if
    end if
  end subroutine handle_result


  subroutine optimize(self, x, opt_f, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(inout) :: x(:)
    real(wp), intent(inout) :: opt_f
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_optimize(self%handle, x, opt_f)
    call handle_result(stat, istat)
  end subroutine optimize


  subroutine set_min_objective(self, f, stat)
    class(nlopt_opt), intent(inout) :: self
    type(nlopt_func), intent(in), target :: f
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_min_objective(self%handle, c_funloc(wrap_func), c_loc(f))
    call handle_result(stat, istat)
  end subroutine set_min_objective

  subroutine set_max_objective(self, f, stat)
    class(nlopt_opt), intent(inout) :: self
    type(nlopt_func), intent(in), target :: f
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_max_objective(self%handle, c_funloc(wrap_func), c_loc(f))
    call handle_result(stat, istat)
  end subroutine set_max_objective

  subroutine set_precond_min_objective(self, f, stat)
    class(nlopt_opt), intent(inout) :: self
    type(nlopt_precond), intent(in), target :: f
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_precond_min_objective(self%handle, c_funloc(wrap_func), &
      & c_funloc(wrap_precond), c_loc(f))
    call handle_result(stat, istat)
  end subroutine set_precond_min_objective

  subroutine set_precond_max_objective(self, f, stat)
    class(nlopt_opt), intent(inout) :: self
    type(nlopt_precond), intent(in), target :: f
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_precond_max_objective(self%handle, c_funloc(wrap_func), &
      & c_funloc(wrap_precond), c_loc(f))
    call handle_result(stat, istat)
  end subroutine set_precond_max_objective


  function get_algorithm(self) result(algorithm)
    class(nlopt_opt), intent(inout) :: self
    integer(nlopt_algorithm) :: algorithm

    algorithm = nlopt_get_algorithm(self%handle)
  end function get_algorithm

  function get_dimension(self) result(n)
    class(nlopt_opt), intent(inout) :: self
    integer :: n

    n = nlopt_get_dimension(self%handle)
  end function get_dimension

  function get_errmsg(self) result(errmsg)
    class(nlopt_opt), intent(inout) :: self
    character(len=:, kind=c_char), pointer :: errmsg

    call c_f_pointer(nlopt_get_errmsg(self%handle), errmsg)
    if (associated(errmsg)) errmsg => errmsg(:strlen(c_loc(errmsg)))
  end function get_errmsg


  ! /* generic algorithm parameters: */
  subroutine set_param(self, name, val, stat)
    class(nlopt_opt), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(wp), intent(in) :: val
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_param(self%handle, as_c_char(name), real(val, c_double))
    call handle_result(stat, istat)
  end subroutine set_param

  function get_param(self, name, defaultval) result(val)
    class(nlopt_opt), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(wp), intent(in) :: defaultval
    real(wp) :: val

    val = nlopt_get_param(self%handle, as_c_char(name), real(defaultval, c_double))
  end function get_param

  function has_param(self, name) result(istat)
    class(nlopt_opt), intent(inout) :: self
    character(len=*), intent(in) :: name
    integer(ik) :: istat

    istat = nlopt_has_param(self%handle, as_c_char(name))
  end function has_param

  function num_params(self) result(n)
    class(nlopt_opt), intent(inout) :: self
    integer(ik) :: n

    n = nlopt_num_params(self%handle)
  end function num_params

  function nth_param(self, n) result(name)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(in) :: n
    character(len=:, kind=c_char), pointer :: name

    call c_f_pointer(nlopt_nth_param(self%handle, int(n, c_int)), name)
    if (associated(name)) name => name(:strlen(c_loc(name)))
  end function nth_param


  ! /* constraints: */
  subroutine set_lower_bounds(self, lb, stat)
    class(nlopt_opt), intent(inout) :: self
    real(c_double), intent(in) :: lb(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_lower_bounds(self%handle, lb)
    call handle_result(stat, istat)
  end subroutine set_lower_bounds

  subroutine set_lower_bounds1(self, lb, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(in) :: lb
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_lower_bounds1(self%handle, real(lb, c_double))
    call handle_result(stat, istat)
  end subroutine set_lower_bounds1

  subroutine set_lower_bound(self, i, lb, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(in) :: i
    real(wp), intent(in) :: lb
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat
    istat = nlopt_set_lower_bound(self%handle, int(i, c_int), real(lb, c_double))
    call handle_result(stat, istat)
  end subroutine set_lower_bound

  subroutine get_lower_bounds(self, lb, stat)
    class(nlopt_opt), intent(inout) :: self
    real(c_double), intent(inout) :: lb(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_get_lower_bounds(self%handle, lb)
    call handle_result(stat, istat)
  end subroutine get_lower_bounds


  subroutine set_upper_bounds(self, ub, stat)
    class(nlopt_opt), intent(inout) :: self
    real(c_double), intent(in) :: ub(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_upper_bounds(self%handle, ub)
    call handle_result(stat, istat)
  end subroutine set_upper_bounds

  subroutine set_upper_bounds1(self, ub, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(in) :: ub
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_upper_bounds1(self%handle, real(ub, c_double))
    call handle_result(stat, istat)
  end subroutine set_upper_bounds1

  subroutine set_upper_bound(self, i, ub, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(in) :: i
    real(wp), intent(in) :: ub
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_upper_bound(self%handle, int(i, c_int), real(ub, c_double))
    call handle_result(stat, istat)
  end subroutine set_upper_bound

  subroutine get_upper_bounds(self, ub, stat)
    class(nlopt_opt), intent(inout) :: self
    real(c_double), intent(inout) :: ub(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat
    istat = nlopt_get_upper_bounds(self%handle, ub)
    call handle_result(stat, istat)
  end subroutine get_upper_bounds


  subroutine remove_inequality_constraints(self, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_remove_inequality_constraints(self%handle)
    call handle_result(stat, istat)
  end subroutine remove_inequality_constraints

  subroutine add_inequality_constraint(self, fc, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    type(nlopt_func), intent(in), target :: fc
    real(wp), intent(in) :: tol
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_add_inequality_constraint(self%handle, c_funloc(wrap_func), c_loc(fc), &
      & real(tol, c_double))
    call handle_result(stat, istat)
  end subroutine add_inequality_constraint

  subroutine add_precond_inequality_constraint(self, fc, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    type(nlopt_precond), intent(in), target :: fc
    real(wp), intent(in) :: tol
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_add_precond_inequality_constraint(self%handle, c_funloc(wrap_func), &
      & c_funloc(wrap_precond), c_loc(fc), real(tol, c_double))
    call handle_result(stat, istat)
  end subroutine add_precond_inequality_constraint

  subroutine add_inequality_mconstraint(self, m, fc, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(in) :: m
    type(nlopt_mfunc), intent(in), target :: fc
    real(c_double), intent(in) :: tol(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_add_inequality_mconstraint(self%handle, int(m, c_int), &
      & c_funloc(wrap_mfunc), c_loc(fc), tol)
    call handle_result(stat, istat)
  end subroutine add_inequality_mconstraint


  subroutine remove_equality_constraints(self, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_remove_equality_constraints(self%handle)
    call handle_result(stat, istat)
  end subroutine remove_equality_constraints

  subroutine add_equality_constraint(self, h, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    type(nlopt_func), intent(in), target :: h
    real(wp), intent(in) :: tol
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_add_equality_constraint(self%handle, c_funloc(wrap_func), c_loc(h), &
      & real(tol, c_double))
    call handle_result(stat, istat)
  end subroutine add_equality_constraint

  subroutine add_precond_equality_constraint(self, h, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    type(nlopt_precond), intent(in), target :: h
    real(wp), intent(in) :: tol
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_add_precond_equality_constraint(self%handle, c_funloc(wrap_func), &
      & c_funloc(wrap_precond), c_loc(h), real(tol, c_double))
    call handle_result(stat, istat)
  end subroutine add_precond_equality_constraint

  subroutine add_equality_mconstraint(self, m, h, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(in) :: m
    type(nlopt_mfunc), intent(in), target :: h
    real(c_double), intent(in) :: tol(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_add_equality_mconstraint(self%handle, int(m, c_int), c_funloc(wrap_mfunc), &
      & c_loc(h), tol)
    call handle_result(stat, istat)
  end subroutine add_equality_mconstraint


  ! /* stopping criteria: */
  subroutine set_stopval(self, stopval, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(in) :: stopval
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_stopval(self%handle, real(stopval, c_double))
    call handle_result(stat, istat)
  end subroutine set_stopval

  function get_stopval(self) result(stopval)
    class(nlopt_opt), intent(inout) :: self
    real(wp) :: stopval

    stopval = nlopt_get_stopval(self%handle)
  end function get_stopval


  subroutine set_ftol_rel(self, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(in) :: tol
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_ftol_rel(self%handle, real(tol, c_double))
    call handle_result(stat, istat)
  end subroutine set_ftol_rel

  function get_ftol_rel(self) result(tol)
    class(nlopt_opt), intent(inout) :: self
    real(wp) :: tol

    tol = nlopt_get_ftol_rel(self%handle)
  end function get_ftol_rel

  subroutine set_ftol_abs(self, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(in) :: tol
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_ftol_abs(self%handle, real(tol, c_double))
    call handle_result(stat, istat)
  end subroutine set_ftol_abs

  function get_ftol_abs(self) result(tol)
    class(nlopt_opt), intent(inout) :: self
    real(wp) :: tol

    tol = nlopt_get_ftol_abs(self%handle)
  end function get_ftol_abs


  subroutine set_xtol_rel(self, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(in) :: tol
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_xtol_rel(self%handle, real(tol, c_double))
    call handle_result(stat, istat)
  end subroutine set_xtol_rel

  function get_xtol_rel(self) result(tol)
    class(nlopt_opt), intent(inout) :: self
    real(wp) :: tol

    tol = nlopt_get_xtol_rel(self%handle)
  end function get_xtol_rel

  subroutine set_xtol_abs1(self, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(in) :: tol
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_xtol_abs1(self%handle, real(tol, c_double))
    call handle_result(stat, istat)
  end subroutine set_xtol_abs1

  subroutine set_xtol_abs(self, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    real(c_double), intent(in) :: tol(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_xtol_abs(self%handle, tol)
    call handle_result(stat, istat)
  end subroutine set_xtol_abs

  subroutine get_xtol_abs(self, tol, stat)
    class(nlopt_opt), intent(inout) :: self
    real(c_double), intent(inout) :: tol(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_get_xtol_abs(self%handle, tol)
    call handle_result(stat, istat)
  end subroutine get_xtol_abs

  subroutine set_x_weights1(self, w, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(in) :: w
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_x_weights1(self%handle, real(w, c_double))
    call handle_result(stat, istat)
  end subroutine set_x_weights1

  subroutine set_x_weights(self, w, stat)
    class(nlopt_opt), intent(inout) :: self
    real(c_double), intent(in) :: w(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_x_weights(self%handle, w)
    call handle_result(stat, istat)
  end subroutine set_x_weights

  subroutine get_x_weights(self, w, stat)
    class(nlopt_opt), intent(inout) :: self
    real(c_double), intent(inout) :: w(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_get_x_weights(self%handle, w)
    call handle_result(stat, istat)
  end subroutine get_x_weights


  subroutine set_maxeval(self, maxeval, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(in) :: maxeval
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_maxeval(self%handle, int(maxeval, c_int))
    call handle_result(stat, istat)
  end subroutine set_maxeval

  function get_maxeval(self) result(maxeval)
    class(nlopt_opt), intent(inout) :: self
    integer(ik) :: maxeval

    maxeval = nlopt_get_maxeval(self%handle)
  end function get_maxeval


  function get_numevals(self) result(numevals)
    class(nlopt_opt), intent(inout) :: self
    integer(ik) :: numevals

    numevals = nlopt_get_numevals(self%handle)
  end function get_numevals


  subroutine set_maxtime(self, maxtime, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(in) :: maxtime
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_maxtime(self%handle, real(maxtime, c_double))
    call handle_result(stat, istat)
  end subroutine set_maxtime

  function get_maxtime(self) result(maxtime)
    class(nlopt_opt), intent(inout) :: self
    real(wp) :: maxtime

    maxtime = nlopt_get_maxtime(self%handle)
  end function get_maxtime


  subroutine force_stop(self, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_force_stop(self%handle)
    call handle_result(stat, istat)
  end subroutine force_stop

  subroutine set_force_stop(self, val, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(in) :: val
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_force_stop(self%handle, int(val, c_int))
    call handle_result(stat, istat)
  end subroutine set_force_stop

  function get_force_stop(self) result(val)
    class(nlopt_opt), intent(inout) :: self
    integer(ik) :: val

    val = nlopt_get_force_stop(self%handle)
  end function get_force_stop


  ! /* more algorithm-specific parameters */
  subroutine set_local_optimizer(self, other, stat)
    class(nlopt_opt), intent(inout) :: self
    type(nlopt_opt), intent(inout) :: other
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_local_optimizer(self%handle, other%handle)
    call handle_result(stat, istat)
  end subroutine set_local_optimizer


  subroutine set_population(self, pop, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(in) :: pop
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_population(self%handle, int(pop, c_int))
    call handle_result(stat, istat)
  end subroutine set_population

  function get_population(self) result(pop)
    class(nlopt_opt), intent(inout) :: self
    integer(ik) :: pop

    pop = nlopt_get_population(self%handle)
  end function get_population


  subroutine set_vector_storage(self, dim, stat)
    class(nlopt_opt), intent(inout) :: self
    integer(ik), intent(in) :: dim
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_vector_storage(self%handle, int(dim, c_int))
    call handle_result(stat, istat)
  end subroutine set_vector_storage

  function get_vector_storage(self) result(dim)
    class(nlopt_opt), intent(inout) :: self
    integer(ik) :: dim

    dim = nlopt_get_vector_storage(self%handle)
  end function get_vector_storage


  subroutine set_default_initial_step(self, x, stat)
    class(nlopt_opt), intent(inout) :: self
    real(c_double), intent(in) :: x(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_default_initial_step(self%handle, x)
    call handle_result(stat, istat)
  end subroutine set_default_initial_step

  subroutine set_initial_step(self, dx, stat)
    class(nlopt_opt), intent(inout) :: self
    real(c_double), intent(in) :: dx(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_initial_step(self%handle, dx)
    call handle_result(stat, istat)
  end subroutine set_initial_step

  subroutine set_initial_step1(self, dx, stat)
    class(nlopt_opt), intent(inout) :: self
    real(wp), intent(in) :: dx
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_set_initial_step1(self%handle, real(dx, c_double))
    call handle_result(stat, istat)
  end subroutine set_initial_step1

  subroutine get_initial_step(self, x, dx, stat)
    class(nlopt_opt), intent(in) :: self
    real(c_double), intent(in) :: x(*)
    real(c_double), intent(in) :: dx(*)
    integer(ik), intent(out), optional :: stat

    integer(nlopt_result) :: istat

    istat = nlopt_get_initial_step(self%handle, x, dx)
    call handle_result(stat, istat)
  end subroutine get_initial_step


  pure function as_c_char(str) result(res)
    character(len=*), intent(in) :: str
    character(kind=c_char) :: res(len(str)+1)
    res = transfer(str // c_null_char, res)
  end function as_c_char

end module nlopt_wrap
