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

#include "nlopt-wrap.F90"

module nlopt_interface
  use, intrinsic :: iso_c_binding, only : c_ptr, c_funptr, c_int, c_double, c_char, c_long, &
     c_null_ptr
  use nlopt_enum, only : nlopt_algorithm, nlopt_result, NLOPT_FAILURE
  implicit none
  private

  public :: &
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

  interface
    ! extern void
    ! nlopt_srand(unsigned long seed);
    subroutine nlopt_srand(seed) bind(c)
      import :: c_long
      implicit none
      integer(c_long), value :: seed
    end subroutine nlopt_srand

    ! extern void
    ! nlopt_srand_time(void);
    subroutine nlopt_srand_time() bind(c)
    end subroutine nlopt_srand_time


    ! extern void
    ! nlopt_version(int *major, int *minor, int *bugfix);
    subroutine nlopt_version(major, minor, bugfix) bind(c)
      import :: c_int
      implicit none
      integer(c_int), intent(out) :: major, minor, bugfix
    end subroutine nlopt_version

    ! extern nlopt_opt
    ! nlopt_create(nlopt_algorithm algorithm, unsigned n);
    function nlopt_create(algorithm, n) result(opt) bind(c)
      import :: nlopt_algorithm, c_int, c_ptr
      implicit none
      integer(nlopt_algorithm), value :: algorithm
      integer(c_int), value :: n
      type(c_ptr) :: opt
    end function

    ! extern void
    ! nlopt_destroy(nlopt_opt opt);
    subroutine nlopt_destroy(opt) bind(c)
      import :: c_ptr
      implicit none
      type(c_ptr), value :: opt
    end subroutine nlopt_destroy

    ! extern nlopt_opt
    ! nlopt_copy(const nlopt_opt opt);
    function nlopt_copy(opt) result(new) bind(c)
      import :: c_ptr
      implicit none
      type(c_ptr), value :: opt
      type(c_ptr) :: new
    end function nlopt_copy


    ! extern nlopt_result
    ! nlopt_optimize(nlopt_opt opt, double *x, double *opt_f);
    function nlopt_optimize(opt, x, opt_f) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(inout) :: x(*)
      real(c_double), intent(inout) :: opt_f
      integer(nlopt_result) :: stat
    end function nlopt_optimize

    ! extern nlopt_result
    ! nlopt_set_min_objective(nlopt_opt opt, nlopt_func f, void *f_data);
    function nlopt_set_min_objective(opt, f, f_data) result(stat) bind(c)
      import :: c_ptr, c_funptr, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      type(c_funptr), value :: f
      type(c_ptr), value :: f_data
      integer(nlopt_result) :: stat
    end function nlopt_set_min_objective
    ! extern nlopt_result
    ! nlopt_set_max_objective(nlopt_opt opt, nlopt_func f, void *f_data);
    function nlopt_set_max_objective(opt, f, f_data) result(stat) bind(c)
      import :: c_ptr, c_funptr, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      type(c_funptr), value :: f
      type(c_ptr), value :: f_data
      integer(nlopt_result) :: stat
    end function nlopt_set_max_objective

#if NLOPT_VERSION >= 20300
    ! extern nlopt_result
    ! nlopt_set_precond_min_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre, void *f_data);
    function nlopt_set_precond_min_objective(opt, f, pre, f_data) result(stat) bind(c)
      import :: c_ptr, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      type(c_ptr), value :: f
      type(c_ptr), value :: pre
      type(c_ptr), value :: f_data
      integer(nlopt_result) :: stat
    end function nlopt_set_precond_min_objective
    ! extern nlopt_result
    ! nlopt_set_precond_max_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre, void *f_data);
    function nlopt_set_precond_max_objective(opt, f, pre, f_data) result(stat) bind(c)
      import :: c_ptr, c_funptr, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      type(c_funptr), value :: f
      type(c_funptr), value :: pre
      type(c_ptr), value :: f_data
      integer(nlopt_result) :: stat
    end function nlopt_set_precond_max_objective
#endif

    ! extern nlopt_algorithm
    ! nlopt_get_algorithm(const nlopt_opt opt);
    function nlopt_get_algorithm(opt) result(algorithm) bind(c)
      import :: c_ptr, nlopt_algorithm
      implicit none
      type(c_ptr), value :: opt
      integer(nlopt_algorithm) :: algorithm
    end function nlopt_get_algorithm
    ! extern unsigned
    ! nlopt_get_dimension(const nlopt_opt opt);
    function nlopt_get_dimension(opt) result(n) bind(c)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: opt
      integer(c_int) :: n
    end function nlopt_get_dimension

#if NLOPT_VERSION >= 20500
    ! extern const char *
    ! nlopt_get_errmsg(nlopt_opt opt);
    function nlopt_get_errmsg(opt) result(errmsg) bind(c)
      import :: c_ptr
      implicit none
      type(c_ptr), value :: opt
      type(c_ptr) :: errmsg
    end function nlopt_get_errmsg
#endif

#if NLOPT_VERSION >= 20700
    ! /* generic algorithm parameters: */
    ! extern nlopt_result
    ! nlopt_set_param(nlopt_opt opt, const char *name, double val);
    function nlopt_set_param(opt, name, val) result(stat) bind(c)
      import :: c_ptr, c_char, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      character(1, c_char), intent(in) :: name(*)
      real(c_double), value :: val
      integer(nlopt_result) :: stat
    end function nlopt_set_param
    ! extern double
    ! nlopt_get_param(const nlopt_opt opt, const char *name, double defaultval);
    function nlopt_get_param(opt, name, defaultval) result(val) bind(c)
      import :: c_ptr, c_char, c_double
      implicit none
      type(c_ptr), value :: opt
      character(1, c_char), intent(in) :: name(*)
      real(c_double), value :: defaultval
      real(c_double) :: val
    end function nlopt_get_param
    ! extern int
    ! nlopt_has_param(const nlopt_opt opt, const char *name);
    function nlopt_has_param(opt, name) result(stat) bind(c)
      import :: c_ptr, c_char, c_int
      implicit none
      type(c_ptr), value :: opt
      character(1, c_char), intent(in) :: name(*)
      integer(c_int) :: stat
    end function nlopt_has_param
    ! extern unsigned
    ! nlopt_num_params(const nlopt_opt opt);
    function nlopt_num_params(opt) result(n) bind(c)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: opt
      integer(c_int) :: n
    end function nlopt_num_params
    ! extern const char *
    ! nlopt_nth_param(const nlopt_opt opt, unsigned n);
    function nlopt_nth_param(opt, n) result(name) bind(c)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: opt
      integer(c_int), value :: n
      type(c_ptr) :: name
    end function nlopt_nth_param
#endif

    ! /* constraints: */

    ! extern nlopt_result
    ! nlopt_set_lower_bounds(nlopt_opt opt, const double *lb);
    function nlopt_set_lower_bounds(opt, lb) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(in) :: lb(*)
      integer(nlopt_result) :: stat
    end function nlopt_set_lower_bounds
    ! extern nlopt_result nlopt_set_lower_bounds1(nlopt_opt opt, double lb);
    function nlopt_set_lower_bounds1(opt, lb) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), value :: lb
      integer(nlopt_result) :: stat
    end function nlopt_set_lower_bounds1
#if NLOPT_VERSION >= 20600
    ! extern nlopt_result
    ! nlopt_set_lower_bound(nlopt_opt opt, int i, double lb);
    function nlopt_set_lower_bound(opt, i, lb) result(stat) bind(c)
      import :: c_ptr, c_int, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(c_int), value :: i
      real(c_double), value :: lb
      integer(nlopt_result) :: stat
    end function nlopt_set_lower_bound
#endif
    ! extern nlopt_result
    ! nlopt_get_lower_bounds(const nlopt_opt opt, double *lb);
    function nlopt_get_lower_bounds(opt, lb) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(inout) :: lb(*)
      integer(nlopt_result) :: stat
    end function nlopt_get_lower_bounds
    ! extern nlopt_result
    ! nlopt_set_upper_bounds(nlopt_opt opt, const double *ub);
    function nlopt_set_upper_bounds(opt, ub) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(in) :: ub(*)
      integer(nlopt_result) :: stat
    end function nlopt_set_upper_bounds
    ! extern nlopt_result nlopt_set_upper_bounds1(nlopt_opt opt, double ub);
    function nlopt_set_upper_bounds1(opt, ub) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), value :: ub
      integer(nlopt_result) :: stat
    end function nlopt_set_upper_bounds1
#if NLOPT_VERSION >= 20600
    ! extern nlopt_result
    ! nlopt_set_upper_bound(nlopt_opt opt, int i, double ub);
    function nlopt_set_upper_bound(opt, i, ub) result(stat) bind(c)
      import :: c_ptr, c_int, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(c_int), value :: i
      real(c_double), value :: ub
      integer(nlopt_result) :: stat
    end function nlopt_set_upper_bound
#endif
    ! extern nlopt_result
    ! nlopt_get_upper_bounds(const nlopt_opt opt, double *ub);
    function nlopt_get_upper_bounds(opt, ub) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(inout) :: ub(*)
      integer(nlopt_result) :: stat
    end function nlopt_get_upper_bounds

    ! extern nlopt_result
    ! nlopt_remove_inequality_constraints(nlopt_opt opt);
    function nlopt_remove_inequality_constraints(opt) result(stat) bind(c)
      import :: c_ptr, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(nlopt_result) :: stat
    end function nlopt_remove_inequality_constraints
    ! extern nlopt_result
    ! nlopt_add_inequality_constraint(nlopt_opt opt, nlopt_func fc, void *fc_data, double tol);
    function nlopt_add_inequality_constraint(opt, fc, fc_data, tol) result(stat) bind(c)
      import :: c_ptr, c_funptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      type(c_funptr), value :: fc
      type(c_ptr), value :: fc_data
      real(c_double), value :: tol
      integer(nlopt_result) :: stat
    end function nlopt_add_inequality_constraint
#if NLOPT_VERSION >= 20300
    ! extern nlopt_result
    ! nlopt_add_precond_inequality_constraint(nlopt_opt opt, nlopt_func fc, nlopt_precond pre, void *fc_data, double tol);
    function nlopt_add_precond_inequality_constraint(opt, fc, pre, fc_data, tol) &
        & result(stat) bind(c)
      import :: c_ptr, c_funptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      type(c_funptr), value :: fc
      type(c_funptr), value :: pre
      type(c_ptr), value :: fc_data
      real(c_double), value :: tol
      integer(nlopt_result) :: stat
    end function nlopt_add_precond_inequality_constraint
#endif
#if NLOPT_VERSION >= 20100
    ! extern nlopt_result
    ! nlopt_add_inequality_mconstraint(nlopt_opt opt, unsigned m, nlopt_mfunc fc, void *fc_data, const double *tol);
    function nlopt_add_inequality_mconstraint(opt, m, fc, fc_data, tol) &
        & result(stat) bind(c)
      import :: c_ptr, c_funptr, c_int, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(c_int), value :: m
      type(c_funptr), value :: fc
      type(c_ptr), value :: fc_data
      real(c_double), intent(in) :: tol(*)
      integer(nlopt_result) :: stat
    end function nlopt_add_inequality_mconstraint
#endif

    ! extern nlopt_result
    ! nlopt_remove_equality_constraints(nlopt_opt opt);
    function nlopt_remove_equality_constraints(opt) result(stat) bind(c)
      import :: c_ptr, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(nlopt_result) :: stat
    end function nlopt_remove_equality_constraints
    ! extern nlopt_result
    ! nlopt_add_equality_constraint(nlopt_opt opt, nlopt_func h, void *h_data, double tol);
    function nlopt_add_equality_constraint(opt, h, h_data, tol) result(stat) bind(c)
      import :: c_ptr, c_funptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      type(c_funptr), value :: h
      type(c_ptr), value :: h_data
      real(c_double), value :: tol
      integer(nlopt_result) :: stat
    end function nlopt_add_equality_constraint
#if NLOPT_VERSION >= 20300
    ! extern nlopt_result
    ! nlopt_add_precond_equality_constraint(nlopt_opt opt, nlopt_func h, nlopt_precond pre, void *h_data, double tol);
    function nlopt_add_precond_equality_constraint(opt, h, pre, h_data, tol) &
        & result(stat) bind(c)
      import :: c_ptr, c_funptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      type(c_funptr), value :: h
      type(c_funptr), value :: pre
      type(c_ptr), value :: h_data
      real(c_double), value :: tol
      integer(nlopt_result) :: stat
    end function nlopt_add_precond_equality_constraint
#endif
#if NLOPT_VERSION > 20100
    ! extern nlopt_result
    ! nlopt_add_equality_mconstraint(nlopt_opt opt, unsigned m, nlopt_mfunc h, void *h_data, const double *tol);
    function nlopt_add_equality_mconstraint(opt, m, h, h_data, tol) result(stat) bind(c)
      import :: c_ptr, c_funptr, c_int, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(c_int), value :: m
      type(c_funptr), value :: h
      type(c_ptr), value :: h_data
      real(c_double), intent(in) :: tol(*)
      integer(nlopt_result) :: stat
    end function nlopt_add_equality_mconstraint
#endif

    ! /* stopping criteria: */

    ! extern nlopt_result
    ! nlopt_set_stopval(nlopt_opt opt, double stopval);
    function nlopt_set_stopval(opt, stopval) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), value :: stopval
      integer(nlopt_result) :: stat
    end function nlopt_set_stopval
    ! extern double
    ! nlopt_get_stopval(const nlopt_opt opt);
    function nlopt_get_stopval(opt) result(stopval) bind(c)
      import :: c_ptr, c_double
      implicit none
      type(c_ptr), value :: opt
      real(c_double) :: stopval
    end function nlopt_get_stopval

    ! extern nlopt_result
    ! nlopt_set_ftol_rel(nlopt_opt opt, double tol);
    function nlopt_set_ftol_rel(opt, tol) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), value :: tol
      integer(nlopt_result) :: stat
    end function nlopt_set_ftol_rel
    ! extern double
    ! nlopt_get_ftol_rel(const nlopt_opt opt);
    function nlopt_get_ftol_rel(opt) result(tol) bind(c)
      import :: c_ptr, c_double
      implicit none
      type(c_ptr), value :: opt
      real(c_double) :: tol
    end function nlopt_get_ftol_rel
    ! extern nlopt_result
    ! nlopt_set_ftol_abs(nlopt_opt opt, double tol);
    function nlopt_set_ftol_abs(opt, tol) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), value :: tol
      integer(nlopt_result) :: stat
    end function nlopt_set_ftol_abs
    ! extern double
    ! nlopt_get_ftol_abs(const nlopt_opt opt);
    function nlopt_get_ftol_abs(opt) result(tol) bind(c)
      import :: c_ptr, c_double
      implicit none
      type(c_ptr), value :: opt
      real(c_double) :: tol
    end function nlopt_get_ftol_abs

    ! extern nlopt_result
    ! nlopt_set_xtol_rel(nlopt_opt opt, double tol);
    function nlopt_set_xtol_rel(opt, tol) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), value :: tol
      integer(nlopt_result) :: stat
    end function nlopt_set_xtol_rel
    ! extern double
    ! nlopt_get_xtol_rel(const nlopt_opt opt);
    function nlopt_get_xtol_rel(opt) result(tol) bind(c)
      import :: c_ptr, c_double
      implicit none
      type(c_ptr), value :: opt
      real(c_double) :: tol
    end function nlopt_get_xtol_rel
    ! extern nlopt
    ! _result nlopt_set_xtol_abs1(nlopt_opt opt, double tol);
    function nlopt_set_xtol_abs1(opt, tol) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), value :: tol
      integer(nlopt_result) :: stat
    end function nlopt_set_xtol_abs1
    ! extern nlopt_result
    ! nlopt_set_xtol_abs(nlopt_opt opt, const double *tol);
    function nlopt_set_xtol_abs(opt, tol) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(in) :: tol(*)
      integer(nlopt_result) :: stat
    end function nlopt_set_xtol_abs
    ! extern nlopt_result
    ! nlopt_get_xtol_abs(const nlopt_opt opt, double *tol);
    function nlopt_get_xtol_abs(opt, tol) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(inout) :: tol(*)
      integer(nlopt_result) :: stat
    end function nlopt_get_xtol_abs
#if NLOPT_VERSION >= 20602
    ! extern nlopt_result nlopt_set_x_weights1(nlopt_opt opt, double w);
    function nlopt_set_x_weights1(opt, w) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), value :: w
      integer(nlopt_result) :: stat
    end function nlopt_set_x_weights1
    ! extern nlopt_result
    ! nlopt_set_x_weights(nlopt_opt opt, const double *w);
    function nlopt_set_x_weights(opt, w) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(in) :: w(*)
      integer(nlopt_result) :: stat
    end function nlopt_set_x_weights
    ! extern nlopt_result
    ! nlopt_get_x_weights(const nlopt_opt opt, double *w);
    function nlopt_get_x_weights(opt, w) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(inout) :: w(*)
      integer(nlopt_result) :: stat
    end function nlopt_get_x_weights
#endif

    ! extern nlopt_result
    ! nlopt_set_maxeval(nlopt_opt opt, int maxeval);
    function nlopt_set_maxeval(opt, maxeval) result(stat) bind(c)
      import :: c_ptr, c_int, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(c_int), value :: maxeval
      integer(nlopt_result) :: stat
    end function nlopt_set_maxeval
    ! extern int
    ! nlopt_get_maxeval(const nlopt_opt opt);
    function nlopt_get_maxeval(opt) result(maxeval) bind(c)
      import :: c_ptr, c_int, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(c_int) :: maxeval
    end function nlopt_get_maxeval

#if NLOPT_VERSION >= 20500
    ! extern int
    ! nlopt_get_numevals(const nlopt_opt opt);
    function nlopt_get_numevals(opt) result(numevals) bind(c)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: opt
      integer(c_int) :: numevals
    end function nlopt_get_numevals
#endif

    ! extern nlopt_result
    ! nlopt_set_maxtime(nlopt_opt opt, double maxtime);
    function nlopt_set_maxtime(opt, maxtime) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), value :: maxtime
      integer(nlopt_result) :: stat
    end function nlopt_set_maxtime
    ! extern double
    ! nlopt_get_maxtime(const nlopt_opt opt);
    function nlopt_get_maxtime(opt) result(maxtime) bind(c)
      import :: c_ptr, c_double
      implicit none
      type(c_ptr), value :: opt
      real(c_double) :: maxtime
    end function nlopt_get_maxtime

    ! extern nlopt_result
    ! nlopt_force_stop(nlopt_opt opt);
    function nlopt_force_stop(opt) result(stat) bind(c)
      import :: c_ptr, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(nlopt_result) :: stat
    end function nlopt_force_stop
    ! extern nlopt_result
    ! nlopt_set_force_stop(nlopt_opt opt, int val);
    function nlopt_set_force_stop(opt, val) result(stat) bind(c)
      import :: c_ptr, c_int, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(c_int), value :: val
      integer(nlopt_result) :: stat
    end function nlopt_set_force_stop
    ! extern int
    ! nlopt_get_force_stop(const nlopt_opt opt);
    function nlopt_get_force_stop(opt) result(val) bind(c)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: opt
      integer(c_int) :: val
    end function nlopt_get_force_stop

    ! /* more algorithm-specific parameters */

    ! extern nlopt_result
    ! nlopt_set_local_optimizer(nlopt_opt opt, const nlopt_opt local_opt);
    function nlopt_set_local_optimizer(opt, local_opt) result(stat) bind(c)
      import :: c_ptr, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      type(c_ptr), value :: local_opt
      integer(nlopt_result) :: stat
    end function nlopt_set_local_optimizer

    ! extern nlopt_result
    ! nlopt_set_population(nlopt_opt opt, unsigned pop);
    function nlopt_set_population(opt, pop) result(stat) bind(c)
      import :: c_ptr, c_int, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(c_int), value :: pop
      integer(nlopt_result) :: stat
    end function nlopt_set_population
    ! extern unsigned
    ! nlopt_get_population(const nlopt_opt opt);
    function nlopt_get_population(opt) result(pop) bind(c)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: opt
      integer(c_int) :: pop
    end function nlopt_get_population

#if NLOPT_VERSION >= 20202
    ! extern nlopt_result
    ! nlopt_set_vector_storage(nlopt_opt opt, unsigned dim);
    function nlopt_set_vector_storage(opt, dim) result(stat) bind(c)
      import :: c_ptr, c_int, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      integer(c_int), value :: dim
      integer(nlopt_result) :: stat
    end function nlopt_set_vector_storage
    ! extern unsigned
    ! nlopt_get_vector_storage(const nlopt_opt opt);
    function nlopt_get_vector_storage(opt) result(dim) bind(c)
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: opt
      integer(c_int) :: dim
    end function nlopt_get_vector_storage
#endif

    ! extern nlopt_result
    ! nlopt_set_default_initial_step(nlopt_opt opt, const double *x);
    function nlopt_set_default_initial_step(opt, x) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(in) :: x(*)
      integer(nlopt_result) :: stat
    end function nlopt_set_default_initial_step
    ! extern nlopt_result
    ! nlopt_set_initial_step(nlopt_opt opt, const double *dx);
    function nlopt_set_initial_step(opt, dx) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(in) :: dx(*)
      integer(nlopt_result) :: stat
    end function nlopt_set_initial_step
    ! extern nlopt_result nlopt_set_initial_step1(nlopt_opt opt, double dx);
    function nlopt_set_initial_step1(opt, dx) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), value :: dx
      integer(nlopt_result) :: stat
    end function nlopt_set_initial_step1
    ! extern nlopt_result
    ! nlopt_get_initial_step(const nlopt_opt opt, const double *x, double *dx);
    function nlopt_get_initial_step(opt, x, dx) result(stat) bind(c)
      import :: c_ptr, c_double, nlopt_result
      implicit none
      type(c_ptr), value :: opt
      real(c_double), intent(in) :: x(*)
      real(c_double), intent(in) :: dx(*)
      integer(nlopt_result) :: stat
    end function nlopt_get_initial_step
  end interface


contains

#if NLOPT_VERSION < 20700
  ! /* generic algorithm parameters: */
  ! extern nlopt_result
  ! nlopt_set_param(nlopt_opt opt, const char *name, double val);
  function nlopt_set_param(opt, name, val) result(stat)
    type(c_ptr), value :: opt
    character(1, c_char), intent(in) :: name(*)
    real(c_double), value :: val
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_set_param
  ! extern double
  ! nlopt_get_param(const nlopt_opt opt, const char *name, double defaultval);
  function nlopt_get_param(opt, name, defaultval) result(val)
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    type(c_ptr), value :: opt
    character(1, c_char), intent(in) :: name(*)
    real(c_double), value :: defaultval
    real(c_double) :: val
    val = ieee_value(val, ieee_quiet_nan)
  end function nlopt_get_param
  ! extern int
  ! nlopt_has_param(const nlopt_opt opt, const char *name);
  function nlopt_has_param(opt, name) result(stat)
    type(c_ptr), value :: opt
    character(1, c_char), intent(in) :: name(*)
    integer(c_int) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_has_param
  ! extern unsigned
  ! nlopt_num_params(const nlopt_opt opt);
  function nlopt_num_params(opt) result(n)
    type(c_ptr), value :: opt
    integer(c_int) :: n
    n = 0_c_int
  end function nlopt_num_params
  ! extern const char *
  ! nlopt_nth_param(const nlopt_opt opt, unsigned n);
  function nlopt_nth_param(opt, n) result(name)
    type(c_ptr), value :: opt
    integer(c_int), value :: n
    type(c_ptr) :: name
    name = c_null_ptr
  end function nlopt_nth_param
#endif

#if NLOPT_VERSION < 20602
  ! extern nlopt_result nlopt_set_x_weights1(nlopt_opt opt, double w);
  function nlopt_set_x_weights1(opt, w) result(stat)
    type(c_ptr), value :: opt
    real(c_double), value :: w
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_set_x_weights1
  ! extern nlopt_result
  ! nlopt_set_x_weights(nlopt_opt opt, const double *w);
  function nlopt_set_x_weights(opt, w) result(stat)
    type(c_ptr), value :: opt
    real(c_double), intent(in) :: w(*)
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_set_x_weights
  ! extern nlopt_result
  ! nlopt_get_x_weights(const nlopt_opt opt, double *w);
  function nlopt_get_x_weights(opt, w) result(stat)
    type(c_ptr), value :: opt
    real(c_double), intent(inout) :: w(*)
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_get_x_weights
#endif

#if NLOPT_VERSION < 20600
  ! extern nlopt_result
  ! nlopt_set_lower_bound(nlopt_opt opt, int i, double lb);
  function nlopt_set_lower_bound(opt, i, lb) result(stat)
    type(c_ptr), value :: opt
    integer(c_int), value :: i
    real(c_double), value :: lb
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_set_lower_bound

  ! extern nlopt_result
  ! nlopt_set_upper_bound(nlopt_opt opt, int i, double ub);
  function nlopt_set_upper_bound(opt, i, ub) result(stat)
    type(c_ptr), value :: opt
    integer(c_int), value :: i
    real(c_double), value :: ub
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_set_upper_bound
#endif


#if NLOPT_VERSION < 20500
  ! extern const char *
  ! nlopt_get_errmsg(nlopt_opt opt);
  function nlopt_get_errmsg(opt) result(errmsg)
    type(c_ptr), value :: opt
    type(c_ptr) :: errmsg
    errmsg = c_null_ptr
  end function nlopt_get_errmsg
#endif

#if NLOPT_VERSION < 20500
  ! extern int
  ! nlopt_get_numevals(const nlopt_opt opt);
  function nlopt_get_numevals(opt) result(numevals)
    type(c_ptr), value :: opt
    integer(c_int) :: numevals
    numevals = -1_c_int
  end function nlopt_get_numevals
#endif

#if NLOPT_VERSION < 20300
  ! extern nlopt_result
  ! nlopt_set_precond_min_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre, void *f_data);
  function nlopt_set_precond_min_objective(opt, f, pre, f_data) result(stat)
    type(c_ptr), value :: opt
    type(c_ptr), value :: f
    type(c_ptr), value :: pre
    type(c_ptr), value :: f_data
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_set_precond_min_objective

  ! extern nlopt_result
  ! nlopt_set_precond_max_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre, void *f_data);
  function nlopt_set_precond_max_objective(opt, f, pre, f_data) result(stat)
    type(c_ptr), value :: opt
    type(c_funptr), value :: f
    type(c_funptr), value :: pre
    type(c_ptr), value :: f_data
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_set_precond_max_objective

  ! extern nlopt_result
  ! nlopt_add_precond_inequality_constraint(nlopt_opt opt, nlopt_func fc, nlopt_precond pre, void *fc_data, double tol);
  function nlopt_add_precond_inequality_constraint(opt, fc, pre, fc_data, tol) &
      & result(stat)
    type(c_ptr), value :: opt
    type(c_funptr), value :: fc
    type(c_funptr), value :: pre
    type(c_ptr), value :: fc_data
    real(c_double), value :: tol
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_add_precond_inequality_constraint

  ! extern nlopt_result
  ! nlopt_add_precond_equality_constraint(nlopt_opt opt, nlopt_func h, nlopt_precond pre, void *h_data, double tol);
  function nlopt_add_precond_equality_constraint(opt, h, pre, h_data, tol) &
      & result(stat)
    type(c_ptr), value :: opt
    type(c_funptr), value :: h
    type(c_funptr), value :: pre
    type(c_ptr), value :: h_data
    real(c_double), value :: tol
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_add_precond_equality_constraint
#endif

#if NLOPT_VERSION < 20202
  ! extern nlopt_result
  ! nlopt_set_vector_storage(nlopt_opt opt, unsigned dim);
  function nlopt_set_vector_storage(opt, dim) result(stat)
    type(c_ptr), value :: opt
    integer(c_int), value :: dim
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_set_vector_storage
  ! extern unsigned
  ! nlopt_get_vector_storage(const nlopt_opt opt);
  function nlopt_get_vector_storage(opt) result(dim)
    type(c_ptr), value :: opt
    integer(c_int) :: dim
    dim = 0_c_int
  end function nlopt_get_vector_storage
#endif

#if NLOPT_VERSION < 20100
  ! extern nlopt_result
  ! nlopt_add_inequality_mconstraint(nlopt_opt opt, unsigned m, nlopt_mfunc fc, void *fc_data, const double *tol);
  function nlopt_add_inequality_mconstraint(opt, m, fc, fc_data, tol) &
      & result(stat)
    type(c_ptr), value :: opt
    integer(c_int), value :: m
    type(c_funptr), value :: fc
    type(c_ptr), value :: fc_data
    real(c_double), intent(in) :: tol(*)
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_add_inequality_mconstraint

  ! extern nlopt_result
  ! nlopt_add_equality_mconstraint(nlopt_opt opt, unsigned m, nlopt_mfunc h, void *h_data, const double *tol);
  function nlopt_add_equality_mconstraint(opt, m, h, h_data, tol) result(stat)
    type(c_ptr), value :: opt
    integer(c_int), value :: m
    type(c_funptr), value :: h
    type(c_ptr), value :: h_data
    real(c_double), intent(in) :: tol(*)
    integer(nlopt_result) :: stat
    stat = NLOPT_FAILURE
  end function nlopt_add_equality_mconstraint
#endif

end module nlopt_interface
