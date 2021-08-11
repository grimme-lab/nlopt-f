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

module nlopt_enum
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: &
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
    NLOPT_LD_CCSAQ, NLOPT_GN_ESCH, NLOPT_GN_AGS, NLOPT_NUM_ALGORITHMS, nlopt_algorithm

  public :: &
    NLOPT_FAILURE, NLOPT_INVALID_ARGS, NLOPT_OUT_OF_MEMORY, NLOPT_ROUNDOFF_LIMITED, &
    NLOPT_FORCED_STOP, NLOPT_SUCCESS, NLOPT_STOPVAL_REACHED, NLOPT_FTOL_REACHED, &
    NLOPT_XTOL_REACHED, NLOPT_MAXEVAL_REACHED, NLOPT_MAXTIME_REACHED, NLOPT_NUM_RESULTS, &
    nlopt_result

  public :: &
    result_to_string, result_from_string, algorithm_to_string, algorithm_from_string, &
    algorithm_name

  enum, bind(c)
    enumerator :: &
      NLOPT_GN_DIRECT = 0, &
      NLOPT_GN_DIRECT_L = 1, &
      NLOPT_GN_DIRECT_L_RAND = 2, &
      NLOPT_GN_DIRECT_NOSCAL = 3, &
      NLOPT_GN_DIRECT_L_NOSCAL = 4, &
      NLOPT_GN_DIRECT_L_RAND_NOSCAL = 5, &
      NLOPT_GN_ORIG_DIRECT = 6, &
      NLOPT_GN_ORIG_DIRECT_L = 7, &
      NLOPT_GD_STOGO = 8, &
      NLOPT_GD_STOGO_RAND = 9, &
      NLOPT_LD_LBFGS_NOCEDAL = 10, &
      NLOPT_LD_LBFGS = 11, &
      NLOPT_LN_PRAXIS = 12, &
      NLOPT_LD_VAR1 = 13, &
      NLOPT_LD_VAR2 = 14, &
      NLOPT_LD_TNEWTON = 15, &
      NLOPT_LD_TNEWTON_RESTART = 16, &
      NLOPT_LD_TNEWTON_PRECOND = 17, &
      NLOPT_LD_TNEWTON_PRECOND_RESTART = 18, &
      NLOPT_GN_CRS2_LM = 19, &
      NLOPT_GN_MLSL = 20, &
      NLOPT_GD_MLSL = 21, &
      NLOPT_GN_MLSL_LDS = 22, &
      NLOPT_GD_MLSL_LDS = 23, &
      NLOPT_LD_MMA = 24, &
      NLOPT_LN_COBYLA = 25, &
      NLOPT_LN_NEWUOA = 26, &
      NLOPT_LN_NEWUOA_BOUND = 27, &
      NLOPT_LN_NELDERMEAD = 28, &
      NLOPT_LN_SBPLX = 29, &
      NLOPT_LN_AUGLAG = 30, &
      NLOPT_LD_AUGLAG = 31, &
      NLOPT_LN_AUGLAG_EQ = 32, &
      NLOPT_LD_AUGLAG_EQ = 33, &
      NLOPT_LN_BOBYQA = 34, &
      NLOPT_GN_ISRES = 35, &
      NLOPT_AUGLAG = 36, &
      NLOPT_AUGLAG_EQ = 37, &
      NLOPT_G_MLSL = 38, &
      NLOPT_G_MLSL_LDS = 39, &
      NLOPT_LD_SLSQP = 40, &
      NLOPT_LD_CCSAQ = 41, &
      NLOPT_GN_ESCH = 42, &
      NLOPT_GN_AGS = 43, &
      NLOPT_NUM_ALGORITHMS
  end enum
  integer, parameter :: nlopt_algorithm = c_int

  enum, bind(c)
    enumerator :: &
      NLOPT_FAILURE = -1, &
      NLOPT_INVALID_ARGS = -2, &
      NLOPT_OUT_OF_MEMORY = -3, &
      NLOPT_ROUNDOFF_LIMITED = -4, &
      NLOPT_FORCED_STOP = -5, &
      NLOPT_SUCCESS = 1, &
      NLOPT_STOPVAL_REACHED = 2, &
      NLOPT_FTOL_REACHED = 3, &
      NLOPT_XTOL_REACHED = 4, &
      NLOPT_MAXEVAL_REACHED = 5, &
      NLOPT_MAXTIME_REACHED = 6, &
      NLOPT_NUM_RESULTS
  end enum
  integer, parameter :: nlopt_result = c_int

  interface
    ! extern const char *
    ! nlopt_algorithm_name(nlopt_algorithm algorithm);
    pure function nlopt_algorithm_name(algorithm) result(name) bind(c)
      import :: nlopt_algorithm, c_ptr
      implicit none
      integer(nlopt_algorithm), value :: algorithm
      type(c_ptr) :: name
    end function nlopt_algorithm_name

#if NLOPT_VERSION >= 20602
    ! extern const char *
    ! nlopt_algorithm_to_string(nlopt_algorithm algorithm);
    pure function nlopt_algorithm_to_string(algorithm) result(name) bind(c)
      import :: nlopt_algorithm, c_ptr
      implicit none
      integer(nlopt_algorithm), value :: algorithm
      type(c_ptr) :: name
    end function nlopt_algorithm_to_string

    ! extern nlopt_algorithm
    ! nlopt_algorithm_from_string(const char *name);
    pure function nlopt_algorithm_from_string(name) result(algorithm) bind(c)
      import :: nlopt_algorithm, c_char
      implicit none
      character(1, c_char), intent(in) :: name(*)
      integer(nlopt_algorithm) :: algorithm
    end function nlopt_algorithm_from_string

    ! extern const char *
    ! nlopt_result_to_string(nlopt_result stat);
    pure function nlopt_result_to_string(stat) result(name) bind(c)
      import :: nlopt_result, c_ptr
      implicit none
      integer(nlopt_result), value :: stat
      type(c_ptr) :: name
    end function nlopt_result_to_string

    ! extern nlopt_result
    ! nlopt_result_from_string(const char *name);
    pure function nlopt_result_from_string(name) result(stat) bind(c)
      import :: nlopt_result, c_char
      implicit none
      character(1, c_char), intent(in) :: name(*)
      integer(nlopt_result) :: stat
    end function nlopt_result_from_string
#endif
  end interface

  interface
    pure function strlen(str) result(len) bind(c)
      import :: c_ptr, c_int
      type(c_ptr), value :: str
      integer(c_int) :: len
    end function strlen
  end interface

contains

  function algorithm_name(algorithm) result(name)
    integer(nlopt_algorithm), intent(in) :: algorithm
    character(len=:, kind=c_char), allocatable :: name

    character(len=1, kind=c_char), pointer :: tmp
    name = arr_to_str(from_c_char(nlopt_algorithm_name(algorithm)))
  end function algorithm_name

  function result_to_string(stat) result(name)
    integer(nlopt_result), intent(in) :: stat
    character(len=:, kind=c_char), allocatable :: name

#if NLOPT_VERSION >= 20602
    name = arr_to_str(from_c_char(nlopt_result_to_string(stat)))
#else
    select case(stat)
    case default; name = "FAILURE"
    case(NLOPT_INVALID_ARGS); name = "INVALID_ARGS"
    case(NLOPT_OUT_OF_MEMORY); name = "OUT_OF_MEMORY"
    case(NLOPT_ROUNDOFF_LIMITED); name = "ROUNDOFF_LIMITED"
    case(NLOPT_FORCED_STOP); name = "FORCED_STOP"
    case(NLOPT_SUCCESS); name = "SUCCESS"
    case(NLOPT_STOPVAL_REACHED); name = "STOPVAL_REACHED"
    case(NLOPT_FTOL_REACHED); name = "FTOL_REACHED"
    case(NLOPT_XTOL_REACHED); name = "XTOL_REACHED"
    case(NLOPT_MAXEVAL_REACHED); name = "MAXEVAL_REACHED"
    case(NLOPT_MAXTIME_REACHED); name = "MAXTIME_REACHED"
    end select
#endif
  end function result_to_string

  function result_from_string(name) result(stat)
    character(len=*), intent(in) :: name
    integer(nlopt_result) :: stat

#if NLOPT_VERSION >= 20602
    stat = nlopt_result_from_string(as_c_char(name))
#else
    select case(name)
    case default; stat = NLOPT_FAILURE
    case("INVALID_ARGS"); stat = NLOPT_INVALID_ARGS
    case("OUT_OF_MEMORY"); stat = NLOPT_OUT_OF_MEMORY
    case("ROUNDOFF_LIMITED"); stat = NLOPT_ROUNDOFF_LIMITED
    case("FORCED_STOP"); stat = NLOPT_FORCED_STOP
    case("SUCCESS"); stat = NLOPT_SUCCESS
    case("STOPVAL_REACHED"); stat = NLOPT_STOPVAL_REACHED
    case("FTOL_REACHED"); stat = NLOPT_FTOL_REACHED
    case("XTOL_REACHED"); stat = NLOPT_XTOL_REACHED
    case("MAXEVAL_REACHED"); stat = NLOPT_MAXEVAL_REACHED
    case("MAXTIME_REACHED"); stat = NLOPT_MAXTIME_REACHED
    end select
#endif
  end function result_from_string

  function algorithm_to_string(algorithm) result(name)
    integer(nlopt_algorithm), intent(in) :: algorithm
    character(len=:, kind=c_char), allocatable :: name

#if NLOPT_VERSION >= 20602
    name = arr_to_str(from_c_char(nlopt_algorithm_to_string(algorithm)))
#else
    select case(algorithm)
    case(NLOPT_GN_DIRECT); name = "GN_DIRECT"
    case(NLOPT_GN_DIRECT_L); name = "GN_DIRECT_L"
    case(NLOPT_GN_DIRECT_L_RAND); name = "GN_DIRECT_L_RAND"
    case(NLOPT_GN_DIRECT_NOSCAL); name = "GN_DIRECT_NOSCAL"
    case(NLOPT_GN_DIRECT_L_NOSCAL); name = "GN_DIRECT_L_NOSCAL"
    case(NLOPT_GN_DIRECT_L_RAND_NOSCAL); name = "GN_DIRECT_L_RAND_NOSCAL"
    case(NLOPT_GN_ORIG_DIRECT); name = "GN_ORIG_DIRECT"
    case(NLOPT_GN_ORIG_DIRECT_L); name = "GN_ORIG_DIRECT_L"
    case(NLOPT_GD_STOGO); name = "GD_STOGO"
    case(NLOPT_GD_STOGO_RAND); name = "GD_STOGO_RAND"
    case(NLOPT_LD_LBFGS_NOCEDAL); name = "LD_LBFGS_NOCEDAL"
    case(NLOPT_LD_LBFGS); name = "LD_LBFGS"
    case(NLOPT_LN_PRAXIS); name = "LN_PRAXIS"
    case(NLOPT_LD_VAR1); name = "LD_VAR1"
    case(NLOPT_LD_VAR2); name = "LD_VAR2"
    case(NLOPT_LD_TNEWTON); name = "LD_TNEWTON"
    case(NLOPT_LD_TNEWTON_RESTART); name = "LD_TNEWTON_RESTART"
    case(NLOPT_LD_TNEWTON_PRECOND); name = "LD_TNEWTON_PRECOND"
    case(NLOPT_LD_TNEWTON_PRECOND_RESTART); name = "LD_TNEWTON_PRECOND_RESTART"
    case(NLOPT_GN_CRS2_LM); name = "GN_CRS2_LM"
    case(NLOPT_GN_MLSL); name = "GN_MLSL"
    case(NLOPT_GD_MLSL); name = "GD_MLSL"
    case(NLOPT_GN_MLSL_LDS); name = "GN_MLSL_LDS"
    case(NLOPT_GD_MLSL_LDS); name = "GD_MLSL_LDS"
    case(NLOPT_LD_MMA); name = "LD_MMA"
    case(NLOPT_LN_COBYLA); name = "LN_COBYLA"
    case(NLOPT_LN_NEWUOA); name = "LN_NEWUOA"
    case(NLOPT_LN_NEWUOA_BOUND); name = "LN_NEWUOA_BOUND"
    case(NLOPT_LN_NELDERMEAD); name = "LN_NELDERMEAD"
    case(NLOPT_LN_SBPLX); name = "LN_SBPLX"
    case(NLOPT_LN_AUGLAG); name = "LN_AUGLAG"
    case(NLOPT_LD_AUGLAG); name = "LD_AUGLAG"
    case(NLOPT_LN_AUGLAG_EQ); name = "LN_AUGLAG_EQ"
    case(NLOPT_LD_AUGLAG_EQ); name = "LD_AUGLAG_EQ"
    case(NLOPT_LN_BOBYQA); name = "LN_BOBYQA"
    case(NLOPT_GN_ISRES); name = "GN_ISRES"
    case(NLOPT_AUGLAG); name = "AUGLAG"
    case(NLOPT_AUGLAG_EQ); name = "AUGLAG_EQ"
    case(NLOPT_G_MLSL); name = "G_MLSL"
    case(NLOPT_G_MLSL_LDS); name = "G_MLSL_LDS"
    case(NLOPT_LD_SLSQP); name = "LD_SLSQP"
    case(NLOPT_LD_CCSAQ); name = "LD_CCSAQ"
    case(NLOPT_GN_ESCH); name = "GN_ESCH"
    case(NLOPT_GN_AGS); name = "GN_AGS"
    case default; name = ""
    end select
#endif
  end function algorithm_to_string

  function algorithm_from_string(name) result(algorithm)
    character(len=*), intent(in) :: name
    integer(nlopt_algorithm) :: algorithm
#if NLOPT_VERSION >= 20602
    algorithm = nlopt_algorithm_from_string(as_c_char(name))
#else
    select case(name)
    case("GN_DIRECT"); algorithm = NLOPT_GN_DIRECT
    case("GN_DIRECT_L"); algorithm = NLOPT_GN_DIRECT_L
    case("GN_DIRECT_L_RAND"); algorithm = NLOPT_GN_DIRECT_L_RAND
    case("GN_DIRECT_NOSCAL"); algorithm = NLOPT_GN_DIRECT_NOSCAL
    case("GN_DIRECT_L_NOSCAL"); algorithm = NLOPT_GN_DIRECT_L_NOSCAL
    case("GN_DIRECT_L_RAND_NOSCAL"); algorithm = NLOPT_GN_DIRECT_L_RAND_NOSCAL
    case("GN_ORIG_DIRECT"); algorithm = NLOPT_GN_ORIG_DIRECT
    case("GN_ORIG_DIRECT_L"); algorithm = NLOPT_GN_ORIG_DIRECT_L
    case("GD_STOGO"); algorithm = NLOPT_GD_STOGO
    case("GD_STOGO_RAND"); algorithm = NLOPT_GD_STOGO_RAND
    case("LD_LBFGS_NOCEDAL"); algorithm = NLOPT_LD_LBFGS_NOCEDAL
    case("LD_LBFGS"); algorithm = NLOPT_LD_LBFGS
    case("LN_PRAXIS"); algorithm = NLOPT_LN_PRAXIS
    case("LD_VAR1"); algorithm = NLOPT_LD_VAR1
    case("LD_VAR2"); algorithm = NLOPT_LD_VAR2
    case("LD_TNEWTON"); algorithm = NLOPT_LD_TNEWTON
    case("LD_TNEWTON_RESTART"); algorithm = NLOPT_LD_TNEWTON_RESTART
    case("LD_TNEWTON_PRECOND"); algorithm = NLOPT_LD_TNEWTON_PRECOND
    case("LD_TNEWTON_PRECOND_RESTART"); algorithm = NLOPT_LD_TNEWTON_PRECOND_RESTART
    case("GN_CRS2_LM"); algorithm = NLOPT_GN_CRS2_LM
    case("GN_MLSL"); algorithm = NLOPT_GN_MLSL
    case("GD_MLSL"); algorithm = NLOPT_GD_MLSL
    case("GN_MLSL_LDS"); algorithm = NLOPT_GN_MLSL_LDS
    case("GD_MLSL_LDS"); algorithm = NLOPT_GD_MLSL_LDS
    case("LD_MMA"); algorithm = NLOPT_LD_MMA
    case("LN_COBYLA"); algorithm = NLOPT_LN_COBYLA
    case("LN_NEWUOA"); algorithm = NLOPT_LN_NEWUOA
    case("LN_NEWUOA_BOUND"); algorithm = NLOPT_LN_NEWUOA_BOUND
    case("LN_NELDERMEAD"); algorithm = NLOPT_LN_NELDERMEAD
    case("LN_SBPLX"); algorithm = NLOPT_LN_SBPLX
    case("LN_AUGLAG"); algorithm = NLOPT_LN_AUGLAG
    case("LD_AUGLAG"); algorithm = NLOPT_LD_AUGLAG
    case("LN_AUGLAG_EQ"); algorithm = NLOPT_LN_AUGLAG_EQ
    case("LD_AUGLAG_EQ"); algorithm = NLOPT_LD_AUGLAG_EQ
    case("LN_BOBYQA"); algorithm = NLOPT_LN_BOBYQA
    case("GN_ISRES"); algorithm = NLOPT_GN_ISRES
    case("AUGLAG"); algorithm = NLOPT_AUGLAG
    case("AUGLAG_EQ"); algorithm = NLOPT_AUGLAG_EQ
    case("G_MLSL"); algorithm = NLOPT_G_MLSL
    case("G_MLSL_LDS"); algorithm = NLOPT_G_MLSL_LDS
    case("LD_SLSQP"); algorithm = NLOPT_LD_SLSQP
    case("LD_CCSAQ"); algorithm = NLOPT_LD_CCSAQ
    case("GN_ESCH"); algorithm = NLOPT_GN_ESCH
    case("GN_AGS"); algorithm = NLOPT_GN_AGS
    case default; algorithm = -1
    end select
#endif
  end function algorithm_from_string

  function from_c_char(ptr) result(tmp)
    type(c_ptr), value :: ptr
    character(len=1, kind=c_char), pointer :: tmp(:)
    character(len=1, kind=c_char), target :: dummy(0)
    if (c_associated(ptr)) then
      call c_f_pointer(ptr, tmp, [strlen(ptr)])
    else
      tmp => dummy
    end if
  end function from_c_char

  pure function arr_to_str(arr) result(str)
    character(len=1, kind=c_char), intent(in) :: arr(:)
    character(len=size(arr), kind=c_char) :: str
    str = transfer(arr, str)
  end function arr_to_str

  pure function as_c_char(str) result(res)
    character(len=*), intent(in) :: str
    character(kind=c_char) :: res(len(str)+1)
    res = transfer(str // c_null_char, res)
  end function as_c_char

end module nlopt_enum
