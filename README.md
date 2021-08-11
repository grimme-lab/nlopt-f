# Fortran bindings for the NLopt library

[![License](https://img.shields.io/badge/license-MIT%7CApache%202.0-blue)](LICENSE-Apache)
[![CI](https://github.com/grimme-lab/nlopt-f/actions/workflows/build.yml/badge.svg)](https://github.com/grimme-lab/nlopt-f/actions/workflows/build.yml)

Fortran bindings for the NLopt library.
While the NLopt library supports Fortran by using implicit interface calling conventions, those are not type-safe.
This project offers an alternative interface to the NLopt C-API.


## Installation

To build this project from the source code in this repository you need to have

- a Fortran compiler supporting Fortran 2008 (GCC 5 or newer or Intel Fortran)
- One of the supported build systems

  - [meson](https://mesonbuild.com) version 0.55 or newer
  - [Fortran package manager (fpm)](https://github.com/fortran-lang/fpm) version 0.2.0 or newer

- [nlopt](https://nlopt.readthedocs.io/en/latest/) version 2.5.0 or later


### Building with meson

Setup a build with

```
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable.
To compile the project run

```
meson compile -C _build
```

You can run the projects testsuite with

```
meson test -C _build --print-errorlogs
```

To include ``nlopt-f`` in your project add the following wrap file to your subprojects directory:

```ini
[wrap-git]
directory = nlopt-f
url = https://github.com/grimme-lab/nlopt-f
revision = head
```

You can retrieve the dependency from the wrap fallback with

```meson
nlopt_dep = dependency('nlopt-f', fallback: ['nlopt-f', 'nlopt_dep'])
```

and add it as dependency to your targets.


### Building with CMake

Alternatively, this project can be build with CMake (in this case ninja 1.10 or newer is required):

```
cmake -B _build -G Ninja
```

To compile the project with CMake run

```
cmake --build _build
```

You can run the project testsuite with

```
pushd _build && ctest && popd
```

Finally, install the project using (the install prefix can be customized with in ``-DCMAKE_INSTALL_PREFIX=/path/to/install`` in the first step)

```
cmake --install _build
```

Now you can use it in your CMake project by finding it again

```cmake
if(NOT "nlopt-f::nlopt-f")
  find_package("nlopt-f" REQUIRED)
endif()
# ...
target_link_libraries("${PROJECT_NAME}-lib" PRIVATE "nlopt-f::nlopt::wrap")
```


### Building with fpm

Invoke fpm in the project root with

```
fpm build
```

To run the testsuite use

```
fpm test
```

To use ``nlopt-f`` include it as dependency in your package manifest

```toml
[dependencies]
nlopt-f.git = "https://github.com/grimme-lab/nlopt-f"
```


## Usage

The Fortran bindings allow an object oriented usage of the NLopt C-API:

```f90
module example_funcs
  implicit none

  integer, parameter :: wp = kind(0.0d0)
  type :: constraint_data
    real(wp) :: d(2)
  end type

contains

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

end module example_funcs


program example
  use example_funcs
  use nlopt_wrap, only : nlopt_opt, nlopt_func, create, destroy
  use nlopt_enum, only : NLOPT_SUCCESS, algorithm_from_string
  implicit none
  type(nlopt_opt) :: opt
  real(wp) :: lb(2), x(2), minf
  integer :: stat
  type(constraint_data), target :: d1, d2
  real(wp), parameter :: xtol = 1.0e-4_wp

  call create(opt, algorithm_from_string("LD_MMA"), 2)

  call opt%get_lower_bounds(lb)

  lb(2) = 0.0_wp
  call opt%set_lower_bounds(lb)

  d1%d = [+2.0_wp, +0.0_wp]
  d2%d = [-1.0_wp, +1.0_wp]
  associate(&
      & f => nlopt_func(myfunc), &
      & fc1 => nlopt_func(myconstraint, d1), &
      & fc2 => nlopt_func(myconstraint, d2))
    call opt%set_min_objective(f)

    call opt%add_inequality_constraint(fc1, 1.0e-8_wp)
    call opt%add_inequality_constraint(fc2, 1.0e-8_wp)

    call opt%set_xtol_rel(xtol)

    x = [1.234_wp, 5.678_wp]
    call opt%optimize(x, minf, stat)
  end associate

  if (stat < NLOPT_SUCCESS) then
    write(*, '(a)') "NLopt failed!"
    stop 1
  endif

  write(*, '(a, *(1x, g0))') "Found minimum at", x
  write(*, '(a, *(1x, g0))') "Minimum value is", minf

  call destroy(opt)
end program example
```


## License

This project is free software: you can redistribute it and/or modify it under the terms of the [Apache License, Version 2.0](LICENSE-Apache) or [MIT license](LICENSE-MIT) at your opinion.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an _as is_ basis, without warranties or conditions of any kind, either express or implied. See the License for the specific language governing permissions and limitations under the License.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in this project by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or conditions.
