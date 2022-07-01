# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.

set(_lib "test-drive")
set(_pkg "TEST_DRIVE")
set(_url "https://github.com/fortran-lang/test-drive")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  if(DEFINED "${PROJECT_NAME}-dependency-method")
    set("${_pkg}_FIND_METHOD" "${${PROJECT_NAME}-dependency-method}")
  else()
    set("${_pkg}_FIND_METHOD" "cmake" "pkgconf" "subproject" "fetch")
  endif()
  set("_${_pkg}_FIND_METHOD")
endif()

foreach(method in ITEMS ${${_pkg}_FIND_METHOD})
  if(TARGET "${_lib}::${_lib}")
    break()
  endif()

  if("${method}" STREQUAL "cmake")
    message(STATUS "${_lib}: Find installed package")
    find_package("${_lib}" CONFIG)
    if("${_lib}_FOUND")
      message(STATUS "${_lib}: Found installed package")
      break()
    endif()
  endif()

  if("${method}" STREQUAL "pkgconf")
    find_package(PkgConfig QUIET)
    pkg_check_modules("${_pkg}" QUIET "${_lib}")
    if("${_pkg}_FOUND")
      message(STATUS "Found ${_lib} via pkg-config")

      add_library("${_lib}::${_lib}" INTERFACE IMPORTED)
      target_link_libraries(
        "${_lib}::${_lib}"
        INTERFACE
        "${${_pkg}_LINK_LIBRARIES}"
        )
      target_include_directories(
        "${_lib}::${_lib}"
        INTERFACE
        "${${_pkg}_INCLUDE_DIRS}"
        )
      break()
    endif()
  endif()

  if("${method}" STREQUAL "subproject")
    set("${_pkg}_SOURCE_DIR" "${PROJECT_SOURCE_DIR}/subprojects/${_lib}")
    set("${_pkg}_BINARY_DIR" "${PROJECT_BINARY_DIR}/subprojects/${_lib}")
    if(EXISTS "${${_pkg}_SOURCE_DIR}/CMakeLists.txt")
      message(STATUS "Include ${_lib} from subprojects")
      add_subdirectory(
        "${${_pkg}_SOURCE_DIR}"
        "${${_pkg}_BINARY_DIR}"
        )

      add_library("${_lib}::${_lib}" INTERFACE IMPORTED)
      target_link_libraries("${_lib}::${_lib}" INTERFACE "${_lib}")

      # We need the module directory in the subproject before we finish the configure stage
      if(NOT EXISTS "${${_pkg}_BINARY_DIR}/include")
        file(MAKE_DIRECTORY "${${_pkg}_BINARY_DIR}/include")
      endif()

      break()
    endif()
  endif()

  if("${method}" STREQUAL "fetch")
    message(STATUS "Retrieving ${_lib} from ${_url}")
    include(FetchContent)
    FetchContent_Declare(
      "${_lib}"
      GIT_REPOSITORY "${_url}"
      GIT_TAG "HEAD"
      )
    FetchContent_MakeAvailable("${_lib}")

    add_library("${_lib}::${_lib}" INTERFACE IMPORTED)
    target_link_libraries("${_lib}::${_lib}" INTERFACE "${_lib}")

    # We need the module directory in the subproject before we finish the configure stage
    FetchContent_GetProperties("${_lib}" BINARY_DIR "${_pkg}_BINARY_DIR")
    if(NOT EXISTS "${${_pkg}_BINARY_DIR}/include")
      file(MAKE_DIRECTORY "${${_pkg}_BINARY_DIR}/include")
    endif()

    break()
  endif()

endforeach()

if(NOT TARGET "${_lib}::${_lib}")
  message(FATAL_ERROR "Could not find dependency ${_lib}")
endif()

if(DEFINED "_${_pkg}_FIND_METHOD")
  unset("${_pkg}_FIND_METHOD")
  unset("_${_pkg}_FIND_METHOD")
endif()
unset(_lib)
unset(_pkg)
unset(_url)
