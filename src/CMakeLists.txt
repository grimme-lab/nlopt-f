# This file is part of nlopt-wrap.
# SPDX-Identifier: Apache-2.0 OR MIT
#
# Licensed under either of Apache License, Version 2.0 or MIT license
# at your option; you may not use this file except in compliance with
# the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/nlopt_callback.f90"
  "${dir}/nlopt_enum.F90"
  "${dir}/nlopt_interface.F90"
  "${dir}/nlopt_wrap.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)
