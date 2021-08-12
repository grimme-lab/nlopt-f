# This file is part of nlopt-f.
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

os = host_machine.system()
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
endif

lib_deps += dependency('nlopt', version: '>=2.0')
if get_option('nlopt_version') == 'auto'
  nlopt_varr = lib_deps.get(-1).version().split('.')
else
  if get_option('nlopt_version').version_compare('>='+lib_deps.get(-1).version())
    error('Requested NLopt version is newer than available NLopt library')
  endif
  nlopt_varr = get_option('nlopt_version').split('.')
endif
nlopt_version = nlopt_varr.get(0, '0').to_int() * 10000 \
    + nlopt_varr.get(1, '0').to_int() * 100 \
    + nlopt_varr.get(2, '0').to_int()
lib_deps += declare_dependency(
  compile_args: '-DNLOPT_VERSION='+nlopt_version.to_string()
)
