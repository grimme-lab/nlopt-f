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

project(
  'nlopt-f',
  'fortran', 'c',
  version: '0.1.0',
  license: 'Apache-2.0 OR MIT',
  meson_version: '>=0.55',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
lib_deps = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
nlopt_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: include_directories('include'),
  install: install,
)

# Export dependency for other projects and test suite
nlopt_inc = nlopt_lib.private_dir_include()
nlopt_dep = declare_dependency(
  link_with: nlopt_lib,
  include_directories: nlopt_inc,
  dependencies: lib_deps,
)

# Package the license files
nlopt_lic = files(
  'LICENSE-Apache',
  'LICENSE-MIT',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    nlopt_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    'config'/'install-mod.py',
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    nlopt_lib,
    description: 'Light-weight tight-binding framework',
    subdirs: ['', module_id],
  )
endif

# add the testsuite
subdir('test')
