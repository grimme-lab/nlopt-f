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

testdrive_dep = dependency('test-drive', required: false)
if not testdrive_dep.found()
  testdrive_prj = subproject(
    'test-drive',
    required: not meson.is_subproject(),
    default_options: [
      'default_library=static',
    ],
  )
  if not testdrive_prj.found()
    subdir_done()
  endif
  testdrive_dep = testdrive_prj.get_variable('testdrive_dep')
endif

tests = [
  'enum',
  'opt',
  'version',
]

test_srcs = files(
  'main.f90',
)
foreach t : tests
  test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [nlopt_dep, testdrive_dep],
)

foreach t : tests
  test(t, tester, args: t)
endforeach
