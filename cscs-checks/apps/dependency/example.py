# Copyright 2016-2020 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

# import reframe as rfm

# @rfm.simple_test
# class A(rfm.CompileOnlyRegressionTest):
#     def __init__(self):
#         self.valid_systems = ['daint:gpu']
#         self.valid_prog_environs = ['PrgEnv-gnu']

# @rfm.simple_test
# class B(rfm.RunOnlyRegressionTest):
#     def __init__(self):
#         self.valid_systems = ['daint:login']
#         self.valid_prog_environs = ['builtin', 'PrgEnv-gnu']

#         self.dep_name = "A"
#         self.depends_on(self.dep_name, how=rfm.DEPEND_EXACT,
#             subdeps={('daint:gpu', 'PrgEnv-gnu'): [('daint:login', self.current_environ.name)]})

import builtins
import os
import reframe as rfm
import reframe.utility.sanity as sn


def parent_part_env(part, env):
    def _parent_part_env(src, dst):
        if dst:
            return dst[0].split(':')[1] == part and dst[1] == env
        return False
    return _parent_part_env


@rfm.simple_test
class A(rfm.CompileOnlyRegressionTest):
    def __init__(self):
        self.valid_systems = ['daint:login', 'dom:login']
        self.valid_prog_environs = ['PrgEnv-gnu']

        self.sourcepath = 'hello_world'
        self.executable = './hello_world'
        self.build_system = 'SingleSource'

        self.sourcesdir = 'src/serial'
        self.sourcepath += '_serial.cpp'
        self.descr += ' Serial'

        self.sanity_patterns = sn.assert_not_found('error', self.stderr)

@rfm.parameterized_test(['builtin'], ['PrgEnv-gnu'])
class B(rfm.RunOnlyRegressionTest):
    def __init__(self, prgenv):
        self.valid_systems = ['daint:gpu', 'dom:gpu']
        self.valid_prog_environs = [prgenv]

        self.dep_name = 'A'
        # self.depends_on(self.dep_name, how=rfm.DEPEND_EXACT,
        #     subdeps={('gpu', prgenv): [('login', 'PrgEnv-gnu')]})
        self.depends_on(self.dep_name, when=parent_part_env('login','PrgEnv-gnu'))
        self.sanity_patterns = sn.assert_not_found('error', self.stderr)

    @rfm.run_after('setup')
    def set_executable(self):
        target = self.getdep(self.dep_name, 'PrgEnv-gnu')
        self.executable = os.path.join(target.stagedir, 'hello_world')
