# Copyright 2016-2020 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import reframe as rfm

@rfm.simple_test
class A(rfm.CompileOnlyRegressionTest):
    def __init__(self):
        self.valid_systems = ['daint:login']
        self.valid_prog_environs = ['PrgEnv-gnu']

@rfm.simple_test
class B(rfm.RunOnlyRegressionTest):
    def __init__(self):
        self.valid_systems = ['daint:login']
        self.valid_prog_environs = ['builtin']

        self.dep_name = "A"
        self.depends_on(self.dep_name, rfm.DEPEND_EXACT, )
