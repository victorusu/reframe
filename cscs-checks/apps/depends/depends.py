# Copyright 2016-2020 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import itertools
import os
import fnmatch
from datetime import datetime

import reframe as rfm
import reframe.utility.sanity as sn
from reframe.core.backends import getlauncher

tests=['A', 'B', 'C', 'D']
# tests=['A', 'B']
# tests2=['C', 'D']
tests2=['E', 'F']

# level 0
@rfm.required_version('>=2.19')
@rfm.parameterized_test(*(tests))
class MyDepTest(rfm.RegressionTest):
    def __init__(self, name):
        # self.name = name
        self.valid_systems = ['daint:login']
        self.valid_prog_environs = ['builtin']
        self.executable = 'echo'
        self.executable_opts = [name]

# level 1
@rfm.required_version('>=2.19')
@rfm.parameterized_test(*([s, v]
                          for s in tests
                          for v in tests))
class MyDepTest(rfm.RegressionTest):
    def __init__(self, name, depends):
        # self.name = name
        self.valid_systems = ['daint:login']
        self.valid_prog_environs = ['builtin']
        self.executable = 'echo'
        self.executable_opts = [name]
        self.depends_on("MyDepTest_" + depends)

# level 1
@rfm.required_version('>=2.19')
@rfm.parameterized_test(*([s, v]
                          for s in tests2
                          for v in tests))
class MyDepTest(rfm.RegressionTest):
    def __init__(self, name, depends):
        # self.name = name
        self.valid_systems = ['daint:login']
        self.valid_prog_environs = ['builtin']
        self.executable = 'echo'
        self.executable_opts = [name]
        for d in tests:
            self.depends_on("MyDepTest_" + d)

# level 2
@rfm.required_version('>=2.19')
@rfm.parameterized_test(*([s, v, t]
                          for s in tests
                          for v in tests
                          for t in tests))
class MyDepTest(rfm.RegressionTest):
    def __init__(self, name, depends, depends2):
        # self.name = name
        self.valid_systems = ['daint:login']
        self.valid_prog_environs = ['builtin']
        self.executable = 'echo'
        self.executable_opts = [name]
        self.depends_on("MyDepTest_" + depends + "_" + depends2)

# level 2
@rfm.required_version('>=2.19')
@rfm.parameterized_test(*([s, v, t]
                          for s in tests2
                          for v in tests
                          for t in tests))
class MyDepTest(rfm.RegressionTest):
    def __init__(self, name, depends, depends2):
        # self.name = name
        self.valid_systems = ['daint:login']
        self.valid_prog_environs = ['builtin']
        self.executable = 'echo'
        self.executable_opts = [name]
        for d in tests:
            self.depends_on("MyDepTest_" + depends + "_" + d)


# level 3
@rfm.required_version('>=2.19')
@rfm.parameterized_test(*([s, v, t, w]
                          for s in tests2
                          for v in tests
                          for t in tests
                          for w in tests2))
class MyDepTest(rfm.RegressionTest):
    def __init__(self, name, depends, depends2, depends3):
        # self.name = name
        self.valid_systems = ['daint:login']
        self.valid_prog_environs = ['builtin']
        self.executable = 'echo'
        self.executable_opts = [name]
        for d in tests:
            self.depends_on("MyDepTest_" + depends + "_" + d + "_" + d)


# @rfm.required_version('>=2.19')
# @rfm.parameterized_test(*([s, v, t, w]
#                           for s in tests
#                           for v in tests
#                           for t in tests
#                           for w in ['A']))
# class MyDepTest(rfm.RegressionTest):
#     def __init__(self, name, depends, depends2, depends3):
#         # self.name = name
#         self.valid_systems = ['daint:login']
#         self.valid_prog_environs = ['builtin']
#         self.executable = 'echo'
#         self.executable_opts = [name]
#         self.depends_on("MyDepTest_" + depends + "_" + depends2 + "_" + depends3)
#         self.depends_on("MyDepTest_" + depends + "_" + depends2 + "_" + depends3)
