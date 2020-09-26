# Copyright 2016-2020 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

from datetime import datetime

import re
import reframe as rfm
import reframe.utility as util
import reframe.utility.sanity as sn
import os

class CompileHelloWorldBaseTest(rfm.CompileOnlyRegressionTest):
    def __init__(self, variant, lang, linkage):
        self.linkage = linkage
        self.variables = {'CRAYPE_LINK_TYPE': linkage}
        self.prgenv_flags = {}
        self.lang_names = {
            'c': 'C',
            'cpp': 'C++',
            'f90': 'Fortran 90'
        }
        self.descr = f'{self.lang_names[lang]} Hello World'
        self.sourcepath = 'hello_world'
        self.executable = './hello_world'
        self.build_system = 'SingleSource'
        # self.valid_systems = ['daint:gpu']
        self.valid_systems = ['daint:gpu', 'daint:mc', 'dom:gpu', 'dom:mc',
                              'kesch:cn', 'tiger:gpu', 'arolla:cn',
                              'arolla:pn', 'tsa:cn', 'tsa:pn']

        self.valid_prog_environs = ['PrgEnv-cray',
                                    'PrgEnv-intel', 'PrgEnv-gnu', 'PrgEnv-pgi',
                                    'PrgEnv-gnu-nocuda', 'PrgEnv-pgi-nocuda']

        if self.current_system.name in ['kesch', 'arolla', 'tsa']:
            self.exclusive_access = True

        # Removing static compilation from kesch
        if (self.current_system.name in ['kesch'] and linkage == 'static'):
            self.valid_prog_environs = []

        self.compilation_time_seconds = None

        result = sn.findall(r'Hello World from thread \s*(\d+) out '
                            r'of \s*(\d+) from process \s*(\d+) out of '
                            r'\s*(\d+)', self.stdout)

        self.sanity_patterns = sn.assert_not_found('error', self.stderr)

        self.perf_patterns = {
            'compilation_time': sn.getattr(self, 'compilation_time_seconds')
        }
        self.reference = {
            '*': {
                'compilation_time': (60, None, 0.1, 's')
            }
        }

        self.maintainers = ['VH', 'EK']
        self.tags = {'production', 'craype'}

    @rfm.run_before('compile')
    def setflags(self):
        envname = re.sub(r'(PrgEnv-\w+).*', lambda m: m.group(1),
                         self.current_environ.name)
        prgenv_flags = self.prgenv_flags[envname]
        self.build_system.cflags = prgenv_flags
        self.build_system.cxxflags = prgenv_flags
        self.build_system.fflags = prgenv_flags

    @rfm.run_before('compile')
    def compile_timer_start(self):
        self.compilation_time_seconds = datetime.now()

    @rfm.run_after('compile')
    def compile_timer_end(self):
        elapsed = datetime.now() - self.compilation_time_seconds
        self.compilation_time_seconds = elapsed.total_seconds()


@rfm.required_version('>=2.14')
@rfm.parameterized_test(*([lang, linkage]
                          for lang in ['cpp', 'c', 'f90']
                          for linkage in ['dynamic', 'static']))
class CompileHelloWorldTestSerial(CompileHelloWorldBaseTest):
    def __init__(self, lang, linkage):
        super().__init__('serial', lang, linkage)
        self.valid_systems += ['kesch:pn', 'arolla:pn', 'tsa:pn']
        self.valid_prog_environs += ['PrgEnv-gnu-nompi', 'PrgEnv-pgi-nompi',
                                     'PrgEnv-gnu-nompi-nocuda',
                                     'PrgEnf-pgi-nompi-nocuda']
        self.sourcesdir = 'src/serial'
        self.sourcepath += '_serial.' + lang
        self.descr += ' Serial ' + linkage.capitalize()
        self.prgenv_flags = {
            'PrgEnv-cray': [],
            'PrgEnv-cray_classic': [],
            'PrgEnv-gnu': [],
            'PrgEnv-intel': [],
            'PrgEnv-pgi': []
        }
        self.num_tasks = 1
        self.num_tasks_per_node = 1
        self.num_cpus_per_task = 1
        if self.current_system.name == 'kesch' and linkage == 'dynamic':
            self.valid_prog_environs += ['PrgEnv-cray-nompi',
                                         'PrgEnv-pgi-nompi',
                                         'PrgEnv-gnu-nompi']
        elif (self.current_system.name in ['arolla', 'tsa'] and
              linkage == 'dynamic'):
            self.valid_prog_environs += ['PrgEnv-pgi-nompi',
                                         'PrgEnv-pgi-nompi-nocuda',
                                         'PrgEnv-gnu-nompi',
                                         'PrgEnv-gnu-nompi-nocuda']


@rfm.required_version('>=2.14')
@rfm.parameterized_test(*([lang, linkage]
                          for lang in ['cpp', 'c', 'f90']
                          for linkage in ['dynamic', 'static']))
class CompileHelloWorldTestOpenMP(CompileHelloWorldBaseTest):
    def __init__(self, lang, linkage):
        super().__init__('openmp', lang, linkage)
        self.valid_systems += ['kesch:pn', 'arolla:pn', 'tsa:pn']
        self.sourcesdir = 'src/openmp'
        self.sourcepath += '_openmp.' + lang
        self.descr += ' OpenMP ' + str.capitalize(linkage)
        self.prgenv_flags = {
            'PrgEnv-cray': ['-homp' if lang == 'F90' else '-fopenmp'],
            'PrgEnv-cray_classic': ['-homp'],
            'PrgEnv-gnu': ['-fopenmp'],
            'PrgEnv-intel': ['-qopenmp'],
            'PrgEnv-pgi': ['-mp']
        }
        if self.current_system.name == 'kesch':
            self.prgenv_flags['PrgEnv-cray'] = ['-homp']

        self.num_tasks = 1
        self.num_tasks_per_node = 1
        self.num_cpus_per_task = 4
        if self.current_system.name == 'kesch' and linkage == 'dynamic':
            self.valid_prog_environs += ['PrgEnv-cray-nompi',
                                         'PrgEnv-pgi-nompi',
                                         'PrgEnv-gnu-nompi']
        elif (self.current_system.name in ['arolla', 'tsa'] and
              linkage == 'dynamic'):
            self.valid_prog_environs += ['PrgEnv-pgi-nompi',
                                         'PrgEnv-pgi-nompi-nocuda',
                                         'PrgEnv-gnu-nompi',
                                         'PrgEnv-gnu-nompi-nocuda']

        # On SLURM there is no need to set OMP_NUM_THREADS if one defines
        # num_cpus_per_task, but adding for completeness and portability
        self.variables = {
            'OMP_NUM_THREADS': str(self.num_cpus_per_task)
        }


class RunHelloWorldBaseTest(rfm.RunOnlyRegressionTest):
    def __init__(self, variant, lang, linkage):
        self.linkage = linkage
        self.lang = lang
        self.variables = {'CRAYPE_LINK_TYPE': linkage}
        self.prgenv_flags = {}
        self.lang_names = {
            'c': 'C',
            'cpp': 'C++',
            'f90': 'Fortran 90'
        }

        self.descr = f'{self.lang_names[lang]} Hello World'
        # self.valid_systems = ['daint:gpu']
        self.valid_systems = ['daint:gpu', 'daint:mc', 'dom:gpu', 'dom:mc',
                              'kesch:cn', 'tiger:gpu', 'arolla:cn',
                              'arolla:pn', 'tsa:cn', 'tsa:pn']

        self.valid_prog_environs = ['PrgEnv-cray',
                                    'PrgEnv-intel', 'PrgEnv-gnu', 'PrgEnv-pgi',
                                    'PrgEnv-gnu-nocuda', 'PrgEnv-pgi-nocuda']

        if self.current_system.name in ['kesch', 'arolla', 'tsa']:
            self.exclusive_access = True

        # Removing static compilation from kesch
        if (self.current_system.name in ['kesch'] and linkage == 'static'):
            self.valid_prog_environs = []

        result = sn.findall(r'Hello World from thread \s*(\d+) out '
                            r'of \s*(\d+) from process \s*(\d+) out of '
                            r'\s*(\d+)', self.stdout)

        num_tasks = sn.getattr(self, 'num_tasks')
        num_cpus_per_task = sn.getattr(self, 'num_cpus_per_task')

        def tid(match):
            return int(match.group(1))

        def num_threads(match):
            return int(match.group(2))

        def rank(match):
            return int(match.group(3))

        def num_ranks(match):
            return int(match.group(4))

        self.sanity_patterns = sn.all(
            sn.chain(
                [sn.assert_eq(sn.count(result), num_tasks*num_cpus_per_task)],
                sn.map(lambda x: sn.assert_lt(tid(x), num_threads(x)), result),
                sn.map(lambda x: sn.assert_lt(rank(x), num_ranks(x)), result),
                sn.map(
                    lambda x: sn.assert_lt(tid(x), num_cpus_per_task), result
                ),
                sn.map(
                    lambda x: sn.assert_eq(num_threads(x), num_cpus_per_task),
                    result
                ),
                sn.map(lambda x: sn.assert_lt(rank(x), num_tasks), result),
                sn.map(
                    lambda x: sn.assert_eq(num_ranks(x), num_tasks), result
                ),
            )
        )

        self.maintainers = ['VH', 'EK']
        self.tags = {'production', 'craype'}


# @rfm.required_version('>=2.14')
# @rfm.parameterized_test(*([lang, linkage]
#                           for lang in ['cpp', 'c', 'f90']
#                           for linkage in ['dynamic', 'static']))
class RunHelloWorldTestSerial(RunHelloWorldBaseTest):
    def __init__(self, lang, linkage):
        super().__init__('serial', lang, linkage)
        self.valid_systems += ['kesch:pn', 'arolla:pn', 'tsa:pn']
        self.valid_prog_environs += ['PrgEnv-gnu-nompi', 'PrgEnv-pgi-nompi',
                                     'PrgEnv-gnu-nompi-nocuda',
                                     'PrgEnf-pgi-nompi-nocuda']
        self.sourcesdir = 'src/serial'
        self.sourcepath += '_serial.' + lang
        self.descr += ' Serial ' + linkage.capitalize()
        self.prgenv_flags = {
            'PrgEnv-cray': [],
            'PrgEnv-cray_classic': [],
            'PrgEnv-gnu': [],
            'PrgEnv-intel': [],
            'PrgEnv-pgi': []
        }
        self.num_tasks = 1
        self.num_tasks_per_node = 1
        self.num_cpus_per_task = 1
        if self.current_system.name == 'kesch' and linkage == 'dynamic':
            self.valid_prog_environs += ['PrgEnv-cray-nompi',
                                         'PrgEnv-pgi-nompi',
                                         'PrgEnv-gnu-nompi']
        elif (self.current_system.name in ['arolla', 'tsa'] and
              linkage == 'dynamic'):
            self.valid_prog_environs += ['PrgEnv-pgi-nompi',
                                         'PrgEnv-pgi-nompi-nocuda',
                                         'PrgEnv-gnu-nompi',
                                         'PrgEnv-gnu-nompi-nocuda']

        self.dep_name = f'CompileHelloWorldTestSerial_{lang}_{linkage}'
        self.depends_on(self.dep_name)


@rfm.simple_test
class RunHelloWorldTestSerial_cpp_dynamic(RunHelloWorldTestSerial):
    def __init__(self):
        super().__init__('cpp', 'dynamic')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestSerial_cpp_dynamic):
        self.executable = os.path.join(
            CompileHelloWorldTestSerial_cpp_dynamic().stagedir,
            'hello_world'
        )


@rfm.simple_test
class RunHelloWorldTestSerial_c_dynamic(RunHelloWorldTestSerial):
    def __init__(self):
        super().__init__('c', 'dynamic')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestSerial_c_dynamic):
        self.executable = os.path.join(
            CompileHelloWorldTestSerial_c_dynamic().stagedir,
            'hello_world'
        )


@rfm.simple_test
class RunHelloWorldTestSerial_f90_dynamic(RunHelloWorldTestSerial):
    def __init__(self):
        super().__init__('f90', 'dynamic')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestSerial_f90_dynamic):
        self.executable = os.path.join(
            CompileHelloWorldTestSerial_f90_dynamic().stagedir,
            'hello_world'
        )


@rfm.simple_test
class RunHelloWorldTestSerial_cpp_static(RunHelloWorldTestSerial):
    def __init__(self):
        super().__init__('cpp', 'static')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestSerial_cpp_static):
        self.executable = os.path.join(
            CompileHelloWorldTestSerial_cpp_static().stagedir,
            'hello_world'
        )


@rfm.simple_test
class RunHelloWorldTestSerial_c_static(RunHelloWorldTestSerial):
    def __init__(self):
        super().__init__('c', 'static')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestSerial_c_static):
        self.executable = os.path.join(
            CompileHelloWorldTestSerial_c_static().stagedir,
            'hello_world'
        )


@rfm.simple_test
class RunHelloWorldTestSerial_f90_static(RunHelloWorldTestSerial):
    def __init__(self):
        super().__init__('f90', 'static')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestSerial_f90_static):
        self.executable = os.path.join(
            CompileHelloWorldTestSerial_f90_static().stagedir,
            'hello_world'
        )


# @rfm.required_version('>=2.14')
# @rfm.parameterized_test(*([lang, linkage]
#                           for lang in ['cpp', 'c', 'f90']
#                           for linkage in ['dynamic', 'static']))
class RunHelloWorldTestOpenMP(RunHelloWorldBaseTest):
    def __init__(self, lang, linkage):
        super().__init__('openmp', lang, linkage)
        self.valid_systems += ['kesch:pn', 'arolla:pn', 'tsa:pn']
        self.sourcesdir = 'src/openmp'
        self.sourcepath += '_openmp.' + lang
        self.descr += ' OpenMP ' + str.capitalize(linkage)
        self.prgenv_flags = {
            'PrgEnv-cray': ['-homp' if lang == 'F90' else '-fopenmp'],
            'PrgEnv-cray_classic': ['-homp'],
            'PrgEnv-gnu': ['-fopenmp'],
            'PrgEnv-intel': ['-qopenmp'],
            'PrgEnv-pgi': ['-mp']
        }
        if self.current_system.name == 'kesch':
            self.prgenv_flags['PrgEnv-cray'] = ['-homp']

        self.num_tasks = 1
        self.num_tasks_per_node = 1
        self.num_cpus_per_task = 4
        if self.current_system.name == 'kesch' and linkage == 'dynamic':
            self.valid_prog_environs += ['PrgEnv-cray-nompi',
                                         'PrgEnv-pgi-nompi',
                                         'PrgEnv-gnu-nompi']
        elif (self.current_system.name in ['arolla', 'tsa'] and
              linkage == 'dynamic'):
            self.valid_prog_environs += ['PrgEnv-pgi-nompi',
                                         'PrgEnv-pgi-nompi-nocuda',
                                         'PrgEnv-gnu-nompi',
                                         'PrgEnv-gnu-nompi-nocuda']

        # On SLURM there is no need to set OMP_NUM_THREADS if one defines
        # num_cpus_per_task, but adding for completeness and portability
        self.variables = {
            'OMP_NUM_THREADS': str(self.num_cpus_per_task)
        }

        self.dep_name = f'CompileHelloWorldTestOpenMP_{lang}_{linkage}'
        self.depends_on(self.dep_name)

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestOpenMP_cpp_dynamic):
        if self.dep_name == f'CompileHelloWorldTestOpenMP_{self.lang}_{self.linkage}':
            self.executable = os.path.join(
                CompileHelloWorldTestOpenMP_cpp_dynamic().stagedir,
                'hello_world'
            )

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestOpenMP_c_dynamic):
        if self.dep_name == f'CompileHelloWorldTestOpenMP_{self.lang}_{self.linkage}':
            self.executable = os.path.join(
                CompileHelloWorldTestOpenMP_c_dynamic().stagedir,
                'hello_world'
            )

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestOpenMP_f90_dynamic):
        if self.dep_name == f'CompileHelloWorldTestOpenMP_{self.lang}_{self.linkage}':
            self.executable = os.path.join(
                CompileHelloWorldTestOpenMP_f90_dynamic().stagedir,
                'hello_world'
            )

@rfm.simple_test
class RunHelloWorldTestOpenMP_cpp_dynamic(RunHelloWorldTestOpenMP):
    def __init__(self):
        super().__init__('cpp', 'dynamic')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestOpenMP_cpp_dynamic):
        self.executable = os.path.join(
            CompileHelloWorldTestOpenMP_cpp_dynamic().stagedir,
            'hello_world'
        )

@rfm.simple_test
class RunHelloWorldTestOpenMP_c_dynamic(RunHelloWorldTestOpenMP):
    def __init__(self):
        super().__init__('c', 'dynamic')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestOpenMP_c_dynamic):
        self.executable = os.path.join(
            CompileHelloWorldTestOpenMP_c_dynamic().stagedir,
            'hello_world'
        )

@rfm.simple_test
class RunHelloWorldTestOpenMP_f90_dynamic(RunHelloWorldTestOpenMP):
    def __init__(self):
        super().__init__('f90', 'dynamic')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestOpenMP_f90_dynamic):
        self.executable = os.path.join(
            CompileHelloWorldTestOpenMP_f90_dynamic().stagedir,
            'hello_world'
        )

@rfm.simple_test
class RunHelloWorldTestOpenMP_cpp_static(RunHelloWorldTestOpenMP):
    def __init__(self):
        super().__init__('cpp', 'static')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestOpenMP_cpp_static):
        self.executable = os.path.join(
            CompileHelloWorldTestOpenMP_cpp_static().stagedir,
            'hello_world'
        )

@rfm.simple_test
class RunHelloWorldTestOpenMP_c_static(RunHelloWorldTestOpenMP):
    def __init__(self):
        super().__init__('c', 'static')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestOpenMP_c_static):
        self.executable = os.path.join(
            CompileHelloWorldTestOpenMP_c_static().stagedir,
            'hello_world'
        )

@rfm.simple_test
class RunHelloWorldTestOpenMP_f90_static(RunHelloWorldTestOpenMP):
    def __init__(self):
        super().__init__('f90', 'static')

    @rfm.require_deps
    def set_executable(self, CompileHelloWorldTestOpenMP_f90_static):
        self.executable = os.path.join(
            CompileHelloWorldTestOpenMP_f90_static().stagedir,
            'hello_world'
        )
