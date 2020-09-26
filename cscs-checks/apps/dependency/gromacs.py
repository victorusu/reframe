# Copyright 2016-2020 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import contextlib
import itertools
import os
import fnmatch
from datetime import datetime

import reframe as rfm
import reframe.utility.sanity as sn
import reframe.utility as util
from reframe.core.backends import getlauncher

# @rfm.parameterized_test(['mpi+ownfftw+shared+gnu'],
#                         ['openmp+ownfftw+shared+gnu'],
#                         ['cuda+ownfftw+shared+gnu'],
#                         ['mpi+cuda+ownfftw+shared+gnu'],
#                         ['mpi+openmp+cuda+ownfftw+shared+gnu'],
#                         ['openmp+cuda+ownfftw+shared+gnu'])
@rfm.parameterized_test(['openmp+ownfftw+shared+gnu'])
class CompileGROMACSMasterTest(rfm.CompileOnlyRegressionTest):
    def __init__(self, variant):

        self.variant = variant

        self.valid_systems = ['daint:gpu','dom:gpu']
        self.modules = ['CMake'] # system's cmake it too old

        if 'gnu' in self.variant:
            self.valid_prog_environs += ['PrgEnv-gnu']
        elif 'intel' in self.variant:
            self.valid_prog_environs += ['PrgEnv-intel']
        elif 'cce' in self.variant:
            self.valid_prog_environs += ['PrgEnv-cray']
        elif 'pgi' in self.variant:
            self.valid_prog_environs += ['PrgEnv-pgi']

        mpi = 'ON' if 'mpi' in self.variant else 'OFF'
        openmp = 'ON' if 'openmp' in self.variant else 'OFF'

        if 'cuda' in self.variant:
            cuda = 'CUDA'
            cuda_opts = '-DCUDA_TOOLKIT_ROOT_DIR=$CUDATOOLKIT_HOME '
            self.modules += ['cudatoolkit']
        else:
            cuda = 'OFF'
            cuda_opts = ''

        if 'ownfftw' in self.variant:
            ownfftw = 'ON'
        else:
            ownfftw = 'OFF'
            self.modules += ['cray-fftw']

        if 'shared' in self.variant:
            self.variables = {'CRAYPE_LINK_TYPE': 'dynamic'}
            shared='-DBUILD_SHARED_LIBS=ON -DGMX_PREFER_STATIC_LIBS=OFF'
        else:
            self.variables = {'CRAYPE_LINK_TYPE': 'static'}
            shared='-DGMXAPI=OFF -DBUILD_SHARED_LIBS=OFF -DGMX_PREFER_STATIC_LIBS=ON'

        self.sourcesdir = 'https://gitlab.com/gromacs/gromacs.git'
        self.build_system = 'CMake'
        self.build_system.builddir = 'build'
        self.build_system.config_opts = [
            # '-DCMAKE_INSTALL_PREFIX=/apps/daint/UES/jenkins/7.0.UP02/gpu/easybuild/software/GROMACS/2020.3-CrayGNU-20.08-cuda ',
            '-DCMAKE_POSITION_INDEPENDENT_CODE=ON ',
            '-DCMAKE_VERBOSE_MAKEFILE=ON ',
            '-DCMAKE_BUILD_TYPE=Release ',
            f'{cuda_opts} ',
            f'-DGMX_BUILD_OWN_FFTW={ownfftw} ',
            f'-DGMX_OPENMP={openmp} ',
            f'-DGMX_GPU={cuda} ',
            '-DGMX_SIMD=AVX2_256 ',
            '-DGMX_CYCLE_SUBCOUNTERS=ON  ',
            f'-DGMX_MPI={mpi} ',
            f'{shared}'
        ]

        self.maintainers = ['VH']
        self.tags = {'ci'}

        self.sanity_patterns = sn.assert_found(r'100%', self.stdout)


class GromacsBaseCheck(rfm.RunOnlyRegressionTest):
    def __init__(self, output_file):

        # Reset sources dir relative to the SCS apps prefix
        self.sourcesdir = os.path.join(self.current_system.resourcesdir,
                                       'Gromacs', 'herflat')
        self.keep_files = [output_file]

        energy = sn.extractsingle(r'\s+Potential\s+Kinetic En\.\s+Total Energy'
                                  r'\s+Conserved En\.\s+Temperature\n'
                                  r'(\s+\S+){2}\s+(?P<energy>\S+)(\s+\S+){2}\n'
                                  r'\s+Pressure \(bar\)\s+Constr\. rmsd',
                                  output_file, 'energy', float, item=-1)
        energy_reference = -3270799.9

        self.sanity_patterns = sn.all([
            sn.assert_found('Finished mdrun', output_file),
            sn.assert_reference(energy, energy_reference, -0.001, 0.001)
        ])

        self.perf_patterns = {
            'perf': sn.extractsingle(r'Performance:\s+(?P<perf>\S+)',
                                     output_file, 'perf', float)
        }

        self.maintainers = ['VH']
        self.strict_check = False
        self.use_multithreading = False
        self.extra_resources = {
            'switches': {
                'num_switches': 1
            }
        }
        self.tags = {'scs', 'external-resources', 'ci'}


# @rfm.required_version('>=2.19')
# @rfm.parameterized_test(*([s, v, d]
#                           for s in ['small', 'large']
#                           for v in ['prod', 'maint']
#                           for d in ['mpi+openmp+cuda+ownfftw+shared+gnu']))
# class GromacsGPUCheck(GromacsBaseCheck):
#     def __init__(self, scale, variant, dependency):
#         super().__init__('md.log')
#         self.valid_systems = ['daint:gpu', 'tiger:gpu']
#         self.descr = 'GROMACS GPU check'
#         self.executable_opts = ['mdrun', '-dlb yes', '-ntomp 1', '-npme 0',
#                                 '-s herflat.tpr']
#         self.variables = {'CRAY_CUDA_MPS': '1'}
#         self.num_gpus_per_node = 1
#         if scale == 'small':
#             self.valid_systems += ['dom:gpu']
#             self.num_tasks = 72
#             self.num_tasks_per_node = 12
#         else:
#             self.num_tasks = 192
#             self.num_tasks_per_node = 12

#         self.depends_on("CompileGROMACSMasterTest_" + util.toalphanum(dependency))

#         references = {
#             'maint': {
#                 'large': {
#                     'daint:gpu': {'perf': (73.4, -0.10, None, 'ns/day')}
#                 }
#             },
#             'prod': {
#                 'small': {
#                     'dom:gpu': {'perf': (37.0, -0.05, None, 'ns/day')},
#                     'daint:gpu': {'perf': (35.0, -0.10, None, 'ns/day')}
#                 },
#                 'large': {
#                     'daint:gpu': {'perf': (63.0, -0.20, None, 'ns/day')}
#                 }
#             },
#         }
#         with contextlib.suppress(KeyError):
#             self.reference = references[variant][scale]

#         self.tags |= {'maintenance' if variant == 'maint' else 'production'}


@rfm.required_version('>=2.19')
@rfm.parameterized_test(['openmp+ownfftw+shared+gnu'])
class GromacsCPUCheck(GromacsBaseCheck):
    def __init__(self, dependency):
        super().__init__('md.log')
        # This is limitation of the dependency per partition
        self.valid_systems = ['daint:gpu','dom:gpu']
        # settings for the openmp case
        self.num_tasks = 1
        self.num_tasks_per_node = 12

        self.descr = 'GROMACS CPU check'
        self.executable_opts = ['mdrun', '-dlb yes',
                                f'-ntomp {self.num_tasks_per_node}', '-npme -1',
                                '-nb cpu', '-s herflat.tpr']

        if 'gnu' in dependency:
            self.valid_prog_environs += ['PrgEnv-gnu']
        elif 'intel' in dependency:
            self.valid_prog_environs += ['PrgEnv-intel']
        elif 'cce' in dependency:
            self.valid_prog_environs += ['PrgEnv-cray']
        elif 'pgi' in dependency:
            self.valid_prog_environs += ['PrgEnv-pgi']

        # this I love
        self.depends_on("CompileGROMACSMasterTest_" + util.toalphanum(dependency))


        self.reference = {
            'dom:gpu': {'perf': (40.0, -0.05, None, 'ns/day')},
            'daint:gpu': {'perf': (38.8, -0.10, None, 'ns/day')}
        }

    @rfm.require_deps
    # this I hate!
    def set_executable(self, CompileGROMACSMasterTest_openmp_ownfftw_shared_gnu):
        self.executable = os.path.join(
            CompileGROMACSMasterTest_openmp_ownfftw_shared_gnu().stagedir,
            'gmx'
        )
