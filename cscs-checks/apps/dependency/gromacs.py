# Copyright 2016-2020 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import contextlib
import itertools
import os
import re
import fnmatch
from datetime import datetime

import reframe as rfm
import reframe.utility.sanity as sn
import reframe.utility as util
from reframe.core.backends import getlauncher

# GROMACS_INSTALLATION_PREFIX = '/scratch/snx3000tds/hvictor/reframe-spack-tests'
GROMACS_INSTALLATION_PREFIX = '/scratch/snx3000/hvictor/reframe-run/gromacs'

PARAMETERISED_TESTS = [
    # ['mpi+openmp+cuda+ownfftw+shared+gnu', '2020.1', 'prod'],
    # ['mpi+openmp+cuda+ownfftw+shared+gnu', '2020.2', 'prod'],
    # ['mpi+openmp+cuda+ownfftw+shared+gnu', '2020.3', 'prod'],
    ['mpi+openmp+cuda+ownfftw+shared+gnu', '2020.4', 'prod'],
    # ['mpi+openmp+cuda+ownfftw+shared+gnu', '2019.6', 'prod'],
    # ['mpi+openmp+cuda+ownfftw+shared+gnu', '2018.8', 'prod'],

    # ['mpi+openmp+cuda+ownfftw+gnu', '2020.4', 'prod'],
    # ['mpi+openmp+cuda+ownfftw+gnu', '2019.6', 'prod'],
    # ['mpi+openmp+cuda+ownfftw+gnu', '2018.8', 'prod'],

    # ['mpi+openmp+cuda+ownfftw+shared+intel', '2020.4', 'prod'],
    # ['mpi+openmp+cuda+ownfftw+shared+intel', '2019.6', 'prod'],
    # ['mpi+openmp+cuda+ownfftw+shared+intel', '2018.8', 'prod'],

    # ['mpi+openmp+cuda+ownfftw+intel', '2020.4', 'prod'],
    # ['mpi+openmp+cuda+ownfftw+intel', '2019.6', 'prod'],
    # ['mpi+openmp+cuda+ownfftw+intel', '2018.8', 'prod'],
]


def parent_part_env(part, env):
    def _parent_part_env(src, dst):
        if dst:
            return dst[0].split(':')[1] == part and dst[1] == env
        return False
    return _parent_part_env


def download_gromacs(version):
    cmds = []
    cmds += [f'wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-{version}.tar.gz']
    cmds += [f'tar xf gromacs-{version}.tar.gz']
    cmds += [f'cd gromacs-{version}']
    return cmds


def get_prgenv(variant):
    prgenv = 'builtin'
    if 'gnu' in variant:
        prgenv= 'PrgEnv-gnu'
    elif 'intel' in variant:
        prgenv = 'PrgEnv-intel'
    elif 'cce' in variant:
        prgenv = 'PrgEnv-cray'
    elif 'pgi' in variant:
        prgenv = 'PrgEnv-pgi'

    return prgenv


def get_installation_path(check):
    install_path = os.path.join(check.variant,
                                f'gromacs-{check.gromacs_version}',
                                f'{util.toalphanum(check.variant)}')
    if check.install_type in ['stage']:
        install_path = check.stagedir
    elif check.install_type in ['output']:
        install_path = check.outputdir
    elif check.install_type in ['prod']:
        install_path = os.path.join(GROMACS_INSTALLATION_PREFIX,
                                    f'gromacs-{check.gromacs_version}',
                                    f'{util.toalphanum(check.variant)}')

    return install_path


def clone_or_download_gromacs(gromacs_version):
    if gromacs_version in ['develop', 'latest'] or len(gromacs_version) > 7:
        cmds = [
            f'git clone https://gitlab.com/gromacs/gromacs.git gromacs-{gromacs_version}',
            f'cd gromacs-{gromacs_version}'
        ]
        if len(gromacs_version) > 7:
            cmds += [
                f'git checkout -b {gromacs_version} {gromacs_version}'
            ]
    else:
        cmds = [
            f'wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-{version}.tar.gz',
            f'tar xf gromacs-{gromacs_version}.tar.gz',
            f'cd gromacs-{gromacs_version}'
        ]
    return cmds


@rfm.parameterized_test(*PARAMETERISED_TESTS)
class CompileGROMACSTest(rfm.CompileOnlyRegressionTest):
    def __init__(self, variant, gromacs_version, installation_type):
        self.variant = variant
        self.gromacs_version = gromacs_version
        self.install_type = installation_type

        self.valid_systems = ['daint:login', 'dom:login']
        self.modules = ['daint-gpu', 'CMake'] # system's cmake it too old
        self.valid_prog_environs = [get_prgenv(variant)]

        mpi = 'ON' if 'mpi' in self.variant else 'OFF'
        openmp = 'ON' if 'openmp' in self.variant else 'OFF'

        if 'cuda' in self.variant:
            cuda = 'CUDA'
            self.modules += ['cudatoolkit']
            cuda_opts = '-DCUDA_TOOLKIT_ROOT_DIR=$CUDATOOLKIT_HOME '
            # self.variables['CUDA_HOME'] = '$CUDATOOLKIT_HOME'
        else:
            cuda = 'OFF'
            cuda_opts = ''

        if 'ownfftw' in self.variant:
            ownfftw = 'ON'
        else:
            ownfftw = 'OFF'
            self.modules += ['cray-fftw']

        if 'shared' in self.variant:
            self.variables['CRAYPE_LINK_TYPE'] = 'dynamic'
            shared='-DBUILD_SHARED_LIBS=ON -DGMX_PREFER_STATIC_LIBS=OFF'
        else:
            self.variables['CRAYPE_LINK_TYPE'] = 'static'
            shared='-DGMXAPI=OFF -DBUILD_SHARED_LIBS=OFF -DGMX_PREFER_STATIC_LIBS=ON'

        if gromacs_version == 'master':
            self.sourcesdir = 'https://gitlab.com/gromacs/gromacs.git'
        else:
            self.prebuild_cmds = download_gromacs(gromacs_version)

        self.postbuild_cmds = ['make install']

        self.build_system = 'CMake'
        self.build_system.sourcesdir = None
        self.build_system.builddir = 'build'
        self.build_system.max_concurrency = 8
        self.build_system.config_opts = [
            '-DCMAKE_POSITION_INDEPENDENT_CODE=ON',
            '-DCMAKE_VERBOSE_MAKEFILE=ON',
            '-DCMAKE_BUILD_TYPE=Release',
            f'{cuda_opts}',
            f'-DGMX_BUILD_OWN_FFTW={ownfftw}',
            f'-DGMX_OPENMP={openmp}',
            f'-DGMX_GPU={cuda}',
            '-DGMX_CYCLE_SUBCOUNTERS=ON ',
            f'-DGMX_MPI={mpi}',
            f'{shared}'
        ]
        if self.current_system.name in ['daint', 'dom']:
            self.build_system.config_opts += ['-DGMX_SIMD=AVX2_256']

        self.sanity_patterns = sn.all([
            sn.assert_found(r'.*gromacs-{0}.tar\.gz.*saved'.format(gromacs_version), self.stderr),
            sn.assert_found(r'Enabling RDTSCP support', self.stdout),
            sn.assert_found(r'A library with LAPACK API found', self.stdout),
            sn.assert_found(r'Build files have been written to', self.stdout),
            sn.assert_found(r'\[100%\] Built target gmx', self.stdout),
            sn.assert_found(r'Install configuration: "Release"', self.stdout),
            # sn.assert_found(r'\[100%\] Built target template', self.stdout),
            # sn.assert_found(r'\[100%\] Built target gmxapi', self.stdout),
            sn.assert_found(r'Install the project', self.stdout)
            # sn.assert_lt(energy_diff, 14.9)
        ])

        self.maintainers = ['VH']
        self.tags = {'ci', 'ci-build', 'external-resources'}


    @rfm.run_after('setup')
    def config_gromacs(self):
        install_path = get_installation_path(self)

        self.build_system.config_opts += [
            f'-DCMAKE_INSTALL_PREFIX={install_path}',
        ]

        self.variables['GROMACS_ROOT'] = f'{install_path}'
        self.variables['PATH'] = r'$GROMACS_ROOT/bin:$PATH'


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
        self.tags = {'ci', 'ci-run', 'external-resources'}


# TODO: PARAMETERISE BASED ON GROMACS CLI FLAGS and input files
@rfm.parameterized_test(*PARAMETERISED_TESTS)
class GromacsCheck(GromacsBaseCheck):
    def __init__(self, variant, gromacs_version, installation_type):
        super().__init__('md.log')
        self.variant = variant

        self.valid_systems = ['daint:gpu','dom:gpu']
        self.valid_prog_environs = ['builtin']

        # generic single node job
        self.num_tasks = 1
        self.num_tasks_per_node = 1
        self.num_cpus_per_task = os.cpu_count()

        if 'cuda' in self.variant:
            self.variables['CRAY_CUDA_MPS'] = '1'
            self.num_gpus_per_node = 1

        self.descr = 'GROMACS RunTime check'

        self.reference = {
            'dom:gpu': {'perf': (40.0, -0.05, None, 'ns/day')},
            'daint:gpu': {'perf': (38.8, -0.10, None, 'ns/day')}
        }

        self.dep_name = re.sub(r'GromacsCheck', r'CompileGROMACSTest', self.name)
        self.depends_on(self.dep_name, when=parent_part_env('login', get_prgenv(self.variant)))

    @rfm.run_after('setup')
    def set_num_tasks(self):
        if self.current_partition.fullname in ['daint:gpu', 'dom:gpu']:
            self.num_tasks = 72
            self.num_tasks_per_node = 12
            self.num_cpus_per_task = 1

    @rfm.run_after('setup')
    def setup_gromacs(self):
        target = self.getdep(self.dep_name, get_prgenv(self.variant))

        self.variables.update(target.variables)
        self.modules += target.modules

        self.executable = 'gmx_mpi' if 'mpi' in self.variant else 'gmx'
        nb_type = 'gpu' if 'cuda' in self.variant else 'cpu'
        self.executable_opts = ['mdrun', '-dlb yes',
                                f'-ntomp {self.num_cpus_per_task}', '-npme -1',
                                f'-nb {nb_type}', '-s herflat.tpr']


# TODO: PARAMETERISE BASED ON GROMACS CLI FLAGS and input files
@rfm.parameterized_test(*PARAMETERISED_TESTS)
class GromacsIndependentCheck(GromacsBaseCheck):
    def __init__(self, variant, gromacs_version, installation_type):
        super().__init__('md.log')
        self.variant = variant
        self.gromacs_version = gromacs_version
        self.install_type = installation_type

        self.valid_systems = ['daint:gpu','dom:gpu']
        self.valid_prog_environs = [get_prgenv(self.variant)]

        # generic single node job
        self.num_tasks = 1
        self.num_tasks_per_node = 1
        self.num_cpus_per_task = os.cpu_count()

        self.descr = 'GROMACS RunTime check'

        self.reference = {
            'dom:gpu': {'perf': (40.0, -0.05, None, 'ns/day')},
            'daint:gpu': {'perf': (38.8, -0.10, None, 'ns/day')}
        }

    @rfm.run_after('setup')
    def set_num_tasks(self):
        if self.current_partition.fullname in ['daint:gpu', 'dom:gpu']:
            self.num_tasks = 72
            self.num_tasks_per_node = 12
            self.num_cpus_per_task = 1
            if 'cuda' in self.variant:
                self.variables['CRAY_CUDA_MPS'] = '1'
                self.num_gpus_per_node = 1
                self.modules = ['cudatoolkit']

    # This has to be defined after set_num_tasks because it depends on self.num_cpus_per_task
    @rfm.run_after('setup')
    def config_gromacs(self):
        install_path = get_installation_path(self)

        self.variables['GROMACS_ROOT'] = f'{install_path}'
        self.variables['PATH'] = r'$GROMACS_ROOT/bin:$PATH'

        self.executable = 'gmx_mpi' if 'mpi' in self.variant else 'gmx'
        nb_type = 'gpu' if 'cuda' in self.variant else 'cpu'
        self.executable_opts = ['mdrun', '-dlb yes',
                                f'-ntomp {self.num_cpus_per_task}', '-npme -1',
                                f'-nb {nb_type}', f'-bonded {nb_type}', '-s herflat.tpr']
