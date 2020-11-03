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


SPACK_ROOT_PROD = '/scratch/snx3000tds/hvictor/reframe-spack-tests'

SPACK = [
    # ['develop', 'stage'],
    ['develop', 'output'],
    # ['develop', 'prod'],
    # ['6fc6f1ea', 'stage'],
    # ['6fc6f1ea', 'prod'],
    # ['0.15.3'],
]

SPECS = [
    ['pkg-config', 'default'],
    ['zlib', 'default'],
    ['gromacs', '~mpi ~shared build_type=Release ^fftw ~mpi']
]

PKG_SPECS = ([d,m] for d in SPECS for m in SPACK)
PKG_SPECS = ([x[0][0], x[0][1], x[1][0], x[1][1]] for x in PKG_SPECS)


def parent_part_env(part, env):
    def _parent_part_env(src, dst):
        if dst:
            return dst[0].split(':')[1] == part and dst[1] == env
        return False
    return _parent_part_env


def clone_or_download_spack(spack_version):
    if spack_version in ['develop', 'latest'] or len(spack_version) > 7:
        cmds = [
            f'git clone https://github.com/spack/spack.git spack-{spack_version}',
            f'git clone https://github.com/eth-cscs/production.git',
            f'cd spack-{spack_version}'
        ]
        if len(spack_version) > 7:
            cmds += [
                f'git checkout -b {spack_version} {spack_version}'
            ]
    else:
        cmds = [
            f'wget https://github.com/spack/spack/releases/download/v{spack_version}/spack-{spack_version}.tar.gz',
            f'git clone https://github.com/eth-cscs/production.git',
            f'tar xf spack-{spack_version}.tar.gz',
            f'cd spack-{spack_version}'
        ]
    return cmds


@rfm.parameterized_test(*SPACK)
class SpackInstallSpack(rfm.RunOnlyRegressionTest):
    def __init__(self, spack_version, spack_path):
        self.spack_version = spack_version
        self.spack_path = spack_path

        self.valid_systems = ['dom:login']
        self.valid_prog_environs = ['builtin']

        self.executable = 'spack'
        # self.executable_opts = ['install', 'zlib']
        # self.executable_opts = ['compiler', 'find', '--scope', 'site']
        self.executable_opts = ['list']

        self.maintainers = ['VH']
        self.tags = {'spack'}

    @rfm.run_after('setup')
    def config_spack(self):
        if self.spack_path in ['stage']:
            self.spack_path = self.stagedir
        elif self.spack_path in ['output']:
            self.spack_path = self.outputdir
        elif self.spack_path in ['prod']:
            self.spack_path = SPACK_ROOT_PROD
        else:
            self.spack_path = spack_path

        self.variables = {
            'SPACK_ROOT': f'{self.spack_path}/spack-{self.spack_version}',
            'PATH': '$SPACK_ROOT/bin:$PATH',
        }

        self.prerun_cmds = [
            f'cd {self.spack_path}'
        ]
        self.prerun_cmds += clone_or_download_spack(self.spack_version)
        self.prerun_cmds += [
            f'spack compiler find --scope site',
            f'spack external find',
            f'mv ~/.spack/packages.yaml {self.spack_path}/production/',
            # https://github.com/eth-cscs/production/tree/master/spack/daint
            # removing recursive link ($spack/lib/spack/docs/_spack_root) from spack folder
            # I hate Spack! Who does that?
            f'rm -rf {self.spack_path}/spack-{self.spack_version}/lib/spack/docs'
            # not using "external find" dom because it breaks my setup
            # since external find ALWAYS write the ~/.spack/packages.yaml
            # Why??
            #f'spack external find',
        ]

    @rfm.run_after('setup')
    def set_sanity_patterns(self):
        patterns = [
            sn.assert_not_found('Error', self.stderr),
            sn.assert_not_found('Errno', self.stderr),
            sn.assert_not_found('See build log for details', self.stderr),
            # sn.assert_found(f'Added \\d+ new compilers to {self.spack_path}/spack-{self.spack_version}/etc/spack/compilers.yaml', self.stdout),
            sn.assert_found('Compilers are defined in the following files', self.stdout)
        ]
        if self.spack_path in [self.stagedir, self.outputdir]:
            patterns += [sn.assert_found(f'Added \\d+ new compilers to {self.spack_path}/spack-{self.spack_version}/etc/spack/compilers.yaml', self.stdout)]
        self.sanity_patterns = sn.all(patterns)



@rfm.parameterized_test(*PKG_SPECS)
# should be rfm.RegressionTest but simplifying for the test
class SpackPackageInstall(rfm.RunOnlyRegressionTest):
    def __init__(self, pkg, spec, spack_version, spack_path):
        self.pkg = pkg
        if spec in [None, 'default']:
            self.spec = ''
        else:
            self.spec = spec

        if spack_path in ['stage']:
            self.spack_path = None
        elif spack_path in ['prod']:
            self.spack_path = SPACK_ROOT_PROD
        else:
            self.spack_path = spack_path

        self.valid_systems = ['dom:gpu']
        self.valid_prog_environs = ['builtin']

        self.executable = 'spack'
        self.executable_opts = ['install', f'{self.pkg}', f'{self.spec}', f'%{compiler}']

        self.num_tasks = 1
        self.num_tasks_per_node = 1
        self.num_cpus_per_task = os.cpu_count()

        self.dep_name = f'SpackInstallSpack_{util.toalphanum(spack_version)}_{util.toalphanum(spack_path)}'
        self.depends_on(self.dep_name, when=parent_part_env('gpu', 'builtin'))

        self.maintainers = ['VH']
        self.tags = {'spack'}

    @rfm.run_after('setup')
    def config_spack(self):
        target = self.getdep(self.dep_name, 'builtin')
        self.variables = target.variables

    @rfm.run_after('setup')
    def set_sanity_patterns(self):
        self.sanity_patterns = sn.all([
            sn.assert_not_found('Error', self.stderr),
            sn.assert_not_found('Errno', self.stderr),
            sn.assert_not_found('See build log for details', self.stderr),
        ])
