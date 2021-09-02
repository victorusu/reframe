# Copyright 2016-2021 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import contextlib
import os

import reframe as rfm
import reframe.utility.sanity as sn
import reframe.utility.udeps as udeps

scientific_systems = ['argon5832', 'nonanol_vacuo']
scientific_systems_energy_refences = {
    # every system has a different reference energy and drift
    'argon': (-443246, 5.0E-05),
    'nonanol': (-234188, 1.0E-04),
}

mpi_code_support = ['mpi', 'nompi']
gpupu_code_support = ['CUDA', 'OpenCL', 'OFF']
# supported_prgenvs = ['cpeGNU', 'PrgEnv-gnu', 'cpeCray', 'PrgEnv-cray']
supported_prgenvs = ['cpeGNU', 'PrgEnv-gnu']
supported_cpu_only_systems = ['daint:mc', 'dom:mc']
supported_gpgpu_systems = ['daint:gpu', 'dom:gpu']

download_files_test_environ_name = 'builtin'
download_files_test_part_name = 'login'

def parent_env_is(name):
    def _parent_env_is(src, dst):
        # print('name: ', name)
        # print('src: ', src)
        # print('dst: ', dst)
        if dst:
            return dst[1] == name

        return False

    return _parent_env_is


@rfm.simple_test
class CompileGROMACSWithNBLIB(rfm.RegressionTest):
    mpi = parameter(mpi_code_support)
    gpgpu = parameter(gpupu_code_support)
    maintainers = ['VH']
    tags = {'benchmark', 'prace', 'nblib'}
    valid_systems = supported_gpgpu_systems + supported_cpu_only_systems
    valid_prog_environs = supported_prgenvs

    # sourcesdir = 'https://gitlab.com/gromacs/gromacs.git'
    sourcesdir = 'https://gitlab.com/gromacs/nb-lib.git'
    build_system = 'CMake'
    build_locally = False
    build_time_limit='1h'

    @run_after('init')
    def filter_compilation_variants(self):
        if self.gpgpu in ['CUDA', 'OpenCL']:
            self.valid_systems = supported_gpgpu_systems
#
    @run_before('compile')
    def change_to_joes_branch(self):
        self.prebuild_cmds = [
            'git checkout -b "argon-script-gpu" "origin/argon-script-gpu"',
        ]

    @run_before('compile')
    def set_modules(self):
        self.modules = ['EasyBuild-custom/cscs', 'CMake/3.19.6']
        if self.gpgpu == 'CUDA':
            self.modules += ['cudatoolkit/11.2.0_3.36-7.0.2.1_2.2__g3ff9ab1', 'gcc/9.3.0']

    @run_before('compile')
    def set_build_job_options(self):
        self.num_cpus = self.current_partition.processor.num_cpus
        if not self.num_cpus:
            self.num_cpus = 12
        self.build_job.options = [f'--cpus-per-task={self.num_cpus}']
        self.build_system.max_concurrency = self.num_cpus

    @run_before('compile')
    def set_compilation_opts(self):
        self.build_system.builddir = f'{self.stagedir}/build_{self.mpi}_{self.gpgpu}'
        self.install_prefix = f'{self.outputdir}/gromacs_install_{self.mpi}_{self.gpgpu}'
        self.bin_prefix = f'{self.install_prefix}/bin'
        self.build_system.config_opts = [
            '-DCMAKE_BUILD_TYPE=Release',
            '-DCMAKE_POSITION_INDEPENDENT_CODE=ON',
            f'-DCMAKE_INSTALL_PREFIX={self.install_prefix}',
            '-DGMX_BUILD_OWN_FFTW=ON',
            '-DGMX_OPENMP=ON',
            f'-DGMX_GPU={self.gpgpu}',
            '-DBUILD_SHARED_LIBS=ON',
            '-DGMX_SIMD=AUTO',
            '-DGMX_CYCLE_SUBCOUNTERS=ON',
            '-DGMX_INSTALL_NBLIB_API=ON',
            '-DGMX_INSTALL_LEGACY_API=ON'
        ]
        if self.mpi == 'mpi':
            self.build_system.config_opts += ['-DGMX_MPI=ON']
        else:
            self.build_system.config_opts += ['-DGMX_MPI=OFF']

        if self.gpgpu == 'CUDA':
            self.build_system.config_opts += ['-DCUDA_TOOLKIT_ROOT_DIR=${CUDATOOLKIT_HOME}']
        # elif self.gpgpu == 'OpenCL':
        #     self.build_system.config_opts += ['-DCUDA_TOOLKIT_ROOT_DIR=${CUDATOOLKIT_HOME}']

    @run_before('compile')
    def set_post_compilation_cmds(self):
        self.postbuild_cmds = [f'make install']
        # nonmpi_build = self.build_system.emit_build_commands(self.current_environ)
        # for line in nonmpi_build:
        #     self.postbuild_cmds += [line.replace('-DGMX_MPI=ON', '-DGMX_MPI=OFF')]
        # self.postbuild_cmds += [f'make install']


    @run_before('run')
    def set_executable_and_opts(self):
        if self.mpi == 'mpi':
            self.executable = f'{self.bin_prefix}/gmx_mpi'
        else:
            self.executable = f'{self.bin_prefix}/gmx'

        self.executable_opts = ['mdrun','--version']

    @sanity_function
    def set_sanity_patterns(self):
        mpi_support = sn.extractsingle(r'MPI library:\s+(?P<mpi>\S+)',
                                self.stdout, 'mpi', str, item=-1)
        openmp_support = sn.extractsingle(r'OpenMP support:\s+(?P<openmp>\S+)',
                                self.stdout, 'openmp', str, item=-1)
        gpgpu_support = sn.extractsingle(r'GPU support:\s+(?P<gpu>\S+)',
                                self.stdout, 'gpu', str, item=-1)
        simd_support = sn.extractsingle(r'SIMD instructions:\s+(?P<simd>\S+)',
                                self.stdout, 'simd', str, item=-1)

        sanities = [
            sn.assert_false('disabled' in openmp_support),
            sn.assert_true('enabled' in openmp_support),
            sn.assert_false('NONE' in simd_support),
            sn.assert_found(f'Set runtime path of .*{self.executable}.*', self.stdout),
            sn.assert_not_found(f'error', self.stderr),
            sn.assert_not_found(f'Error', self.stderr),
            sn.assert_not_found(f'ERROR', self.stderr),
        ]
        if self.gpgpu == 'CUDA':
            sanities += [sn.assert_eq(gpgpu_support, 'CUDA')]
        elif self.gpgpu == 'OpenCL':
            sanities += [sn.assert_eq(gpgpu_support, 'OpenCL')]

        if self.mpi == 'mpi':
            sanities += [sn.assert_eq(mpi_support, 'MPI')]
        else:
            sanities += [sn.assert_eq(mpi_support, 'thread_mpi')]

        return sn.all([sanities])


@rfm.simple_test
class DownloadTestSystems(rfm.RunOnlyRegressionTest):
    benchmark = parameter(scientific_systems)
    maintainers = ['VH']
    tags = {'benchmark', 'prace', 'nblib'}
    valid_systems = [f'daint:{download_files_test_part_name}',
                     f'dom:{download_files_test_part_name}']
    valid_prog_environs = [download_files_test_environ_name]
    executable = 'echo Done'

    @run_after('init')
    def download_files(self):
        benchmark_url = 'https://gitlab.com/gromacs/gromacs/-/raw/master/src/testutils/simulationdatabase'
        self.keep_files = [
            f'{self.benchmark}.gro',
            f'{self.benchmark}.ndx',
            f'{self.benchmark}.top',
        ]

        self.prerun_cmds = []
        for files_to_keep in self.keep_files:
            # using curl instead of wget because it should, in principle, be present everywhere
            self.prerun_cmds += [f'curl -LJO {benchmark_url}/{files_to_keep}']

    @sanity_function
    def set_sanity_patterns(self):
        allfiles = os.listdir(f'{self.stagedir}')

        sanities = []
        for files_to_keep in self.keep_files:
            sanities += [sn.assert_true(files_to_keep in allfiles)]

        return sn.all([sanities])


@rfm.simple_test
class RunGROMACSWithNBLIB(rfm.RunOnlyRegressionTest):
    benchmark = parameter(scientific_systems)
    mpi = parameter(mpi_code_support)
    gpgpu = parameter(gpupu_code_support)
    multithreading = parameter(['no_multithreading', 'multithreading'])
    pin = parameter(['pin_on', 'pin_off'])
    nsteps = parameter([100, 1000, 10000, 100000])

    ener_ref = scientific_systems_energy_refences
    output_file = 'md.log'

    maintainers = ['VH']
    tags = {'benchmark', 'prace', 'nblib'}
    valid_systems = supported_gpgpu_systems + supported_cpu_only_systems
    valid_prog_environs = supported_prgenvs

    extra_resources = {
        'switches': {
            'num_switches': 1
        }
    }
    strict_check = False

    @run_after('init')
    def filter_compilation_variants(self):
        if self.gpgpu in ['CUDA', 'OpenCL']:
            self.valid_systems = supported_gpgpu_systems

    @run_after('init')
    def set_dependencies(self):
        self.download_files_test_name = f'DownloadTestSystems_{self.benchmark}'
        self.depends_on(self.download_files_test_name,
                        how=parent_env_is(download_files_test_environ_name))

        self.compiler_test_name = f'CompileGROMACSWithNBLIB_{self.mpi}_{self.gpgpu}'
        self.depends_on(self.compiler_test_name)

        if self.mpi == 'mpi':
            self.nonmpi_compiler_test_name = f'CompileGROMACSWithNBLIB_nompi_{self.gpgpu}'
            self.depends_on(self.nonmpi_compiler_test_name)

    @run_after('init')
    def set_files_to_keep(self):
        self.keep_files = [self.output_file]

    @run_before('setup')
    def set_multithreading(self):
        self.use_multithreading = True if self.multithreading == 'multithreading' else False

    @run_after('setup')
    def set_executable_cmds(self):
        compiler_test = self.getdep(self.compiler_test_name)
        self.executable = compiler_test.executable

        grompp_executable = self.executable
        if self.mpi == 'mpi':
            nonmpi_compiler_test = self.getdep(self.nonmpi_compiler_test_name)
            grompp_executable = nonmpi_compiler_test.executable

        grompp_cmdline, tprfile = self.construct_grompp_cmdline()
        # grompp executable options
        # this should include the `gmx` executable
        self.prerun_cmds = [
            grompp_executable + ' ' + grompp_cmdline,
        ]

        # mdrun executable options
        # this should NNOT include the `gmx` executable because it is defined by
        # self.executable
        self.executable_opts = ['mdrun', '-s', tprfile, '-g', self.output_file]
        if self.pin == 'pin_on':
            self.executable_opts += ['-pin', 'on']
        else:
            self.executable_opts += ['-pin', 'off']

    @sanity_function
    def set_sanity_patterns(self):
        # Conserved energy drift: -6.70e-06 kJ/mol/ps per atom
        energy = sn.extractsingle(r'Conserved energy drift:\s+(?P<energy>\S+)',
                                  self.output_file, 'energy', float, item=-1)
        return sn.all([
            sn.assert_found('Finished mdrun', self.output_file),
            sn.assert_lt(sn.abs(energy), 1.0e-3)
        ])

    @performance_function('ns/day')
    def min_perf(self, nid=None):
        '''Lowest performance recorded.'''
        return sn.extractsingle(r'Performance:\s+(?P<perf>\S+)',
                                     self.output_file, 'perf', float)

    def construct_grompp_cmdline(self):
        download_files_test = self.getdep(self.download_files_test_name,
                                  environ=download_files_test_environ_name,
                                  part=download_files_test_part_name)
        mdp_file_name = self.create_mdp_input_file()
        tprfile = os.path.join(download_files_test.stagedir,
                               f'{self.benchmark}.tpr')

        grompp_cmdline = ['grompp', '-o', tprfile, '-f', mdp_file_name]
        for inputfile in download_files_test.keep_files:
            if inputfile.endswith('.gro'):
                grompp_cmdline += ['-c',
                                   os.path.join(download_files_test.stagedir,
                                                inputfile)]
            elif inputfile.endswith('.ndx'):
                grompp_cmdline += ['-n',
                                   os.path.join(download_files_test.stagedir,
                                                inputfile)]
            elif inputfile.endswith('.top'):
                grompp_cmdline += ['-p',
                                   os.path.join(download_files_test.stagedir,
                                                inputfile)]

        return ' '.join(grompp_cmdline), tprfile

    def create_mdp_input_file(self):
        mdp_file_name = os.path.join(self.outputdir, 'grompp.mdp')
        mdp_input = [('integrator', 'md'),
                     ('dt', 0.0001),
                     ('cutoff-scheme', 'Verlet'),
                     ('nsteps', self.nsteps),
                     ('nstxout', 0),
                     ('nstvout', 0),
                     ('nstfout', 0)]

        mdp_input = '\n'.join([' = '.join([str(item) for item in kvpair]) for kvpair in mdp_input])
        with open(mdp_file_name, 'w') as fh:
            fh.write(mdp_input)
            fh.write('\n')

        return mdp_file_name

# @rfm.simple_test
class RunNBLIBBenchmark(rfm.RegressionTest):
    benchmark = parameter(scientific_systems)
    gpgpu = parameter(gpupu_code_support)
    multithreading = parameter(['no_multithreading', 'multithreading'])
    nsteps = parameter([100, 1000, 10000, 100000])

    ener_ref = scientific_systems_energy_refences

    maintainers = ['VH']
    tags = {'benchmark', 'prace', 'nblib'}
    valid_systems = supported_gpgpu_systems + supported_cpu_only_systems
    valid_prog_environs = supported_prgenvs

    extra_resources = {
        'switches': {
            'num_switches': 1
        }
    }
    strict_check = False

    @run_after('init')
    def filter_compilation_variants(self):
        if self.gpgpu in ['CUDA', 'OpenCL']:
            self.valid_systems = supported_gpgpu_systems

    @run_after('init')
    def set_dependencies(self):
        self.download_files_test_name = f'DownloadTestSystems_{self.benchmark}'
        self.depends_on(self.download_files_test_name,
                        how=parent_env_is(download_files_test_environ_name))

        self.nonmpi_compiler_test_name = f'CompileGROMACSWithNBLIB_nompi_{self.gpgpu}'
        self.depends_on(self.nonmpi_compiler_test_name)

    @run_before('setup')
    def set_multithreading(self):
        self.use_multithreading = True if self.multithreading == 'multithreading' else False

    @run_after('setup')
    def set_executable_cmds(self):
        compiler_test = self.getdep(self.compiler_test_name)
        self.executable = compiler_test.executable

        grompp_executable = self.executable
        if self.mpi == 'mpi':
            nonmpi_compiler_test = self.getdep(self.nonmpi_compiler_test_name)
            grompp_executable = nonmpi_compiler_test.executable

        grompp_cmdline, tprfile = self.construct_grompp_cmdline()
        # grompp executable options
        # this should include the `gmx` executable
        self.prerun_cmds = [
            grompp_executable + ' ' + grompp_cmdline,
        ]

        # mdrun executable options
        # this should NNOT include the `gmx` executable because it is defined by
        # self.executable
        self.executable_opts = ['mdrun', '-s', tprfile, '-g', self.output_file]
        if self.pin == 'pin_on':
            self.executable_opts += ['-pin', 'on']
        else:
            self.executable_opts += ['-pin', 'off']

        self.keep_files = [tprfile]

    @sanity_function
    def set_sanity_patterns(self):
        # Conserved energy drift: -6.70e-06 kJ/mol/ps per atom
        energy = sn.extractsingle(r'Conserved energy drift:\s+(?P<energy>\S+)',
                                  self.output_file, 'energy', float, item=-1)
        return sn.all([
            sn.assert_found('Finished mdrun', self.output_file),
            sn.assert_lt(sn.abs(energy), 1.0e-3)
        ])

    @performance_function('ns/day')
    def min_perf(self, nid=None):
        '''Lowest performance recorded.'''
        return sn.extractsingle(r'Performance:\s+(?P<perf>\S+)',
                                     self.output_file, 'perf', float)

    def construct_grompp_cmdline(self):
        download_files_test = self.getdep(self.download_files_test_name,
                                  environ=download_files_test_environ_name,
                                  part=download_files_test_part_name)
        mdp_file_name = self.create_mdp_input_file()
        # tprfile = os.path.join(self.stagedir, f'{self.benchmark}.tpr')
        tprfile = f'{self.benchmark}.tpr'

        grompp_cmdline = ['grompp', '-o', tprfile, '-f', mdp_file_name]
        for inputfile in download_files_test.keep_files:
            if inputfile.endswith('.gro'):
                grompp_cmdline += ['-c',
                                   os.path.join(download_files_test.stagedir,
                                                inputfile)]
            elif inputfile.endswith('.ndx'):
                grompp_cmdline += ['-n',
                                   os.path.join(download_files_test.stagedir,
                                                inputfile)]
            elif inputfile.endswith('.top'):
                grompp_cmdline += ['-p',
                                   os.path.join(download_files_test.stagedir,
                                                inputfile)]

        return ' '.join(grompp_cmdline), tprfile

    def create_mdp_input_file(self):
        mdp_file_name = os.path.join(self.outputdir, 'grompp.mdp')
        mdp_input = [('integrator', 'md'),
                     ('dt', 0.0001),
                     ('cutoff-scheme', 'Verlet'),
                     ('nsteps', self.nsteps),
                     ('nstxout', 0),
                     ('nstvout', 0),
                     ('nstfout', 0)]

        mdp_input = '\n'.join([' = '.join([str(item) for item in kvpair]) for kvpair in mdp_input])
        with open(mdp_file_name, 'w') as fh:
            fh.write(mdp_input)
            fh.write('\n')

        return mdp_file_name
