# Copyright 2016-2021 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import os

import reframe as rfm
import reframe.utility.sanity as sn

from reframe.core.deferrable import deferrable, _DeferredExpression

@deferrable
def invert(x):
    return 1.0 / x

@rfm.simple_test
class NAMDCPUExplorationCheck(rfm.RunOnlyRegressionTest):
    valid_systems = ['eiger:mc', 'pilatus:mc']
    valid_prog_environs = ['cpeGNU']
    modules = ['NAMD']
    executable = 'namd2'
    time_limit = '20m'
    #num_nodes = paramater([1, 6, 16])
    multithreading = parameter([True, False])
    nnodes = parameter([1, 6, 16])
    ntasks_per_node = parameter([1, 8, 16, 128, 256])
    affinity = parameter([True, False])
    idlepoll = parameter([True, False])
    threads = parameter([True, False])
    pemap = parameter([True, False])
    commap = parameter([True, False])
    pppn = parameter([True, False])
    # ppppn = parameter([True, False])
    isomalloc = parameter([True, False])
    #+setcpuaffinity +pemap 0-127:32.31 +commap 31-127:32'
    maintainers = ['CB', 'LM']
    tags = {'scs', 'external-resources'}
    extra_resources = {
        'switches': {
            'num_switches': 1
        }
    }

    @run_after('init')
    def set_resources_dir(self):
        self.sourcesdir = os.path.join(self.current_system.resourcesdir, 'NAMD', 'prod')
        self.readonly_files = ['par_all27_prot_na.inp','stmv.namd','stmv.pdb','stmv.psf']

    @run_before('run')
    def set_execution_opts(self):
        self.executable_opts = []
        if self.affinity:
            self.executable_opts.append('+setcpuaffinity')

        if self.idlepoll:
            self.executable_opts.append('+idlepoll')

        if self.isomalloc:
            self.executable_opts.append('+isomalloc_sync')

        self.use_multithreading = self.multithreading
        if self.multithreading:
            if self.threads:
                self.job.options = [f'--threads-per-core=2']
                self.num_tasks_per_core = None
            else:
                self.num_tasks_per_core = 2
            if self.current_system.name in ['eiger', 'pilatus']:
                self.max_tasks_per_node = 256
            else:
                self.max_tasks_per_node = 72
        else:
            self.num_tasks_per_core = 1
            if self.current_system.name in ['eiger', 'pilatus']:
                self.max_tasks_per_node = 128
            else:
                self.max_tasks_per_node = 36

        self.num_nodes = self.nnodes
        self.num_tasks_per_node = self.ntasks_per_node

        # skip it is not possible to accomodate the number of tasks per node
        self.skip_if(self.num_tasks_per_node > self.max_tasks_per_node)

        self.num_tasks = self.num_nodes * self.num_tasks_per_node
        self.num_cpus_per_task = self.max_tasks_per_node // self.num_tasks_per_node

        min_cpus_per_task = max(self.num_cpus_per_task-1, 1)
        if self.pppn:
            self.executable_opts += ['+ppn', f'{min_cpus_per_task}']

        # if self.ppppn:
        #     self.executable_opts += ['++ppn', f'{self.num_tasks_per_node}']

        if self.pemap:
            self.executable_opts += [
                '+pemap',
                f'0-{self.num_tasks_per_node}:{self.num_cpus_per_task}.{min_cpus_per_task}'
            ]

        if self.commap:
            self.executable_opts += [
                '+commap',
                f'{min_cpus_per_task}-{self.num_tasks_per_node}:{self.num_cpus_per_task}.1'
            ]

        self.executable_opts.append('stmv.namd')

    @run_after('init')
    def set_sanity_patterns(self):
        energy = sn.avg(sn.extractall(
            r'ENERGY:([ \t]+\S+){10}[ \t]+(?P<energy>\S+)',
            self.stdout, 'energy', float)
        )
        energy_reference = -2451359.5
        energy_diff = sn.abs(energy - energy_reference)
        self.sanity_patterns = sn.all([
            sn.assert_eq(sn.count(sn.extractall(
                         r'TIMING: (?P<step_num>\S+)  CPU:',
                         self.stdout, 'step_num')), 50),
            sn.assert_lt(energy_diff, 2720)
        ])

    @run_after('init')
    def set_generic_perf_references(self):
        self.reference.update({'*': {
            'ns_days': (0, None, None, 'ns/day')
        }})

    @run_after('init')
    def set_perf_patterns(self):
        self.perf_patterns = {
            'ns_days': invert(sn.avg(sn.extractall(
                r'Info: Benchmark time: \S+ CPUs \S+ '
                r's/step (?P<days_ns>\S+) days/ns \S+ MB memory',
                self.stdout, 'days_ns', float)))
        }

@rfm.simple_test
class NAMDCrayXCCPUExplorationCheck(NAMDCPUExplorationCheck):
    valid_systems = ['daint:mc', 'dom:mc']
    valid_prog_environs = ['builtin']
    ntasks_per_node = parameter([1, 4, 9, 36, 72])

