# # Copyright 2016-2020 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# # ReFrame Project Developers. See the top-level LICENSE file for details.
# #
# # SPDX-License-Identifier: BSD-3-Clause

# import contextlib
# import itertools
# import os

# import reframe as rfm
# import reframe.utility.sanity as sn


# class GREASYGROMACSBaseCheck(rfm.RunOnlyRegressionTest):
#     def __init__(self, output_file, nnodes,  num_greasy_tasks, nworkes_per_node,
#                  nranks_per_worker, ncpus_per_worker):
#         self.valid_prog_environs = ['builtin']
#         self.executable = 'gmx_mpi'

#         # Reset sources dir relative to the SCS apps prefix
#         self.sourcesdir = os.path.join(self.current_system.resourcesdir,
#                                        'Gromacs', 'herflat')
#         self.keep_files = [output_file]


#         self.tasks_file = 'tasks-gromacs.txt'
#         self.executable_opts = [self.tasks_file]
#         self.greasy_logfile = 'greasy.log'
#         self.keep_files = [self.tasks_file, self.greasy_logfile]
#         self.use_multithreading = False
#         self.num_greasy_tasks = num_greasy_tasks
#         self.nworkes_per_node = nworkes_per_node
#         self.nranks_per_worker = nranks_per_worker
#         self.num_tasks_per_node = nranks_per_worker * nworkes_per_node
#         self.num_tasks = self.num_tasks_per_node * nnodes
#         self.num_cpus_per_task = ncpus_per_worker
#         self.sanity_patterns = self.eval_sanity()

#         # Reference value is system agnostic
#         # Adding 10 secs of slowdown per greasy tasks
#         # this is to compensate for whenever the systems are full and srun gets
#         # slightly slower
#         refperf = (
#             (self.sleep_time+10)*num_greasy_tasks / nworkes_per_node / nnodes
#         )
#         self.reference = {
#             '*': {
#                 'time': (refperf, None, 0.5, 's')
#             }
#         }
#         self.perf_patterns = {
#             'time': sn.extractsingle(r'Total time: (?P<perf>\S+)',
#                                      self.greasy_logfile,
#                                      'perf', to_seconds)
#         }
#         # On SLURM there is no need to set OMP_NUM_THREADS if one defines
#         # num_cpus_per_task, but adding for completeness and portability
#         self.variables = {
#             'OMP_NUM_THREADS': str(self.num_cpus_per_task),
#             'GREASY_NWORKERS_PER_NODE': str(nworkes_per_node),
#             'GREASY_LOGFILE': self.greasy_logfile
#         }


#         self.modules = ['GROMACS']
#         self.maintainers = ['VH', 'SO']
#         self.strict_check = False
#         self.use_multithreading = False
#         self.extra_resources = {
#             'switches': {
#                 'num_switches': 1
#             }
#         }
#         self.tags = {'scs', 'external-resources'}

#     @rfm.run_before('run')
#     def generate_tasks_file(self):
#         with open(os.path.join(self.stagedir, self.tasks_file), 'w') as fp:
#             for i in range(self.num_greasy_tasks):
#                 fp.write(f'./{self.executable} output-{i}\n')

#     @rfm.run_before('run')
#     def daint_dom_gpu_specific_workaround(self):
#         if self.current_partition.fullname in ['daint:gpu', 'dom:gpu']:
#             self.variables['CRAY_CUDA_MPS'] = '1'
#             self.variables['CUDA_VISIBLE_DEVICES'] = '0'
#             self.variables['GPU_DEVICE_ORDINAL'] = '0'
#             self.extra_resources = {
#                 'gres': {
#                     'gres': 'gpu:0,craynetwork:4'
#                 }
#             }
#         elif self.current_partition.fullname in ['dom:mc']:
#             self.extra_resources = {
#                 'gres': {
#                     'gres': 'craynetwork:72'
#                 }
#             }
#         elif self.current_partition.fullname in ['daint:mc']:
#             if self.variant != 'serial':
#                 self.extra_resources = {
#                     'gres': {
#                         'gres': 'craynetwork:72'
#                     }
#                 }

#     @rfm.run_before('run')
#     def change_executable_name(self):
#         # After compiling the code we can change the executable to be
#         # the greasy one
#         self.executable = 'greasy'

#     @rfm.run_before('run')
#     def set_launcher(self):
#         # The job launcher has to be changed to local since greasy
#         # make calls to srun
#         self.job.launcher = getlauncher('local')()

#     @sn.sanity_function
#     def eval_sanity(self):
#         output_files = []
#         output_files = [file for file in os.listdir(self.stagedir)
#                         if file.startswith('md-output')]
#         num_greasy_tasks = len(output_files)
#         failure_msg = (f'Requested {self.num_greasy_tasks} task(s), but '
#                        f'executed only {num_greasy_tasks} tasks(s)')
#         sn.evaluate(sn.assert_eq(num_greasy_tasks, self.num_greasy_tasks,
#                                  msg=failure_msg))

#         energy = sn.extractsingle(r'\s+Potential\s+Kinetic En\.\s+Total Energy'
#                                   r'\s+Conserved En\.\s+Temperature\n'
#                                   r'(\s+\S+){2}\s+(?P<energy>\S+)(\s+\S+){2}\n'
#                                   r'\s+Pressure \(bar\)\s+Constr\. rmsd',
#                                   output_file, 'energy', float, item=-1)
#         energy_reference = -3270799.9

#         self.sanity_patterns = sn.all([
#             sn.assert_found('Finished mdrun', output_file),
#             sn.assert_reference(energy, energy_reference, -0.001, 0.001)
#         ])

#         self.perf_patterns = {
#             'perf': sn.extractsingle(r'Performance:\s+(?P<perf>\S+)',
#                                      output_file, 'perf', float)
#         }

#         num_tasks = sn.getattr(self, 'nranks_per_worker')
#         num_cpus_per_task = sn.getattr(self, 'num_cpus_per_task')

#         def tid(match):
#             return int(match.group(1))

#         def num_threads(match):
#             return int(match.group(2))

#         def rank(match):
#             return int(match.group(3))

#         def num_ranks(match):
#             return int(match.group(4))

#         for output_file in output_files:
#             result = sn.findall(r'Hello, World from thread \s*(\d+) out '
#                                 r'of \s*(\d+) from process \s*(\d+) out of '
#                                 r'\s*(\d+)', output_file)

#             failure_msg = (f'Found {sn.count(result)} Hello, World... '
#                            f'pattern(s) but expected '
#                            f'{num_tasks * num_cpus_per_task} pattern(s) '
#                            f'inside the output file {output_file}')
#             sn.evaluate(sn.assert_eq(sn.count(result),
#                                      num_tasks * num_cpus_per_task,
#                                      msg=failure_msg))

#             sn.evaluate(sn.all(
#                 sn.chain(
#                     sn.map(
#                         lambda x: sn.assert_lt(
#                             tid(x), num_threads(x),
#                             msg=(f'Found {tid(x)} threads rather than '
#                                  f'{num_threads(x)}')
#                         ), result
#                     ),
#                     sn.map(
#                         lambda x: sn.assert_lt(
#                             rank(x), num_ranks(x),
#                             msg=(f'Rank id {rank(x)} is not lower than the '
#                                  f'number of ranks {self.nranks_per_worker} '
#                                  f'in output file')
#                         ), result
#                     ),
#                     sn.map(
#                         lambda x: sn.assert_lt(
#                             tid(x), self.num_cpus_per_task,
#                             msg=(f'Rank id {tid(x)} is not lower than the '
#                                  f'number of cpus per task '
#                                  f'{self.num_cpus_per_task} in output '
#                                  f'file {output_file}')
#                         ), result
#                     ),
#                     sn.map(
#                         lambda x: sn.assert_eq(
#                             num_threads(x), num_cpus_per_task,
#                             msg=(f'Found {num_threads(x)} threads rather than '
#                                  f'{self.num_cpus_per_task} in output file '
#                                  f'{output_file}')
#                         ), result
#                     ),
#                     sn.map(
#                         lambda x: sn.assert_lt(
#                             rank(x), num_tasks,
#                             msg=(f'Found {rank(x)} threads rather than '
#                                  f'{self.num_cpus_per_task} in output file '
#                                  f'{output_file}')
#                         ), result
#                     ),
#                     sn.map(
#                         lambda x: sn.assert_eq(
#                             num_ranks(x), num_tasks,
#                             msg=(f'Number of ranks {num_ranks(x)} is not '
#                                  f'equal to {self.nranks_per_worker} in '
#                                  f'output file {output_file}')
#                         ), result
#                     )
#                 )
#             ))
#         sn.evaluate(sn.assert_found(r'Finished greasing', self.greasy_logfile))
#         sn.evaluate(sn.assert_found(
#             (f'INFO: Summary of {self.num_greasy_tasks} '
#              f'tasks: '
#              f'{self.num_greasy_tasks} OK, '
#              f'0 FAILED, '
#              f'0 CANCELLED, '
#              fr'0 INVALID\.'), self.greasy_logfile
#         ))

#         return True


# @rfm.required_version('>=2.19')
# @rfm.parameterized_test(*([s, v]
#                           for s in ['small', 'large']
#                           for v in ['prod', 'maint']))
# class GREASYGROMACSGPUCheck(GREASYGROMACSBaseCheck):
#     def __init__(self, scale, variant):
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


# @rfm.required_version('>=2.19')
# @rfm.parameterized_test(*([s, v]
#                           for s in ['small', 'large']
#                           for v in ['prod']))
# class GromacsCPUCheck(GromacsBaseCheck):
#     def __init__(self, scale, variant):
#         super().__init__('md.log')
#         self.valid_systems = ['daint:mc']
#         self.descr = 'GROMACS CPU check'
#         self.executable_opts = ['mdrun', '-dlb yes', '-ntomp 1', '-npme -1',
#                                 '-nb cpu', '-s herflat.tpr']

#         if scale == 'small':
#             self.valid_systems += ['dom:mc']
#             self.num_tasks = 216
#             self.num_tasks_per_node = 36
#         else:
#             self.num_tasks = 576
#             self.num_tasks_per_node = 36

#         references = {
#             'prod': {
#                 'small': {
#                     'dom:mc': {'perf': (40.0, -0.05, None, 'ns/day')},
#                     'daint:mc': {'perf': (38.8, -0.10, None, 'ns/day')}
#                 },
#                 'large': {
#                     'daint:mc': {'perf': (68.0, -0.20, None, 'ns/day')}
#                 }
#             },
#         }
#         self.reference = references[variant][scale]
#         self.tags |= {'maintenance' if variant == 'maint' else 'production'}
