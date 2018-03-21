import itertools
import os

import reframe.utility.sanity as sn
from reframe.core.pipeline import RunOnlyRegressionTest


class PROMDBaseCheck(RunOnlyRegressionTest):
    def __init__(self, name, output_file, **kwargs):
        super().__init__(name, os.path.dirname(__file__), **kwargs)

        self.valid_prog_environs = ['PrgEnv-gnu']
        self.executable = 'promd'

        # Reset sources dir relative to the SCS apps prefix
        #self.sourcesdir = os.path.join(self.current_system.resourcesdir,
        #                               'PROMD')
        self.keep_files = [output_file]

        self.sanity_patterns = sn.assert_found('Finished mdrun', output_file)

        self.perf_patterns = {
            'perf': sn.extractsingle(r'Performance:\s+(?P<perf>\S+)',
                                     output_file, 'perf', float)
        }

        self.modules = ['GROMACS']
        self.maintainers = ['VH']
        self.strict_check = False
        self.use_multithreading = False
        self.extra_resources = {
            'switches': {
                'num_switches': 1
            }
        }



class PROMDCPUCheck(PROMDBaseCheck):
    def __init__(self, variant, **kwargs):
        super().__init__('gromacs_cpu_%s_check' % variant,
                         'md.log', **kwargs)

        self.valid_systems = ['daint:mc', 'dom:mc']
        self.descr = 'GROMACS CPU check'
        self.executable_opts = ('mdrun -dlb yes -ntomp 1 -npme -1 '
                                '-nb cpu -s herflat.tpr ').split()

        if self.current_system.name == 'dom':
            self.num_tasks = 216
            self.num_tasks_per_node = 36
        else:
            self.num_tasks = 576
            self.num_tasks_per_node = 36


class PROMDCPUProdCheck(PROMDCPUCheck):
    def __init__(self, **kwargs):
        super().__init__('prod', **kwargs)
        self.tags |= {'production'}
        self.reference = {
            'dom:mc': {
                'perf': (38.0, -0.05, None)
            },
            'daint:mc': {
                'perf': (73.0, -0.50, None)
            },
        }

def _get_checks(**kwargs):
    return [PROMDCPUProdCheck(**kwargs)]
