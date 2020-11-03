import functools
import reframe.core.fields as fields

import reframe.libstest.systems.daint as daint
import reframe.libstest.systems.daint as tsa

def _template_param(namespace, param, param_type):
    namespace[param] = fields.TypedField(param, param_type)
    print('setting template param on class')

# Proposal to use directives for the parameterization
# This introduces the concept of directives.
# Where classes can be templated based on sets of parameters
class Meta(type):
    @classmethod
    def __prepare__(cls, name, bases, **kwargs):
        namespace = super().__prepare__(cls, name, bases, **kwargs)
        namespace['template_param'] = functools.partial(_template_param, namespace)
        return namespace

# This is a minimified example of ReFrame's rfm.RegressionTest class
# Still needs to add the decorator magic to instantiate multiple classes
class RegressionTest(metaclass=Meta):
    parameterize('cores', int)

    def __init__(self):
        print('RegressionTest init')


# Library Test
class GROMACSBaseCheck(RegressionTest):
    parameterize('inputfile', ['heartflat', 'adenosine', 'membrane'])
    parameterize('num_tasks', [1, 2, 4, 6, 8, 12, 24])
    def __init__(self):

        self.valid_systems = ['daint:gpu','dom:gpu']
        self.valid_prog_environs = ['builtin']

        self.num_cpus_per_task = 1

        self.perf_patterns = {
            'perf': sn.extractsingle(r'Performance:\s+(?P<perf>\S+)',
                                     output_file, 'perf', float)
        }

        self.sourcesdir = 'https://github.com/victorusu/gromacs-inputfiles.git'

        self.strict_check = False
        self.use_multithreading = False

        self.sanity_patterns = sn.all([
            sn.assert_found('Finished mdrun', output_file),
        ])

        self.executable = 'gmx_mpi'
        self.executable_opts = ['mdrun', '-dlb yes',
                                f'-ntomp {self.num_cpus_per_task}', '-npme -1',
                                f'-s {inputfile}.tpr']

        self.tags = {'gromacs', 'external-resources'}


@rfm.simple_test
class GROMACSGPUCheck(GROMACSBaseCheck):
    def __init__(self):
        super().__init__()
        self.inputfile = self.inputfile + '.tpr'
        self.valid_systems = ['daint:gpu', 'dom:gpu']

        if self.num_tasks >= 12:
            self.num_tasks_per_node = 12
        self.num_cpus_per_task = 1

        self.executable_opts += [f'-nb {nb_type}']

@rfm.simple_test
class GROMACSCPUCheck(GROMACSBaseCheck):
    def __init__(self):
        self.inputfile = self.inputfile + '.tpr'
        self.valid_systems = ['daint:gpu', 'dom:gpu']

        if self.num_tasks >= 12:
            self.num_tasks_per_node = 12


@rfm.simple_test
class DDaintCheck(RegressionTest,DaintScalingTest,MPIPrgEnvs):
    # parameterize('num_tasks', [1])
    def __init__(self):
        # self.num_tasks = num_tasks
        # self.valid_systems = ['tsa:cn']


@rfm.simple_test
class DDaintCheck(RegressionTest,daint.DaintScalingTest,daint.MPIPrgEnvs):
    # parameterize('num_tasks', [1])
    def __init__(self):

    #     # self.num_tasks = num_tasks
    #     self.valid_systems = ['tsa:cn']

@rfm.simple_test
class DTsaCheck(RegressionTest,tsa.DaintScalingTest,tsa.MPIPrgEnvs):
    # parameterize('num_tasks', [1])
    def __init__(self):
        # self.num_tasks = num_tasks
        self.valid_systems = ['tsa:cn']


@rfm.parameterized_test(['PrgEnv-gnu'])
class DDaintCheck(RegressionTest,DaintScalingTest):
    parameterize('prg_env', ['PrgEnv-gnu'])
    def __init__(self):
        # self.num_tasks = num_tasks
        self.valid_systems = ['tsa:cn']


@rfm.simple_test
class DDomCheck(RegressionTest,DomScalingTest):
    # parameterize('num_tasks', [1])
    def __init__(self):
        # self.num_tasks = num_tasks


@rfm.simple_test
class DVictorCheck(RegressionTest):
    # parameterize('num_tasks', [1])
    def __init__(self):
        dom = DomScalingTest()
        diant = DaintScalingTest()
        self.num_tasks = dom.num_task + 2.6 * daint.num_task
        # self.num_tasks = num_tasks



reframe -G num_task=[10, 12, 24] -G nsteps=[100, 500] -n Gromacs

if __name__ == '__main__':
    # c = C()
    # c.cores = 10

    d = DCheck()
    d.eirini = 'Hi!'
    d.cores = 10
    # d.eirini = 10
#    c.cores = 'foo'
