import reframe as rfm
import reframe.utility.sanity as sn
import reframe.utility.typecheck as typ


class StreamTest(rfm.RegressionTest, abstract=True):
    template_param('cores', int)
    template_param('flags', typ.List[str])
    perf_param('Copy', 'MB/s')
    perf_param('Scale', 'MB/s')
    perf_param('Add', 'MB/s')
    perf_param('Triad', 'MB/s')

    def __init__(self):
        self.prebuild_cmds = [
            'wget http://www.cs.virginia.edu/stream/FTP/Code/stream.c',
        ]
        self.build_system = 'SingleSource'
        self.sourcepath = 'stream.c'
        self.build_system.cppflags = ['-DSTREAM_ARRAY_SIZE=$((1 << 25))']
        self.sanity_patterns = sn.assert_found(r'Solution Validates',
                                               self.stdout)
        self.perf_patterns = {
            'Copy': sn.extractsingle(r'Copy:\s+(\S+)\s+.*',
                                     self.stdout, 1, float),
            'Scale': sn.extractsingle(r'Scale:\s+(\S+)\s+.*',
                                      self.stdout, 1, float),
            'Add': sn.extractsingle(r'Add:\s+(\S+)\s+.*',
                                    self.stdout, 1, float),
            'Triad': sn.extractsingle(r'Triad:\s+(\S+)\s+.*',
                                      self.stdout, 1, float)
        }

    @rfm.run_before('compile')
    def setflags(self):
        environ = self.current_environ.name
        self.build_system.cflags = self.flags

    @rfm.run_before('run')
    def set_num_threads(self):
        self.num_cpus_per_task = self.cores
        self.variables = {
            'OMP_NUM_THREADS': str(self.cores),
            'OMP_PLACES': 'cores'
        }


class SpecializedStreamTest(StreamTest):
    def __init__(self):
        self.valid_systems = ['daint:gpu']
        self.valid_prg_environs = ['PrgEnv-gnu']

    @rfm.run_after('setup')
    def prepare_test(self):
        self.set_param('cores',  4, when=curr_sys('catalina:default'))
        self.set_param('cores', 12, when=curr_sys('daint:gpu'))
        self.set_param('cores', 36, when=curr_sys('daint:mc'))
        self.set_param('cores', 10, when=curr_sys('daint:login'))

        self.set_param( 'flags', ['-fopenmp', '-O3', '-Wall'], when=curr_env('cray'))
        self.set_param('flags', ['-fopenmp', '-O3', '-Wall'], when=curr_env('gnu'))
        self.set_param('flags', ['-qopenmp', '-O3', '-Wall'], when=curr_env('intel'))
        self.set_param('flags', ['-mp', '-O3'], when=curr_env('pgi'))

        self.add_reference('Copy',  25200, -0.05, 0.05, when=curr_sys('catalina'))
        self.add_reference('Scale', 25200, -0.05, 0.05, when=curr_sys('catalina'))
        self.add_reference('Add',   25200, -0.05, 0.05, when=curr_sys('catalina'))
        self.add_reference('Triad', 25200, -0.05, 0.05, when=curr_sys('catalina'))


class CompileStreamTest(StreamTest):
    def __init__(self):
        self.valid_systems = ['daint:login']
        self.valid_prg_environs = ['PrgEnv-gnu']
        self.skip_phase = 'run'

class RunStreamTest(StreamTest):
    def __init__(self):
        self.valid_systems = ['daint:gpu']
        self.valid_prg_environs = ['PrgEnv-gnu']
        self.skip_phase = 'compile'
