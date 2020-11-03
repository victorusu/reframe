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


def mkdirp(*paths, **kwargs):
    """Creates a directory, as well as parent directories if needed.
    Arguments:
        paths (str): paths to create with mkdirp
    Keyword Aguments:
        mode (permission bits or None, optional): optional permissions to set
            on the created directory -- use OS default if not provided
        group (group name or None, optional): optional group for permissions of
            final created directory -- use OS default if not provided. Only
            used if world write permissions are not set
        default_perms ('parents' or 'args', optional): The default permissions
            that are set for directories that are not themselves an argument
            for mkdirp. 'parents' means intermediate directories get the
            permissions of their direct parent directory, 'args' means
            intermediate get the same permissions specified in the arguments to
            mkdirp -- default value is 'args'
    """
    mode = kwargs.get('mode', None)
    group = kwargs.get('group', None)
    default_perms = kwargs.get('default_perms', 'args')

    for path in paths:
        if not os.path.exists(path):
            try:
                # detect missing intermediate folders
                intermediate_folders = []
                last_parent = ''

                intermediate_path = os.path.dirname(path)

                while intermediate_path:
                    if os.path.exists(intermediate_path):
                        last_parent = intermediate_path
                        break

                    intermediate_folders.append(intermediate_path)
                    intermediate_path = os.path.dirname(intermediate_path)

                # create folders
                os.makedirs(path)

                # leaf folder permissions
                if mode is not None:
                    os.chmod(path, mode)
                if group:
                    chgrp_if_not_world_writable(path, group)
                    if mode is not None:
                        os.chmod(path, mode)  # reset sticky grp bit post chgrp

                # for intermediate folders, change mode just for newly created
                # ones and if mode_intermediate has been specified, otherwise
                # intermediate folders list is not populated at all and default
                # OS mode will be used
                if default_perms == 'args':
                    intermediate_mode = mode
                    intermediate_group = group
                elif default_perms == 'parents':
                    stat_info = os.stat(last_parent)
                    intermediate_mode = stat_info.st_mode
                    intermediate_group = stat_info.st_gid
                else:
                    msg = "Invalid value: '%s'. " % default_perms
                    msg += "Choose from 'args' or 'parents'."
                    raise ValueError(msg)

                for intermediate_path in reversed(intermediate_folders):
                    if intermediate_mode is not None:
                        os.chmod(intermediate_path, intermediate_mode)
                    if intermediate_group is not None:
                        chgrp_if_not_world_writable(intermediate_path,
                                                    intermediate_group)
                        os.chmod(intermediate_path,
                                 intermediate_mode)  # reset sticky bit after

            except OSError as e:
                if e.errno != errno.EEXIST or not os.path.isdir(path):
                    raise e
        elif not os.path.isdir(path):
            raise OSError(errno.EEXIST, "File already exists", path)


TEMPLATE_SPECS = {
    #
    # GROMACS 2020
    #
    ('GROMACS@2020', 'CrayGNU@20.08') : {
        'spec' : 'mpi,openmp,ownfftw,~shared,cuda,~plumed,~pat',
        'builddependencies' : {
            ('CMake@3.14.5', 'system'): {},
            ('PLUMED@2.6.1', 'CrayGNU@20.08') : {
                'when' : 'plumed'
            }
        },
        'dependencies' : {
            ('zlib@1.2.11', 'CrayGNU@20.08'): {},
            ('GSL@2.5', 'CrayGNU@20.08'): {},
            # ('PLUMED@2.6.1', 'CrayGNU@20.08') : {
            #     'when' : 'plumed'
            # }
        },
        'variants' : ['~cuda', '~cuda,pat', 'pat', 'plumed,pat', '~cuda,plumed,pat']
    },
    ('GROMACS@2020.1', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2020', 'CrayGNU@20.08')
    },
    ('GROMACS@2020.2', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2020', 'CrayGNU@20.08'),
        # 'variants' : ['~cuda', '~cuda,pat', 'pat', 'plumed,pat', '~cuda,plumed,pat']
    },
    ('GROMACS@2020.3', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2020', 'CrayGNU@20.08')
    },
    ('GROMACS@2020.4', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2020', 'CrayGNU@20.08')
    },
    #
    # GROMACS 2019
    #
    ('GROMACS@2018', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2020', 'CrayGNU@20.08')
    },
    ('GROMACS@2018.1', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2018', 'CrayGNU@20.08')
    },
    ('GROMACS@2018.2', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2018', 'CrayGNU@20.08'),
    },
    ('GROMACS@2018.3', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2018', 'CrayGNU@20.08')
    },
    ('GROMACS@2018.4', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2018', 'CrayGNU@20.08')
    },
    ('GROMACS@2018.5', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2018', 'CrayGNU@20.08')
    },
    ('GROMACS@2018.6', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2018', 'CrayGNU@20.08')
    },
    ('GROMACS@2018.7', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2018', 'CrayGNU@20.08')
    },
    ('GROMACS@2018.8', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2018', 'CrayGNU@20.08')
    },
    #
    # GROMACS 2019
    #
    ('GROMACS@2019', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2020', 'CrayGNU@20.08')
    },
    ('GROMACS@2019.1', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2019', 'CrayGNU@20.08')
    },
    ('GROMACS@2019.2', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2019', 'CrayGNU@20.08'),
    },
    ('GROMACS@2019.3', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2019', 'CrayGNU@20.08')
    },
    ('GROMACS@2019.4', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2019', 'CrayGNU@20.08')
    },
    ('GROMACS@2019.5', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2019', 'CrayGNU@20.08')
    },
    ('GROMACS@2019.6', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2019', 'CrayGNU@20.08')
    },
    ('PLUMED@2.6.0', 'CrayGNU@20.08') : {
        'spec' : 'mpi,openmp',
        'dependencies' : {
            ('zlib@1.2.11', 'CrayGNU@20.08'): {},
            ('GSL@2.5', 'CrayGNU@20.08'): {}
        },
    },
    ('PLUMED@2.6.1', 'CrayGNU@20.08') : {
        'from' : ('PLUMED@2.6.0', 'CrayGNU@20.08')
    },
    ('PLUMED@2.5.0', 'CrayGNU@20.08') : {
        'spec' : 'mpi,openmp',
        'dependencies' : {
            ('zlib@1.2.11', 'CrayGNU@20.08'): {},
            ('GSL@2.5', 'CrayGNU@20.08'): {}
        },
    },
    ('PLUMED@2.5.1', 'CrayGNU@20.08') : {
        'from' : ('PLUMED@2.5.0', 'CrayGNU@20.08')
    },
    ('PLUMED@2.5.2', 'CrayGNU@20.08') : {
        'from' : ('PLUMED@2.5.0', 'CrayGNU@20.08')
    },
    ('PLUMED@2.5.3', 'CrayGNU@20.08') : {
        'from' : ('PLUMED@2.5.0', 'CrayGNU@20.08')
    },
    ('PLUMED@2.5.4', 'CrayGNU@20.08') : {
        'from' : ('PLUMED@2.5.0', 'CrayGNU@20.08')
    },
    ('PLUMED@2.5.5', 'CrayGNU@20.08') : {
        'from' : ('PLUMED@2.5.0', 'CrayGNU@20.08')
    },
    ('GSL@2.5', 'CrayGNU@20.08') : {
        'spec' : 'openmp,optarch,unroll',
    },
    ('zlib@1.2.11', 'CrayGNU@20.08') : {
        'spec' : 'shared',
    },
}


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



class Pkg():
    def __init__(self, name_and_toolchain, template_specs=None):
        if template_specs:
            self.template_specs = template_specs
        else:
            self.template_specs = TEMPLATE_SPECS

        self.name_and_toolchain = name_and_toolchain
        self.pkg_info = self._extract_pkg_info(self.name_and_toolchain)

        self.spec = self.pkg_info[name_and_toolchain]['spec'] if 'spec' in self.pkg_info[name_and_toolchain] else None
        self.dependencies = self.pkg_info[name_and_toolchain]['dependencies'] if 'dependencies' in self.pkg_info[name_and_toolchain] else None
        self.builddependencies = self.pkg_info[name_and_toolchain]['builddependencies'] if 'builddependencies' in self.pkg_info[name_and_toolchain] else None
        self.variants = self.pkg_info[name_and_toolchain]['variants'] if 'variants' in self.pkg_info[name_and_toolchain] else None
        self.template_file = self.pkg_info[name_and_toolchain]['template_file'] if 'template_file' in self.pkg_info[name_and_toolchain] else None

        self.name, self.version, self.toolchain, self.toolchain_version = self._get_info_from_name_and_toolchain(name_and_toolchain)

        if not self.template_file:
            self.template_file = f'{self.name}.j2'

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret  = f'Pkg({self.name_and_toolchain})'
        # ret  = f'Pkg({self.name_and_toolchain}):\n'
        # ret += f'                 spec: {self.spec}\n'
        # ret += f'         dependencies: {self.dependencies}\n'
        # ret += f'    builddependencies: {self.builddependencies}\n'
        # ret += f'             variants: {self.variants}\n'
        # ret += f'        template_file: {self.template_file}'
        return ret

    def _get_info_from_name_and_toolchain(self, name_and_toolchain):
        name = ''
        version = ''
        toolchain = ''
        toolchain_version = ''

        if not isinstance(name_and_toolchain, tuple):
            raise ValueError(f'Unable to retrieve info from {name_and_toolchain}')

        name_and_version = name_and_toolchain[0].split('/')
        if len(name_and_version) > 1:
            if name_and_toolchain[1] == 'external':
                toolchain = 'external'
                toolchain_version = 'external'
            else:
                raise ValueError(f'Unable to retrieve info from {name_and_toolchain}')

            name, version = name_and_version
            return name, version, toolchain, toolchain_version

        name_and_version = name_and_toolchain[0].split('@')
        name, version = name_and_version

        if len(name_and_version) < 2:
            raise ValueError(f'Unable to retrieve info from {name_and_toolchain}')

        try:
            if name_and_toolchain[1] == 'system':
                toolchain = 'system'
                toolchain_version = 'system'
            else:
                toolchain_and_version = name_and_toolchain[1].split('@')
                toolchain = toolchain_and_version[0]
                toolchain_version = toolchain_and_version[1]
        except:
            raise ValueError(f'Unable to retrieve info from {name_and_toolchain}')

        return name, version, toolchain, toolchain_version

    @classmethod
    def from_pkg_info(self, pkg_info):
        key = list(pkg_info.keys())[0]
        return Pkg(key, pkg_info)

    def _expand_pkg_info(self, name_and_toolchain):
        ret = {}
        if name_and_toolchain in self.template_specs:
            ret = self.template_specs[name_and_toolchain]
            if 'from' in ret:
                new = self._expand_pkg_info(ret['from'])
                new.update(ret)
                ret = new
            return ret
        else:
            return {}

    def _extract_pkg_info(self, name_and_toolchain):
        pkg_info = {}
        pkg_info[name_and_toolchain] = self._expand_pkg_info(name_and_toolchain)
        if 'from' in pkg_info[name_and_toolchain]:
            del pkg_info[name_and_toolchain]['from']
        return pkg_info

    def _replace_spec_with_variant(self, spec, variant):
        specs = spec.split(',')
        ret = ''
        for var in variant.split(','):
            if var.startswith('~'):
                opposite = var.replace('~', '')
            else:
                opposite = '~' + var

            specs = list(map(lambda x: x if x != opposite else var, specs))
        ret = ','.join(specs)
        return ret

    def _get_variants(self, spec):
        import copy
        pkg_info = copy.deepcopy(self.pkg_info)
        ret = [pkg_info]
        if 'variants' in pkg_info[spec] and 'spec' in pkg_info[spec]:
            for variant in pkg_info[spec]['variants']:
                new_pkg_info = copy.deepcopy(pkg_info)
                new_pkg_info[spec]['spec'] = self._replace_spec_with_variant(new_pkg_info[spec]['spec'], variant)
                ret += [new_pkg_info]

        return ret

    def get_variants(self):
        return self._get_variants(self.name_and_toolchain)

    def _get_deps(self, spec):
        node = self._extract_pkg_info(spec)
        ret = []
        if node:
            ret = [node]
        if 'builddependencies' in node[spec]:
            for builddep in node[spec]['builddependencies']:
                ret += self._get_deps(builddep)
        if 'dependencies' in node[spec]:
            for dep in node[spec]['dependencies']:
                ret += self._get_deps(dep)
        return ret

    def get_build_tree(self):
        return self._get_deps(self.name_and_toolchain)

    def _find_all_dependencies_with_name(self, name, dependencies):
        ret = []
        for key, value in dependencies.items():
            depname, _, _, _ = self._get_info_from_name_and_toolchain(key)
            if depname == name:
                ret += [(key, value)]

        return ret

    def _prune_dependencies(self, dependencies, spec_list):
        import copy
        deps = copy.deepcopy(dependencies)
        dep_list = list(deps.items() if deps else {})
        for dep in dep_list:
            key, value = dep
            dep_name, dep_version, _, _ = self._get_info_from_name_and_toolchain(key)
            add_dep = True
            if 'when' in value:
                for cond in value['when'].split(','):
                    if cond not in spec_list:
                        add_dep = False

            if not add_dep:
                del deps[key]

        return deps

    def expand_dependencies(self, dependencies):
        ret = []
        if not dependencies:
            return []

        for dep, value in dependencies.items():
            name, version, toolchain, toolchain_version = self._get_info_from_name_and_toolchain(dep)
            suffix = ''
            if 'suffix' in value:
                suffix = value['suffix']

            if 'system' == toolchain:
                suffix = suffix if suffix else ''
                ret += [f"('{name}', '{version}', '{suffix}', True),"]
            elif 'external' == toolchain:
                ret += [f"'({name}/{version}', EXTERNAL_MODULE),"]
            else:
                if suffix:
                    ret += [f"('{name}', '{version}', '{suffix}'),"]
                else:
                    ret += [f"('{name}', '{version}'),"]

        return ret

    def get_version_suffix_and_dependencies(self, spec_list, builddependencies, dependencies):
        version_suffix = ''
        builddeps = self._prune_dependencies(builddependencies, spec_list)
        deps = self._prune_dependencies(dependencies, spec_list)
        if 'plumed' in spec_list:
            plumed_version = ''
            if builddependencies:
                # builddeps = self._prune_dependencies('PLUMED', builddependencies, spec_list)
                plumed_build = self._find_all_dependencies_with_name('PLUMED', builddeps)
                if plumed_build:
                    _, plumed_version, _, _ = self._get_info_from_name_and_toolchain(plumed_build[0][0])

            if dependencies:
                # deps = self._prune_dependencies('PLUMED', dependencies, spec_list)
                plumed = self._find_all_dependencies_with_name('PLUMED', deps)
                if plumed:
                    _, plumed_version, _, _ = self._get_info_from_name_and_toolchain(plumed[0][0])

            if plumed_build or plumed:
                version_suffix += '-PLUMED'

            if plumed_version:
                version_suffix += f'-{plumed_version}'

        if 'cuda' in spec_list:
            version_suffix += '-cuda'

        if 'pat' in spec_list:
            version_suffix += '-pat'

        return version_suffix, builddeps, deps

    def to_easyconfig(self):
        builddeps = self.builddependencies
        deps = self.dependencies

        ec = {
            'name': self.name,
            'version': self.version,
            'toolchain': self.toolchain,
            'toolchain_version': self.toolchain_version,
        }

        version_suffix = ''
        if self.spec:
            spec_list = self.spec.split(',')
            version_suffix, builddeps, deps = self.get_version_suffix_and_dependencies(spec_list, builddeps, deps)

            for spec in spec_list:
                ec[spec] = spec

        ec['version_suffix'] = version_suffix
        ec['builddependencies'] = self.expand_dependencies(builddeps)
        ec['dependencies'] = self.expand_dependencies(deps)

        if ec['toolchain_version'] == 'system':
            filename = f"{ec['name']}-{ec['version']}.eb"
        else:
            filename = f"{ec['name']}-{ec['version']}-{ec['toolchain']}-{ec['toolchain_version']}{ec['version_suffix']}.eb"

        return filename, ec


def parent_part_env(part, env):
    def _parent_part_env(src, dst):
        if dst:
            return dst[0].split(':')[1] == part and dst[1] == env
        return False
    return _parent_part_env


EASYBUILD_PKGS = [
    [Pkg(spec)] for spec in TEMPLATE_SPECS
]


@rfm.parameterized_test(*EASYBUILD_PKGS)
class GenerateEBFiles(rfm.CompileOnlyRegressionTest):
    def __init__(self, pkg):
        self.pgk_list = [pkg]
        self.name = f'GenerateEBFiles_{pkg.name}_{pkg.version}_{pkg.toolchain}_{pkg.toolchain_version}'
        self.descr = f'Generate easyconfig files for {pkg.name}_{pkg.version}_{pkg.toolchain}_{pkg.toolchain_version}'
        self.valid_systems = ['daint:login', 'dom:login']
        self.valid_prog_environs = ['builtin']

        self.num_tasks = 1
        self.num_tasks_per_node = 1
        self.num_cpus_per_task = 1

        self.sourcesdir = 'src/eb_templates'
        self.sourcepath = 'deleteme.c'
        self.executable = 'deleteme'
        self.build_system = 'SingleSource'

        self.prebuild_cmds = [
            'echo "int main(){return 0;}" > deleteme.c'
        ]

        self.sanity_patterns = sn.assert_not_found('ERROR', self.stderr)

        self.maintainers = ['VH']
        self.tags = {'gromacs', 'easybuild'}

    def _process_pkg(self, pkg):
        for spec in pkg.get_build_tree():
            template_pkg = Pkg.from_pkg_info(spec)
            for variant in template_pkg.get_variants():
                pkg_variant = Pkg.from_pkg_info(variant)
                filename, ec_dict = pkg_variant.to_easyconfig()

                with open(os.path.join(self.stagedir, f'{pkg_variant.template_file}'), 'r') as fp:
                    ec_template = fp.read()

                import jinja2
                if pkg_variant.template_file:
                    easyconfig = jinja2.Template(ec_template)
                else:
                    raise ValueError(f'Unable to open template file {pkg_variant.template_file}')

                easybuild_dir = os.path.join(self.outputdir, 'easybuild','easyconfigs',f'{pkg_variant.name[0].lower()}',f'{pkg_variant.name}')
                mkdirp(easybuild_dir)

                with open(os.path.join(easybuild_dir, filename), 'w') as fp:
                    fp.write(easyconfig.render(ec_dict))
                    fp.write('\n')

    @rfm.run_after('compile')
    def generate_ebfile(self):
        for pkg in self.pgk_list:
            self._process_pkg(pkg)


class GromacsCompileAndRunBaseCheck(rfm.RegressionTest):
    def __init__(self, pkg, inputfile):
        self.pkg = pkg
        self.spec_list = pkg.spec.split(',')
        spec = ','.join([x for x in self.spec_list if not x.startswith('~') and  x != 'cpu'])
        self.inputfile = inputfile

        self.filename, _ = pkg.to_easyconfig()

        self.runtype = 'gpu'
        if 'cuda' in self.spec_list:
            if 'cpu' in self.spec_list:
                self.runtype = 'cpu'
        else:
            self.runtype = 'cpu'

        self.sourcesdir = 'src/eb_templates'
        self.sourcepath = 'deleteme.c'
        self.executable = 'deleteme'
        self.build_system = 'SingleSource'

        self.prebuild_cmds = [
            'echo "int main(){return 0;}" > deleteme.c'
        ]


        self.name = f'{pkg.name}_{pkg.version}_{pkg.toolchain}_{pkg.toolchain_version}_{util.toalphanum(spec)}_{self.runtype.upper()}Check'
        self.descr = f'{self.runtype.upper()} check for {pkg.name}_{pkg.version}_{pkg.toolchain}_{pkg.toolchain_version}_{util.toalphanum(spec)}'

        # Reset sources dir relative to the SCS apps prefix
        self.sourcesdir = os.path.join(self.current_system.resourcesdir,
                                       'Gromacs', self.inputfile)

        output_file = 'md.log'
        self.keep_files = [output_file]

        # TODO: add support for other input files here
        energy = sn.extractsingle(r'\s+Potential\s+Kinetic En\.\s+Total Energy'
                                  r'\s+Conserved En\.\s+Temperature\n'
                                  r'(\s+\S+){2}\s+(?P<energy>\S+)(\s+\S+){2}\n'
                                  r'\s+Pressure \(bar\)\s+Constr\. rmsd',
                                  output_file, 'energy', float, item=-1)
        energy_reference = -3270799.9

        self.sanity_patterns = sn.all([
            sn.assert_not_found('FAILED: Installation ended unsuccessfully', self.stdout),
            sn.assert_not_found('ERROR', self.stderr),
            sn.assert_found('COMPLETED: Installation ended successfully', self.stdout),
            sn.assert_found('Build succeeded', self.stdout),
            sn.assert_found('Finished mdrun', output_file),
            sn.assert_reference(energy, energy_reference, -0.001, 0.001)
        ])

        self.perf_patterns = {
            'perf': sn.extractsingle(r'Performance:\s+(?P<perf>\S+)',
                                     output_file, 'perf', float)
        }

        self.strict_check = False
        self.use_multithreading = False
        self.extra_resources = {
            'switches': {
                'num_switches': 1
            }
        }

        self.tags = {'gromacs', 'easybuild', 'external-resources'}

        self.valid_systems = ['daint:gpu','dom:gpu']
        self.valid_prog_environs = ['builtin']

        # generic single node job
        self.num_tasks = 1
        self.num_tasks_per_node = 1
        self.num_cpus_per_task = os.cpu_count()

        # if 'cuda' in self.variant:
        #     self.variables['CRAY_CUDA_MPS'] = '1'
        #     self.num_gpus_per_node = 1

        self.reference = {
            'dom:gpu': {'perf': (40.0, -0.05, None, 'ns/day')},
            'daint:gpu': {'perf': (38.8, -0.10, None, 'ns/day')}
        }

        # self.dep_name = re.sub(r'GromacsCheck', r'CompileGROMACSTest', self.name)
        # self.depends_on(self.dep_name, when=parent_part_env('login', 'builtin'))

        self.maintainers = ['VH']
        self.tags = {'gromacs', 'easybuild'}

    @rfm.run_after('setup')
    def config_eb(self):
        install_path = self.outputdir

        self.variables['EASYBUILD_PREFIX'] = f'{install_path}'
        self.variables['EB_CUSTOM_REPOSITORY'] = f'{self.stagedir}/production/easybuild'
        self.variables['EASYBUILD_TMPDIR'] = f'{self.stagedir}/tmpdir'
        self.variables['EASYBUILD_BUILDPATH'] = f'{self.stagedir}/builddir'
        self.variables['EASYBUILD_OPTARCH'] = r'${CRAY_CPU_TARGET}'
        self.variables['EASYBUILD_RECURSIVE_MODULE_UNLOAD'] = '0'
        self.prerun_cmds = [
            f'eb {self.filename} -r',
            f'module load GROMACS'
        ]

    @rfm.run_after('setup')
    def set_num_tasks(self):
        if self.current_partition.fullname in ['daint:gpu', 'dom:gpu']:
            self.num_tasks = 72
            self.num_tasks_per_node = 12
            self.num_cpus_per_task = 1

    @rfm.run_after('compile')
    def setup_gromacs(self):
        target = self.getdep(self.dep_name, 'builtin')

        self.variables.update(target.variables)
        self.modules += target.modules

        self.executable = 'gmx_mpi' if 'mpi' in self.spec_list else 'gmx'
        nb_type = 'gpu' if 'cuda' in self.spec_list else 'cpu'
        self.executable_opts = ['mdrun', '-dlb yes',
                                f'-ntomp {self.num_cpus_per_task}', '-npme -1',
                                f'-nb {nb_type}', '-s herflat.tpr']


def generate_gromacs_variants(inputfiles):
    gromacs_pkgs = []
    for spec in TEMPLATE_SPECS:
        if spec[0].startswith('GROMACS'):
            gromacs_pkgs += [Pkg(spec)]

    ret = []
    for infile in inputfiles:
        for pkg in gromacs_pkgs:
            for variant in pkg.get_variants():
                pkg = Pkg.from_pkg_info(variant)
                # duplicating GROMACS test when cuda in spec
                # this allows to test the GPU, but also the CPU of the build
                if 'cuda' in pkg.spec.split(','):
                    pkg.spec += ',cpu'
                    ret.append([pkg, infile])
                ret.append([Pkg.from_pkg_info(variant), infile])

    return ret

GROMACS_PKGS = generate_gromacs_variants(['herflat'])


@rfm.parameterized_test(*GROMACS_PKGS)
class GromacsCompileAndRunCheck(GromacsCompileAndRunBaseCheck):
    def __init__(self, pkg, inputfile):
        super().__init__(pkg, inputfile)
        spec = ','.join([x for x in self.spec_list if not x.startswith('~') and  x != 'cpu'])

        self.valid_systems = ['daint:gpu','dom:gpu']
        self.valid_prog_environs = ['builtin']

        # generic single node job
        self.num_tasks = 1
        self.num_tasks_per_node = 1
        self.num_cpus_per_task = os.cpu_count()

        if 'cuda' in self.spec_list:
            self.variables['CRAY_CUDA_MPS'] = '1'
            self.num_gpus_per_node = 1

        # f'GenerateEBFiles_{pkg.name}_{pkg.version}_{pkg.toolchain}_{pkg.toolchain_version}'
        # self.dep_name = re.sub(r'GromacsCheck', r'CompileGROMACSTest', self.name)
        self.dep_name = f'GenerateEBFiles_{self.pkg.name}_{self.pkg.version}_{self.pkg.toolchain}_{self.pkg.toolchain_version}'
        # self.depends_on(self.dep_name, when=parent_part_env('login', 'builtin'))
        self.depends_on(self.dep_name, when=parent_part_env('login', 'builtin'))

        self.maintainers = ['VH']
        self.tags = {'gromacs', 'easybuild'}


# # TODO: PARAMETERISE BASED ON GROMACS CLI FLAGS and input files
# @rfm.parameterized_test(*EASYBUILD_PKGS)
# class GromacsIndependentCheck(GromacsBaseCheck):
#     def __init__(self, variant, gromacs_version, installation_type):
#         super().__init__('md.log')
#         self.variant = variant
#         self.gromacs_version = gromacs_version
#         self.install_type = installation_type

#         self.valid_systems = ['daint:gpu','dom:gpu']
#         self.valid_prog_environs = ['builtin']

#         # generic single node job
#         self.num_tasks = 1
#         self.num_tasks_per_node = 1
#         self.num_cpus_per_task = os.cpu_count()

#         self.descr = 'GROMACS RunTime check'

#         self.reference = {
#             'dom:gpu': {'perf': (40.0, -0.05, None, 'ns/day')},
#             'daint:gpu': {'perf': (38.8, -0.10, None, 'ns/day')}
#         }

#     @rfm.run_after('setup')
#     def set_num_tasks(self):
#         if self.current_partition.fullname in ['daint:gpu', 'dom:gpu']:
#             self.num_tasks = 72
#             self.num_tasks_per_node = 12
#             self.num_cpus_per_task = 1
#             if 'cuda' in self.variant:
#                 self.variables['CRAY_CUDA_MPS'] = '1'
#                 self.num_gpus_per_node = 1
#                 self.modules = ['cudatoolkit']

#     # This has to be defined after set_num_tasks because it depends on self.num_cpus_per_task
#     @rfm.run_after('setup')
#     def config_gromacs(self):
#         install_path = self.outputdir

#         self.variables['GROMACS_ROOT'] = f'{install_path}'
#         self.variables['PATH'] = r'$GROMACS_ROOT/bin:$PATH'

#         self.executable = 'gmx_mpi' if 'mpi' in self.variant else 'gmx'
#         nb_type = 'gpu' if 'cuda' in self.variant else 'cpu'
#         self.executable_opts = ['mdrun', '-dlb yes',
#                                 f'-ntomp {self.num_cpus_per_task}', '-npme -1',
#                                 f'-nb {nb_type}', f'-bonded {nb_type}', '-s herflat.tpr']
