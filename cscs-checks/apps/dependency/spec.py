import contextlib
import itertools
import os
import re
import fnmatch
from datetime import datetime

TEMPLATE_SPECS = {
    ('GROMACS@2020.4', 'CrayGNU@20.08') : {
        'spec' : 'mpi,openmp,cuda,ownfftw,~shared,~plumed,~pat',
        'builddependencies' : {
            ('CMake@3.14.5', 'system'): {}
        },
        'dependencies' : {
            ('zlib@1.2.11', 'CrayGNU@20.08'): {},
            ('GSL@2.5', 'CrayGNU@20.08'): {}
        },
        'variants' : ['pat', '~cuda,pat']
    },
    ('GROMACS@2020.3', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2020.4', 'CrayGNU@20.08')
    },
    ('GROMACS@2020.2', 'CrayGNU@20.08') : {
        'from' : ('GROMACS@2020.4', 'CrayGNU@20.08'),
        'dependencies' : {
            ('PLUMED@2.6.1', 'CrayGNU@20.08') : {
                'when' : 'plumed'
            }
        },
        'variants' : ['pat', '~cuda,pat', 'plumed,pat', '~cuda,plumed,pat']
    },
    ('PLUMED@2.6.1', 'CrayGNU@20.08') : {
        'spec' : 'mpi,openmp',
        'dependencies' : {
            ('zlib@1.2.11', 'CrayGNU@20.08'): {},
            ('GSL@2.5', 'CrayGNU@20.08'): {}
        },
    },
    ('GSL@2.5', 'CrayGNU@20.08') : {
        'spec' : 'openmp,optarch,unroll',
    },
    ('zlib@1.2.11', 'CrayGNU@20.08') : {
        'spec' : 'shared',
    },
}

class Pkg():
    # def __init__(self, name_version, toolchain_version):
    def __init__(self, name_and_toolchain, template_specs=None):
        if template_specs:
            self.template_specs = template_specs
        else:
            self.template_specs = TEMPLATE_SPECS
        # self.name_and_version = name_version
        # self.toolchain_and_version = toolchain_version
        # self.pkg_info = self._extract_pkg_info((self.name_and_version, self.toolchain_and_version))
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
        # import pprint
        ret  = f'Pkg({self.name_and_toolchain}):\n'
        # ret += f'             pkg_info: {self.pkg_info}\n'
        ret += f'                 spec: {self.spec}\n'
        ret += f'         dependencies: {self.dependencies}\n'
        ret += f'    builddependencies: {self.builddependencies}\n'
        ret += f'             variants: {self.variants}\n'
        ret += f'        template_file: {self.template_file}'
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
        base = copy.deepcopy(self.pkg_info)
        ret = [base]
        if 'variants' in base[spec] and 'spec' in base[spec]:
            for variant in base[spec]['variants']:
                new_base = copy.deepcopy(base)
                new_base[spec]['spec'] = self._replace_spec_with_variant(new_base[spec]['spec'], variant)
                ret += [new_base]

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

    def _remove_dependency(self, name, dependencies, spec_list):
        import copy
        deps = copy.deepcopy(dependencies)
        dep_list = self._find_all_dependencies_with_name(name, deps)
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
        if 'cuda' in spec_list:
            version_suffix += '-cuda'

        if 'pat' in spec_list:
            version_suffix += '-pat'

        builddeps = None
        deps = None
        if 'plumed' in spec_list:
            if builddependencies:
                builddeps = self._remove_dependency('PLUMED', builddependencies, spec_list)
                plumed_build = self._find_all_dependencies_with_name('PLUMED', builddeps)

            if dependencies:
                deps = self._remove_dependency('PLUMED', dependencies, spec_list)
                plumed = self._find_all_dependencies_with_name('PLUMED', deps)

            if plumed_build or plumed:
                version_suffix += '-plumed'

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

        filename = f"{ec['name']}-{ec['version']}-{ec['toolchain']}-{ec['toolchain_version']}{ec['version_suffix']}.eb"

        return filename, ec


            # with open(os.path.join(easybuild_dir, filename), 'w') as fp:
            #     fp.write(template.render(ec))
            #     fp.write('\n')

        # for variant in self.pgk_list:
        #     gromacs_filename, ec = gromacs_variant_to_dict(variant, self.current_system.name)

        #     if gromacs_filename:
        #         easybuild_dir = os.path.join(self.outputdir, 'easybuild','easyconfigs','g','GROMACS')
        #         mkdirp(easybuild_dir)

        #         with open(os.path.join(easybuild_dir, gromacs_filename), 'w') as fp:
        #             fp.write(template.render(ec))
        #             fp.write('\n')

GROMACS_SPECS = [
    # 2020
    Pkg(('GROMACS@2020.4', 'CrayGNU@20.08')),
    Pkg(('GROMACS@2020.3', 'CrayGNU@20.08')),
    Pkg(('GROMACS@2020.2', 'CrayGNU@20.08')),
    Pkg(('GROMACS@2020.1', 'CrayGNU@20.08')),
    Pkg(('GROMACS@2020',   'CrayGNU@20.08')),
    Pkg(('GROMACS@2018',   'CrayGNU@20.08'))
]

def process_pkg(pkg):
    for spec in pkg.get_build_tree():
        template_pkg = Pkg.from_pkg_info(spec)
        for variant in template_pkg.get_variants():
            pkg_variant = Pkg.from_pkg_info(variant)
            filename, ec = pkg_variant.to_easyconfig()

            # with open(os.path.join(easybuild_dir, filename), 'w') as fp:
            #     fp.write(template.render(ec))
            #     fp.write('\n')

        # for variant in self.pgk_list:
        #     gromacs_filename, ec = gromacs_variant_to_dict(variant, self.current_system.name)

        #     if gromacs_filename:
        #         easybuild_dir = os.path.join(self.outputdir, 'easybuild','easyconfigs','g','GROMACS')
        #         mkdirp(easybuild_dir)

        #         with open(os.path.join(easybuild_dir, gromacs_filename), 'w') as fp:
        #             fp.write(template.render(ec))
        #             fp.write('\n')

    print(f'pkg.name: {pkg.name}')


def generate_ebfile():
    pgk_list = GROMACS_SPECS
    for pkg in pgk_list:
        process_pkg(pkg)


from math import floor, sqrt
try:
    long
except NameError:
    long = int

def fac(n):
    step = lambda x: 1 + (x<<2) - ((x>>1)<<1)
    maxq = long(floor(sqrt(n)))
    d = 1
    q = 2 if n % 2 == 0 else 3
    while q <= maxq and n % q != 0:
        q = step(d)
        d += 1
    return [q] + fac(n // q) if q <= maxq else [n]

def difference(list1, list2):
    # return list(list(set(list1)-set(list2)) + list(set(list2)-set(list1)))
    li_dif = [i for i in list1 + list2 if i not in list1 or i not in list2]
    return li_dif

if __name__ == '__main__':
    # spec = Specs('GROMACS', '2020.2', 'CrayGNU', '20.08')
    # print('original spec: ', spec)
    # for var in spec.get_variants():
    #     print('modified spec: ', var)

    # for spec in TEMPLATE_SPECS:
    #     ret = TEMPLATE_SPECS[spec]
    #     print(f'ret.name: {ret}')

    # for tspec in TEMPLATE_SPECS:
    # pkg = Pkg(('GROMACS@2020.2', 'CrayGNU@20.08'))
    # print(pkg)

    # pkg = Pkg(('CMake@3.14.5', 'system'))
    # print(pkg)

    # generate_ebfile()
    print(fac(16))
    print(fac(24))
    print(fac(15*30))
    print(fac(17))

    prod = 1
    ret = []
    primes = fac(15*30*13)
    for x in primes:
        if (prod * x <= 12):
            prod *= x
            ret += [x]

    print('ret: ', ret)
    print('primes: ', primes)
    print('diff: ', difference(primes, ret))



    # for s in spec.get_build_tree():
    #     print('\n',s)
    #     # ret = TEMPLATE_SPECS[spec]
        # print(f'ret.name: {ret}')
