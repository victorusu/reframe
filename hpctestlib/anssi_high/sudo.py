# Copyright 2016-2023 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import os

import reframe as rfm
import reframe.utility.sanity as sn
import reframe.utility.typecheck as typ

from reframe.core.exceptions import SanityError


class SkipIfNotRoot(rfm.RegressionMixin):
    @run_after('init')
    def skip_if_not_root(self):
        self.skip_if(os.getuid() != 0,
                     msg='Skipping test because it has be executed as root')


class SkipIfNotLocal(rfm.RegressionMixin):
    @run_before('run')
    def skip_if_not_local(self):
        self.skip_if(not self.is_local,
                     msg="Skipping the test because it is not local")


@rfm.simple_test
class sudo_config_check(rfm.RunOnlyRegressionTest, SkipIfNotRoot):
    '''Ensure sudo is functional'''

    #: sudoers file
    #:
    #: :type: :class:`str`
    #: :default: ``/etc/sudoers``
    sudoers_file = variable(str, value='/etc/sudoers')

    executable = '/usr/sbin/visudo'
    executable_opts = ['-cf']
    tags = {'system', 'anssi', 'sudo'}

    @run_before('run')
    def set_sudoers_file(self):
        self.executable_opts += [self.sudoers_file]

    @sanity_function
    def assert_checks(self):

        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_found(rf'{self.sudoers_file}:\s+parsed OK', self.stdout)]
            )


@rfm.simple_test
class sudoers_config_check(rfm.RunOnlyRegressionTest,
                           SkipIfNotRoot, SkipIfNotLocal):
    '''Ensure sudo has the right configurations'''

    #: sudoers file
    #:
    #: :type: :class:`str`
    #: :default: ``/etc/sudoers``
    sudoers_file = variable(str, value='/etc/sudoers')

    #: eBPF program tags
    #:
    #: :type: `List[str]`. The list should contain all the whitelist of
    #: ebpf program tags
    sudo_defaults = variable(typ.List[str])

    executable = 'echo done'
    tags = {'system', 'anssi', 'sudo'}

    def assert_all_defaults(self):
        unset_defaults = set()
        for defaults in self.sudo_defaults:
            if defaults == '':
                continue
            try:
                sn.evaluate(sn.assert_found(rf'^[\s]*Defaults.*\b'
                                            rf'{defaults}\b.*$',
                                            self.sudoers_file))
            except SanityError:
                unset_defaults.add(defaults)

        return sn.assert_eq(len(unset_defaults), 0,
                            msg=f'defaults {unset_defaults} are not configured '
                                 'in the system')

    @sanity_function
    def assert_checks(self):
        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            self.assert_all_defaults(),
            sn.assert_not_found(rf'(^(?!#).*[\s]+\!authenticate.*$)',
                                self.sudoers_file,
                                msg='!authenticated found in file '
                                    f'{self.sudoers_file}'),
            sn.assert_not_found(rf'(^(?!#).*[\s]+NOPASSWD[\s]*\:.*$)',
                                self.sudoers_file,
                                msg=f'NOPASSWD found in file {self.sudoers_file}')
            ])


@rfm.simple_test
class sudoersd_config_check(rfm.RunOnlyRegressionTest,
                            SkipIfNotRoot, SkipIfNotLocal):
    '''Ensure that additional sudoers configs have the right configurations'''

    #: sudoers.d directory
    #:
    #: :type: :class:`str`
    #: :default: ``/etc/sudoers``
    sudoersd = variable(str, value='/etc/sudoers.d')

    executable = 'echo done'
    tags = {'system', 'anssi', 'sudo'}

    def get_all_files(self):
        files = []
        for (dirpath, _, filenames) in os.walk(self.sudoersd):
            for f in filenames:
                files.append(os.path.join(dirpath, f))

        return files

    @sanity_function
    def assert_checks(self):
        files = self.get_all_files()

        return sn.all(sn.chain(
            sn.map(lambda x: sn.assert_not_found(
                rf'(^(?!#).*[\s]+\!authenticate.*$)', x,
                msg=f'!authenticated found in file {x}'), files),
            sn.map(lambda x: sn.assert_not_found(
                rf'(^(?!#).*[\s]+NOPASSWD[\s]*\:.*$)', x,
                msg=f'NOPASSWD found in file {x}'), files),
        ))
