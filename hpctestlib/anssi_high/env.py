# Copyright 2016-2023 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import os

import reframe as rfm
import reframe.utility.sanity as sn
import reframe.utility.osext as osext


class SkipIfNotRoot(rfm.RegressionMixin):
    @run_after('init')
    def skip_if_not_root(self):
        self.skip_if(os.getuid() != 0,
                     msg='Skipping test because it has be executed as root')


class SkipIfRoot(rfm.RegressionMixin):
    @run_after('init')
    def skip_if_root(self):
        self.skip_if(os.getuid() == 0,
                     msg='Skipping test because it has been executed as root')


class SkipIfNotLocal(rfm.RegressionMixin):
    @run_before('run')
    def skip_if_not_local(self):
        self.skip_if(not self.is_local,
                     msg="Skipping the test because it is not local")


@rfm.simple_test
class tmpout_check(rfm.RunOnlyRegressionTest,
                   SkipIfRoot, SkipIfNotLocal):
    '''Check if the TMOUT variable is correctly configured'''

    #: TMOUT value
    #:
    #: :type: :class:`int`
    #: :default: ``3600``
    tmout_value = variable(int, value=3600)

    executable = 'echo'
    executable_opts = ['TMOUT=$TMOUT']
    tags = {'system', 'anssi', 'env'}

    @run_after('init')
    def get_readonly_tmout(self):
        self.readonly = osext.run_command('bash -ic "readonly -p"')

    @sanity_function
    def assert_checks(self):
        value = sn.extractsingle('^TMOUT=(?P<value>[0-9]+)$', self.stdout,
                              'value', int)
        tmout = sn.extractsingle_s(r'declare[\s]*-rx[\s]*TMOUT="'
                                   r'(?P<value>[0-9]*)"', self.readonly.stdout,
                                   'value', int)

        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_eq(value, tmout,
                         msg=f'TMOUT is not set to be readonly'),
            sn.assert_eq(value, self.tmout_value,
                         msg=f'TMOUT value ({value}) is different from the '
                             f'reference ({self.tmout_value})'),
            ])


@rfm.simple_test
class umask_bashrc_check(rfm.RunOnlyRegressionTest):
    '''Check if the umask is correctly configured'''

    #: Umask used by login on non-PAM enabled systems
    #:
    #: :type: :class:`str`
    #: :default: ``0700``
    login_defs_umask = variable(str)

    executable = 'cat'
    executable_opts = ['/etc/bashrc']
    skip_if_empty = False
    tags = {'system', 'anssi', 'env'}

    @sanity_function
    def assert_checks(self):
        umasks = sn.extractall(r'^[\s]*umask[\s]+(?P<value>\d+)', self.stdout,
                              'value')

        files = ','.join(self.executable_opts)

        if self.skip_if_empty:
            self.skip_if(umasks == [], msg=f'Umask not defined in {files}')

        return sn.all(sn.chain([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_ne(sn.count(umasks), 0,
                         msg=f'umask not defined in {files}'),
            sn.assert_eq(sn.count(umasks), 1,
                         msg=f'umask defined multiple times in {files}')],
            sn.map(lambda x: sn.assert_eq(x, self.login_defs_umask,
                         msg=f'umask is not set to {self.login_defs_umask} '
                             f'in {files}'),
                         umasks),
            ))


@rfm.simple_test
class umask_profile_check(umask_bashrc_check):
    '''Check if the umask is correctly configured'''

    executable_opts = ['/etc/profile']


@rfm.simple_test
class umask_profiled_check(umask_bashrc_check):
    '''Check if the umask is correctly configured'''

    executable_opts = ['/etc/profile.d/*.sh', '/etc/profile.d/*.local']
    skip_if_empty = True
