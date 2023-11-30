# Copyright 2016-2023 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import os
import stat

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
class munge_folder_permissions_check(rfm.RunOnlyRegressionTest):
    '''Check if the munged folders have right permissions'''

    #: Parameter pack encoding the munged folders and its permissions.
    #:
    #: The first element of the tuple refers to the folder name,
    #: the second is the desired permission of that folder.
    #:
    #: :type: `Tuple[str, str]`
    #: :values:
    munge_folders = parameter([('/etc/munge', '0700'),
                                ('/var/lib/munge', '0701'),
                                ('/var/log/munge', '0700'),
                                ('/run/munge', '0755')])

    #: munge user and group
    #:
    #: :type: :class:`str`
    #: :default: ``munge``
    munge_user = variable(str, value='munge')

    executable = 'getfacl'
    tags = {'system', 'anssi', 'munge'}

    @run_before('run')
    def set_executable_opts(self):
        self.executable_opts += [self.munge_folders[0]]

    @run_before('sanity')
    def skip_if_munge_not_installed(self):
        stderr = os.path.join(self.stagedir, sn.evaluate(self.stderr))
        try:
            sn.evaluate(sn.assert_not_found('No such file or directory',
                                            self.stderr))
        except SanityError:
            self.skip('munge is not installed')

    @deferrable
    def assert_option(self, regex):
        return sn.extractsingle(regex, self.stdout, 'value')

    @deferrable
    def assert_owner(self):
        return self.assert_option(r'^#[\s]*owner:[\s]*(?P<value>.*)')

    @deferrable
    def assert_group(self):
        return self.assert_option(r'^#[\s]*group:[\s]*(?P<value>.*)')

    @deferrable
    def assert_permission_user(self):
        return self.assert_option(r'^user::[\s]*(?P<value>.*)')

    @deferrable
    def assert_permission_group(self):
        return self.assert_option(r'^group::[\s]*(?P<value>.*)')

    @deferrable
    def assert_permission_other(self):
        return self.assert_option(r'^other::[\s]*(?P<value>.*)')

    @sanity_function
    def assert_checks(self):
        value = ('?' +
                 self.assert_permission_user() +
                 self.assert_permission_group() +
                 self.assert_permission_other)

        folder_permission = stat.filemode(int(self.munge_folders[1], 8))

        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_eq(self.assert_owner(), self.munge_user,
                         msg=f'folder {self.munge_folders[0]} is not owned by '
                             f'{self.munge_user}'),
            sn.assert_eq(self.assert_group(), self.munge_user,
                         msg=f'folder {self.munge_folders[0]} is not in group '
                             f'{self.munge_user}'),
            sn.assert_eq(value, folder_permission,
                         msg=f'File permissions {value} are not correct '
                             f'{folder_permission}'),
            ])
