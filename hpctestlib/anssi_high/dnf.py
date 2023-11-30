# Copyright 2016-2023 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import os

import reframe as rfm
import reframe.utility.sanity as sn


class SkipIfNotRoot(rfm.RegressionMixin):
    @run_after('init')
    def skip_if_not_root(self):
        self.skip_if(os.getuid() != 0,
                     msg='Skipping test because it has be executed as root')


@rfm.simple_test
class packages_updated_check(rfm.RunOnlyRegressionTest, SkipIfNotRoot):
    '''Ensure there are no security updates'''

    executable = 'dnf'
    executable_opts = ['update', '--security', '<<<N']
    tags = {'system', 'anssi', 'packages', 'dnf'}

    @sanity_function
    def assert_checks(self):
        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_found('Dependencies resolved.', self.stdout),
            sn.assert_found('Nothing to do.', self.stdout),
            sn.assert_found('Complete!', self.stdout),
            ])


@rfm.simple_test
class needs_reboot_check(rfm.RunOnlyRegressionTest, SkipIfNotRoot):
    '''Check if the system needs rebooting'''

    executable = 'needs-restarting'
    executable_opts = ['-r']
    tags = {'system', 'anssi', 'packages', 'dnf'}

    @sanity_function
    def assert_checks(self):
        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_found('Reboot should not be necessary.', self.stdout,
                            msg='System requires a reboot after update'),
            ])


@rfm.simple_test
class gpg_enabled_check(rfm.RunOnlyRegressionTest, SkipIfNotRoot):
    '''Check if all repos have GPG enabled'''

    executable = 'cat'
    executable_opts = ['/etc/yum.conf', '/etc/yum.repos.d/*']
    tags = {'system', 'anssi', 'packages', 'dnf'}

    @sanity_function
    def assert_checks(self):
        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_not_found(r'^gpgcheck[\s]*=[\s]*0', self.stdout,
                                msg='Found repo with gpgpcheck=0'),
            ])
