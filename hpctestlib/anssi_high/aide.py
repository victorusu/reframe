# Copyright 2016-2023 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import reframe as rfm
import reframe.utility.sanity as sn


class SudoCmd(rfm.RegressionMixin):
    @run_before('run')
    def set_sudo_cmd(self):
        self.executable = f'sudo {self.executable}'


@rfm.simple_test
class aide_database_check(rfm.RunOnlyRegressionTest, SudoCmd):
    '''Check if aide is installed and is minimally configured'''

    executable = 'aide'
    executable_opts = ['--check']
    tags = {'system', 'anssi', 'aide'}

    @run_before('run')
    def set_sudo_cmd(self):
        self.executable = f'{self.executable}'

    @sanity_function
    def assert_checks(self):

        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_not_found(r'Couldn\'t open file .* for reading', self.stdout),
            sn.assert_found('Summary', self.stdout),
            sn.assert_found('aide.db.gz', self.stdout),
            sn.assert_found('End timestamp', self.stdout)]
            )


@rfm.simple_test
class aide_configured_check(rfm.RunOnlyRegressionTest, SudoCmd):
    '''Check if aide was configured with proper rule options'''

    #: Parameter listing the aide rules that must have
    #: the aide_rule_opts configured
    #:
    #: :type: :class:`str`
    #: :values: ``['NORMAL', 'DIR', 'PERMS', 'LOG',
    #              'CONTENT_EX', 'DATAONLY']``
    aide_rules = parameter(['NORMAL', 'DIR', 'PERMS', 'LOG',
                            'CONTENT_EX', 'DATAONLY'])

    #: Parameter listing the aide rules options that must have
    #: be configured to the
    #:
    #: :type: :class:`str`
    #: :values: ``['NORMAL', 'DIR', 'PERMS', 'LOG',
    #              'CONTENT_EX', 'DATAONLY']``
    aide_rules_opts = parameter(['ACL', 'SELinux', 'XATTRS'])

    executable = 'cat'
    executable_opts = ['/etc/aide.conf']
    tags = {'system', 'anssi', 'aide'}

    @sanity_function
    def assert_checks(self):

        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_not_found(r'Couldn\'t open file .* for reading',
                                self.stdout),
            sn.assert_found(rf'^{self.aide_rules}\s?=\s?.*'
                            rf'{self.aide_rules_opts.lower()}.*', self.stdout)]
            )

