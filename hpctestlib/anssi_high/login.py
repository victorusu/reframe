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
class authselect_check(rfm.RunOnlyRegressionTest):
    '''Check if authselect is correctly configured'''

    executable = 'authselect'
    executable_opts = ['check']
    tags = {'system', 'anssi', 'authselect'}

    @sanity_function
    def assert_checks(self):
        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_found('Current configuration is valid.', self.stdout)
            ])


@rfm.simple_test
class authselect_profile_check(rfm.RunOnlyRegressionTest):
    '''Check if authselect has the right profile'''

    #: Profile to be checked via authselect
    #:
    #: :type: :class:`str`
    #: :default: ``None``
    authselect_profile = variable(str)

    executable = 'authselect'
    executable_opts = ['current']
    tags = {'system', 'anssi', 'authselect'}

    @sanity_function
    def assert_checks(self):
        cur = sn.extractsingle(r'Profile ID:[\s]*(?P<profile>\S*)',
                               self.stdout, 'profile')

        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_eq(cur, self.authselect_profile),
            ])


@rfm.simple_test
class authselect_features_check(rfm.RunOnlyRegressionTest):
    '''Check if authselect has the right features enabled'''

    #: List of features to enable in authselect
    #:
    #: :type: `List[str]`. The list should contain all the features enabled in
    #: authselect
    authselect_features = variable(typ.List[str])

    executable = 'authselect'
    executable_opts = ['current']
    tags = {'system', 'anssi', 'authselect'}

    def assert_all_features(self):
        unset_features = set()
        for feature in self.authselect_features:
            if feature == '':
                continue
            try:
                sn.evaluate(sn.assert_found(rf'-[\s]*{feature}$', self.stdout))
            except SanityError:
                unset_features.add(feature)

        return sn.assert_eq(len(unset_features), 0,
                            msg=f'Features {unset_features} are not configured '
                                 'in the system')

    @sanity_function
    def assert_checks(self):
        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_not_found('error', self.stderr),
            sn.assert_not_found('Unable to get profile information', self.stderr),
            self.assert_all_features()
            ])


@rfm.simple_test
class pam_login_check(rfm.RunOnlyRegressionTest):
    '''Check if login has the right configuration'''

    executable = 'cat'
    executable_opts = ['/etc/pam.d/login']
    tags = {'system', 'anssi', 'pam'}

    def assert_namespace(self):
        return sn.assert_found('session[\s]*required[\s]*pam_namespace.so',
                               self.stdout),

    def assert_selinux(self):
        return sn.assert_found(r'session[\s]*required[\s]*pam_selinux.so[\s]*'
                               r'close.*(\n.*)+session[\s]*required[\s]*'
                               r'pam_loginuid.so(\n.*)+session[\s]*required'
                               r'[\s]*pam_selinux.so[\s]*open', self.stdout),

    def assert_loginuid(self):
        cur = sn.extractall(r'session[\s]*required[\s]*pam_loginuid.so',
                            self.stdout)

        return sn.assert_eq(sn.count(cur), 1),

    def assert_selinux_close_first_session(self):
        return sn.assert_not_found(r'^session[\s]*(required|included)*[\s]*.*'
                                   r'(\n.*)+session[\s]*required[\s]*'
                                   r'pam_selinux.so[\s]*close.*', self.stdout),

    @sanity_function
    def assert_checks(self):
        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            self.assert_namespace,
            self.assert_selinux,
            self.assert_loginuid,
            self.assert_selinux_close_first_session,
            ])


@rfm.simple_test
class pam_system_check(rfm.RunOnlyRegressionTest):
    '''Check if system-auth has the right configurations'''

    executable = 'cat'
    executable_opts = ['/etc/pam.d/system-auth']
    tags = {'system', 'anssi', 'pam'}

    def assert_pam_env(self):
        return sn.assert_found('^auth[\s]*required[\s]*pam_env.so',
                               self.stdout),

    def assert_limits(self):
        return sn.assert_found(r'session[\s]*required[\s]*pam_limits.so',
                               self.stdout),

    def assert_unix_account(self):
        cur = sn.extractall(r'account[\s]*required[\s]*pam_unix.so',
                            self.stdout)

        return sn.assert_eq(sn.count(cur), 1),

    def assert_unix(self):
        cur = sn.extractall(r'session[\s]*required[\s]*pam_unix.so',
                            self.stdout)

        return sn.assert_eq(sn.count(cur), 1),

    @sanity_function
    def assert_checks(self):
        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            self.assert_pam_env,
            self.assert_limits,
            self.assert_unix_account,
            self.assert_unix,
            ])


@rfm.simple_test
class login_defs_check(rfm.RunOnlyRegressionTest, SkipIfNotLocal):
    '''Check if the min password len is properly set'''

    #: Umask used by login on non-PAM enabled systems
    #:
    #: :type: :class:`str`
    #: :default: ``0700``
    login_defs_umask = variable(str, value='0700')

    #: Umask used to create home folders
    #:
    #: :type: :class:`str`
    #: :default: ``0700``
    login_defs_home_mode = variable(str, value='0700')

    #: Encrypt method used to encrypting passwords
    #:
    #: :type: :class:`str`
    #: :default: ``0700``
    login_defs_encrypt_method = variable(str, value='SHA512')

    #: Define a default home when one can't cd to the home directory?
    #:
    #: :type: :class:`str`
    #: :default: ``no``
    login_defs_default_home = variable(str, value='no')

    #: Allows userdel to remove user groups if no members exist
    #:
    #: :type: :class:`str`
    #: :default: ``yes``
    login_defs_usergrps_enabled = variable(str, value='yes')

    #: Create user home directories by default when running useradd?
    #:
    #: :type: :class:`str`
    #: :default: ``yes``
    login_defs_create_home = variable(str, value='yes')

    executable = 'echo done'
    tags = {'system', 'anssi', 'login'}

    @run_after('init')
    def set_login_defs_file(self):
        self.login_defs_file = '/etc/login.defs'

    def assert_configuration(self, name, regex, ref):
        found = sn.extractall(regex,
                              self.login_defs_file,
                              'value')

        return (sn.assert_true(found, msg=f'{name} was not found in the '
                                           'configuration file '
                                          f'{self.login_defs_file}')
                and
                sn.assert_eq(found[-1], ref,
                             msg=f'Defined {name} ({found[-1]}) is different '
                                 f'from the expected value ({ref})'))


    @sanity_function
    def assert_checks(self):
        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            self.assert_configuration(name='UMASK',
                                      regex=r'^UMASK[\s]*(?P<value>[0-9]*)',
                                      ref=self.login_defs_umask),
            self.assert_configuration(name='HOME_MODE',
                                      regex=r'^HOME_MODE[\s]*(?P<value>[0-9]*)',
                                      ref=self.login_defs_home_mode),
            self.assert_configuration(name='ENCRYPT_METHOD',
                                      regex=r'^ENCRYPT_METHOD[\s]*(?P<value>\S*)',
                                      ref=self.login_defs_encrypt_method),
            self.assert_configuration(name='DEFAULT_HOME',
                                      regex=r'^DEFAULT_HOME[\s]*(?P<value>\S*)',
                                      ref=self.login_defs_default_home),
            self.assert_configuration(name='USERGROUPS_ENAB',
                                      regex=r'^USERGROUPS_ENAB[\s]*(?P<value>\S*)',
                                      ref=self.login_defs_usergrps_enabled),
            self.assert_configuration(name='CREATE_HOME',
                                      regex=r'^CREATE_HOME[\s]*(?P<value>\S*)',
                                      ref=self.login_defs_create_home),
            ])
