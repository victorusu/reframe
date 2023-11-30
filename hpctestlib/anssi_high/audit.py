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
class audit_rules_persistency_check(rfm.RunOnlyRegressionTest,
                                    SkipIfNotRoot, SkipIfNotLocal):
    '''Check if the audit rules cannot be changed'''

    executable = 'echo'
    executable_opts = ['done']
    audit_rules_file = '/etc/audit/rules.d/audit.rules'
    tags = {'system', 'anssi', 'audit'}

    @sanity_function
    def assert_checks(self):
        value = sn.extractall(r'^-e[\s]+(?P<value>\d+)', self.audit_rules_file,
                              'value', int)

        return sn.all([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_true(value, msg='-e 2 not found in the '
                                      f'{self.audit_rules_file} file'),
            sn.assert_eq(sn.count(value), 1,
                         msg=f'There are multiple entries to -e flag'),
            sn.assert_eq(value[0], 2,
                         msg=f'The -e flag should be set to 2'),
            ])


@rfm.simple_test
class audit_immutable_rules_check(audit_rules_persistency_check):
    audit_rules_file = '/etc/audit/rules.d/immutable.rules'


@rfm.simple_test
class audit_monitor_path_rules_check(rfm.RunOnlyRegressionTest,
                                    SkipIfNotRoot, SkipIfNotLocal):
    '''Check if the audit rules cannot be changed'''

    #: Audit path rules to be monitored
    #:
    #: :type: `Dict[str, str]`. The key should be the path to a file or folder.
    #:   and the value should be a string with the audit rule options.
    #:   E.g., '{"/etc/selinux": "key=MAC-policy,perm=rwa"}'
    audit_path_rules = variable(typ.Dict[str,str])

    executable = 'echo'
    executable_opts = ['done']
    audit_rules_file = '/etc/audit/rules.d'
    tags = {'system', 'anssi', 'audit'}

    def get_all_files(self):
        files = []
        for (dirpath, _, filenames) in os.walk(self.audit_rules_file):
            for f in filenames:
                files.append(os.path.join(dirpath, f))

        return files

    def explode_rules_str(self, opts_str):
        result = {}
        for opt in opts_str.split(','):
            if opt == '':
                continue

            opt_parts = opt.split('=', maxsplit=2)
            keystr = opt_parts[0]
            valstr = opt_parts[1] if len(opt_parts) > 1 else ''
            result[keystr] = valstr

        return result

    def is_key_ne_item_dict(self, key, item, dictionary):
        if key and key in dictionary:
            if item.group(key) == dictionary[key]:
                return False
        return True

    def incr_counter(self, key, dictionary):
        if key:
            if key in dictionary:
                dictionary[key] += 1
            else:
                dictionary[key] = 1

    def process_rules(self, found_rules, error_msgs, monitored_files,
                            tobe_monitored_files):
        for found_rule in found_rules:
            self.incr_counter(found_rule.group('path'), monitored_files)

            for audit_item, audit_rules in self.audit_path_rules.items():
                rules = self.explode_rules_str(audit_rules)

                if audit_item == found_rule.group('path'):
                    self.incr_counter(found_rule.group('path'),
                                      tobe_monitored_files)
                    if self.is_key_ne_item_dict('key', found_rule, rules):
                        error_msgs.add('audit keyname not correct for path '
                                       f'{audit_item}')
                    if self.is_key_ne_item_dict('perm', found_rule, rules):
                        error_msgs.add('audit permission not correct for path '
                                       f'{audit_item}')


    def evaluate_sanity(self, error_msgs, monitored_files, tobe_monitored_files):
        result = []
        nl = '\n'

        all_sanities = sn.chain([
            sn.assert_not_found('command not found', self.stderr),
            sn.assert_not_found('Permission denied', self.stderr),
            sn.assert_ne(sn.count(error_msgs), 0,
                         msg=f'{self.audit_path_rules} is not being monitored '
                             'via audit'),
            sn.assert_eq(sn.count(error_msgs), 0,
                         msg=f'{nl.join(error_msgs)}')],
            sn.map(lambda x: sn.assert_in(x, tobe_monitored_files,
                                          msg=f'File {x} is not being '
                                              'monitored'),
                                              self.audit_path_rules),
            sn.map(lambda x: sn.assert_eq(x[1], 1,
                                          msg=f'File {x[0]} is monitored {x[1]}'
                                              ' times'),
                                              monitored_files.items()),
            )

        for sanity in all_sanities:
            try:
                sn.evaluate(sanity)
            except SanityError as e:
                result.append(str(e))

        return result


    @sanity_function
    def assert_checks(self):
        error_msgs = set()
        monitored_files = {}
        tobe_monitored_files = {}
        files = self.get_all_files()

        for file in files:
            found_rules = sn.evaluate(
                sn.findall(
                    r'^(?!#)'
                    r'(?=.*(?:-F key=|-k\s+)(?P<key>[\S+]*))?'
                    r'(?=.*(?:-F perm=|-p\s+)(?P<perm>[\S+]*))?'
                    r'(?=.*(?:-F dir=|-F path=|-w\s+)(?P<path>[\S+]*))?'
                    r'.*$',file))

            self.process_rules(found_rules, error_msgs,
                               monitored_files, tobe_monitored_files)

        result = self.evaluate_sanity(error_msgs, monitored_files, tobe_monitored_files)
        return sn.assert_false(result, msg='\n'+'\n'.join(result))
