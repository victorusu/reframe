# Copyright 2016-2023 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import json
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


@rfm.simple_test
class ebpf_prog_list_check(rfm.RunOnlyRegressionTest, SkipIfNotRoot):
    '''Ensure ebpf is functional'''

    #: eBPF program tags
    #:
    #: :type: `List[str]`. The list should contain all the whitelist of
    #: ebpf program tags
    ebpf_tags = variable(typ.List[str])

    executable = 'bpftool'
    executable_opts = ['prog', 'list', '--json', '--pretty']
    tags = {'system', 'anssi', 'ebpf'}

    @run_before('sanity')
    def skip_if_command_not_found(self):
        stderr = os.path.join(self.stagedir, sn.evaluate(self.stderr))
        try:
            sn.evaluate(sn.assert_not_found('command not found', stderr))
        except SanityError as e:
            self.skip('bpftool is not installed')

    @run_before('sanity')
    def read_output_json(self):
        self.ebpf_json = None
        stdout = os.path.join(self.stagedir, sn.evaluate(self.stdout))
        with open(stdout, 'r') as fp:
            self.ebpf_json = json.load(fp)

    @sanity_function
    def assert_checks(self):
        return sn.all(sn.chain(
            [sn.assert_not_found('Permission denied', self.stderr)],
            sn.map(lambda x: sn.assert_in(x['tag'], self.ebpf_tags,
                                          msg='found unknown ebpf program '
                                              f'"{x["name"]}" with '
                                              f'id "{x["id"]}" and '
                                              f'tag "{x["tag"]}"'),
                                          self.ebpf_json),
            ))
