# Copyright 2016-2020 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

#
# Test case graph functionality
#

import collections
import itertools

import reframe as rfm
import reframe.utility as util
from reframe.core.exceptions import DependencyError


def build_deps(cases, default_cases=None):
    '''Build dependency graph from test cases.

    The graph is represented as an adjacency list in a Python dictionary
    holding test cases. The dependency information is also encoded inside each
    test case.
    '''

    # Index cases for quick access
    def build_index(cases):
        if cases is None:
            return {}

        ret = {}
        for c in cases:
            cname = c.check.name
            ret.setdefault(cname, [])
            ret[cname].append(c)

        return ret

    all_cases_map = build_index(cases)
    default_cases_map = build_index(default_cases)

    def resolve_dep(src, dst):
        errmsg = f'could not resolve dependency: {src!r} -> {dst!r}'
        try:
            ret = all_cases_map[dst]
        except KeyError:
            # try to resolve the dependency in the fallback map
            try:
                ret = default_cases_map[dst]
            except KeyError:
                raise DependencyError(errmsg) from None

        if not ret:
            raise DependencyError(errmsg)

        return ret

    # NOTE on variable names
    #
    # c stands for check or case depending on the context
    # p stands for partition
    # e stands for environment
    # t stands for target

    # We use an ordered dict here, because we need to keep the order of
    # partitions and environments
    graph = collections.OrderedDict()
    for c in cases:
        psrc = c.partition.name
        esrc = c.environ.name
        for dep in c.check.user_deps():
            tname, when = dep
            for d in resolve_dep(c, tname):
                pdst = d.partition.name
                edst = d.environ.name
                if when((psrc, esrc), (pdst, edst)):
                    c.deps.append(d)

        graph[c] = util.OrderedSet(c.deps)

    # Calculate in-degree of each node
    for u, adjacent in graph.items():
        for v in adjacent:
            v.in_degree += 1

    return graph


def print_deps(graph):
    for c, deps in graph.items():
        ret = []
        for d in deps:
            ret.append(d.check.name)
        print(c, '->', ret)


def _reduce_deps(graph):
    '''Reduce test case graph to a test-only graph.'''
    ret = {}
    for case, deps in graph.items():
        test_deps = util.OrderedSet(d.check.name for d in deps)
        try:
            ret[case.check.name] |= test_deps
        except KeyError:
            ret[case.check.name] = test_deps

    return ret


def validate_deps(graph):
    '''Validate dependency graph.'''

    # Reduce test case graph to a test name only graph; this disallows
    # pseudo-dependencies as follows:
    #
    # (t0, e1) -> (t1, e1)
    # (t1, e0) -> (t0, e0)
    #
    test_graph = _reduce_deps(graph)

    # Check for cyclic dependencies in the test name graph
    visited = set()
    sources = set(test_graph.keys())
    path = []

    # Since graph may comprise multiple not connected subgraphs, we search for
    # cycles starting from all possible sources
    while sources:
        unvisited = [(sources.pop(), None)]
        while unvisited:
            node, parent = unvisited.pop()
            while path and path[-1] != parent:
                path.pop()

            adjacent = test_graph[node]
            path.append(node)
            for n in adjacent:
                if n in path:
                    cycle_str = '->'.join(path + [n])
                    raise DependencyError(
                        'found cyclic dependency between tests: ' + cycle_str)

                if n not in visited:
                    unvisited.append((n, node))

            visited.add(node)

        sources -= visited


# FIXME This is a replica of the toposort method.
# There must be a cleaverer way of getting this info
def toposortdepth(graph, is_subgraph=False):
    '''Return a list of the graph nodes topologically sorted, with info about
    the node's depht and dependencies.

    If ``is_subgraph`` is ``True``, graph will by treated a subgraph, meaning
    that any dangling edges will be ignored.
    '''
    visited = collections.OrderedDict()
    max_depth = 0

    def retrieve(d, key, default):
        try:
            return d[key]
        except KeyError:
            if is_subgraph:
                return default
            else:
                raise

    def visit(node, path, level):
        # We assume an acyclic graph
        assert node not in path

        path.add(node)
        nonlocal max_depth

        # Do a DFS visit of all the adjacent nodes
        depth = 0
        for adj in graph[node]:
            if adj not in visited:
                visit(adj, path, level+1)
            else:
                depth = max(depth, visited[adj]['depth'] + 1)

        path.pop()
        visited[node] = {
            'depth': max(level, depth),
            'depends' : [adj for adj in graph[node]]
        }
        max_depth = max(level, depth, max_depth)

    for node in graph.keys():
        if node not in visited:
            visit(node, util.OrderedSet(), 0)

    return visited, max_depth


def toposort(graph, is_subgraph=False):
    '''Return a list of the graph nodes topologically sorted.

    If ``is_subgraph`` is ``True``, graph will by treated a subgraph, meaning
    that any dangling edges will be ignored.
    '''
    test_deps = _reduce_deps(graph)
    visited = util.OrderedSet()

    def retrieve(d, key, default):
        try:
            return d[key]
        except KeyError:
            if is_subgraph:
                return default
            else:
                raise

    def visit(node, path):
        # We assume an acyclic graph
        assert node not in path

        path.add(node)

        # Do a DFS visit of all the adjacent nodes
        for adj in retrieve(test_deps, node, []):
            if adj not in visited:
                visit(adj, path)

        path.pop()
        visited.add(node)

    for r in test_deps.keys():
        if r not in visited:
            visit(r, util.OrderedSet())

    # Index test cases by test name
    cases_by_name = {}
    for c in graph.keys():
        try:
            cases_by_name[c.check.name].append(c)
        except KeyError:
            cases_by_name[c.check.name] = [c]

    return list(itertools.chain(*(retrieve(cases_by_name, n, [])
                                  for n in visited)))
