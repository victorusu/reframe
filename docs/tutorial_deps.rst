===============================================
Tutorial 2: Using Dependencies in ReFrame Tests
===============================================

.. versionadded:: 3.3


A ReFrame test may define dependencies to other tests.
An example scenario is to test different runtime configurations of a benchmark that you need to compile, or run a scaling analysis of a code.
In such cases, you don't want to download and rebuild your test for each runtime configuration.
You could have a test where you fetches the sources, which all build tests would depend on.
And similarly, all the runtime tests would depend on their corresponding build test.
This is the approach we take with the following example, that fetches, builds and runs several `OSU benchmarks <http://mvapich.cse.ohio-state.edu/benchmarks/>`__.
We first create a basic run-only test, that fetches the benchmarks:

.. literalinclude:: ../tutorials/deps/osu_benchmarks.py
   :lines: 112-123

This test doesn't need any specific programming environment, we pick the `gnu` environment in the `login` partition.
The build tests would then copy the benchmark and build them for the different programming environments:

.. literalinclude:: ../tutorials/deps/osu_benchmarks.py
   :lines: 93-109

The new part comes in with the :class:`OSUBuildTest` test in the following line:

.. literalinclude:: ../tutorials/deps/osu_benchmarks.py
   :lines: 99

Here we tell ReFrame that this test depends on a test named ``OSUDownloadTest``.
This test may or may not be defined in the same test file; all ReFrame needs is the test name.
The :func:`depends_on() <reframe.core.pipeline.RegressionTest.depends_on>` function will create dependencies between the individual test cases of the :class:`OSUBuildTest` and the :class:`OSUDownloadTest`, such that all the instances of :class:`OSUBuildTest` will depend on the outcome of the :class:`OSUDownloadTest` using ``gnu``, but not on the outcome of the other :class:`OSUBuildTest` instances.
This behaviour can be changed, but it is covered in detail in :doc:`dependencies`.
You can create arbitrary test dependency graphs, but they need to be acyclic.
If ReFrame detects cyclic dependencies, it will refuse to execute the set of tests and will issue an error pointing out the cycle.

A ReFrame test with dependencies will execute, i.e., enter its `setup` stage, only after `all` of its dependencies have succeeded.
If any of its dependencies fails, the current test will be marked as failure as well.

The next step for the :class:`OSUBuildTest` is to set its sourcesdir to point to the benchmarks that were fetched by the :class:`OSUDownloadTest`.
This is achieved with the following specially decorated function:

.. literalinclude:: ../tutorials/deps/osu_benchmarks.py
   :lines: 104-109

The :func:`@require_deps <reframe.core.decorators.require_deps>` decorator will bind the arguments passed to the decorated function to the result of the dependency that each argument names.
In this case, it binds the ``OSUDownloadTest`` function argument to the result of a dependency named ``OSUDownloadTest``.
In order for the binding to work correctly the function arguments must be named after the target dependencies.
However, referring to a dependency only by the test's name is not enough, since a test might be associated with multiple programming environments.
For this reason, a dependency argument is actually bound to a function that accepts as argument the name of a target programming environment.
If no arguments are passed to that function, as in this example, the current programming environment is implied, such that ``OSUDownloadTest()`` is equivalent to ``OSUDownloadTest(self.current_environ.name, self.current_partition.name)``.
This call returns the actual test case of the dependency that has been executed.
This allows you to access any attribute from the target test, as we do in this example by accessing the target test's stage directory, which we use to construct the sourcesdir of the test.

For the next test we need to use the OSU benchmark binaries that we just built, so as to run the MPI ping-pong benchmark.
Here is the relevant part:

.. literalinclude:: ../tutorials/deps/osu_benchmarks.py
   :lines: 13-44

First, since we will have multiple similar benchmarks, we move all the common functionality to the :class:`OSUBenchmarkTestBase` base class.
Again nothing new here; we are going to use two nodes for the benchmark and we set :attr:`sourcesdir <reframe.core.pipeline.RegressionTest.sourcesdir>` to ``None``, since none of the benchmark tests will use any additional resources.
Similar to before, we will define the dependencies with the the following line:

.. literalinclude:: ../tutorials/deps/osu_benchmarks.py
   :lines: 32

Here we tell ReFrame that this test depends on a test named ``OSUBuildTest`` ``by_env``.
This means that it will depend only on the instances of this test that have the same environment.
The partition will be different in this case.

The next step for the :class:`OSULatencyTest` is to set its executable to point to the binary produced by the :class:`OSUBuildTest`.
This is achieved with the following specially decorated function:

.. literalinclude:: ../tutorials/deps/osu_benchmarks.py
   :lines: 38-44

This concludes the presentation of the :class:`OSULatencyTest` test. The :class:`OSUBandwidthTest` is completely analogous.

The :class:`OSUAllreduceTest` shown below is similar to the other two, except that it is parameterized.
It is essentially a scalability test that is running the ``osu_allreduce`` executable created by the :class:`OSUBuildTest` for 2, 4, 8 and 16 nodes.

.. literalinclude:: ../tutorials/deps/osu_benchmarks.py
   :lines: 70-90

The full set of OSU example tests is shown below:

.. literalinclude:: ../tutorials/deps/osu_benchmarks.py

Notice that the order in which dependencies are defined in a test file is irrelevant.
In this case, we define :class:`OSUBuildTest` at the end.
ReFrame will make sure to properly sort the tests and execute them.

Here is the output when running the OSU tests with the asynchronous execution policy:

.. code-block:: none

   [ReFrame Setup]
     version:           3.3-dev1 (rev: 734a53df)
     command:           './bin/reframe -c tutorials/deps/osu_benchmarks.py -C tutorials/config/settings.py -r'
     launched by:       user@daint103
     working directory: '/path/to/reframe'
     settings file:     'tutorials/config/settings.py'
     check search path: '/path/to/reframe/tutorials/deps/osu_benchmarks.py'
     stage directory:   '/path/to/reframe/stage'
     output directory:  '/path/to/reframe/output'

   [==========] Running 8 check(s)
   [==========] Started on Tue Nov  3 13:20:28 2020

   [----------] started processing OSUDownloadTest (OSU benchmarks download sources)
   [ RUN      ] OSUDownloadTest on daint:login using gnu
   [----------] finished processing OSUDownloadTest (OSU benchmarks download sources)

   [----------] started processing OSUBuildTest (OSU benchmarks build test)
   [ RUN      ] OSUBuildTest on daint:login using gnu
   [      DEP ] OSUBuildTest on daint:login using gnu
   [ RUN      ] OSUBuildTest on daint:login using intel
   [      DEP ] OSUBuildTest on daint:login using intel
   [ RUN      ] OSUBuildTest on daint:login using pgi
   [      DEP ] OSUBuildTest on daint:login using pgi
   [----------] finished processing OSUBuildTest (OSU benchmarks build test)

   [----------] started processing OSULatencyTest (OSU latency test)
   [ RUN      ] OSULatencyTest on daint:gpu using gnu
   [      DEP ] OSULatencyTest on daint:gpu using gnu
   [ RUN      ] OSULatencyTest on daint:gpu using intel
   [      DEP ] OSULatencyTest on daint:gpu using intel
   [ RUN      ] OSULatencyTest on daint:gpu using pgi
   [      DEP ] OSULatencyTest on daint:gpu using pgi
   [----------] finished processing OSULatencyTest (OSU latency test)

   [----------] started processing OSUBandwidthTest (OSU bandwidth test)
   [ RUN      ] OSUBandwidthTest on daint:gpu using gnu
   [      DEP ] OSUBandwidthTest on daint:gpu using gnu
   [ RUN      ] OSUBandwidthTest on daint:gpu using intel
   [      DEP ] OSUBandwidthTest on daint:gpu using intel
   [ RUN      ] OSUBandwidthTest on daint:gpu using pgi
   [      DEP ] OSUBandwidthTest on daint:gpu using pgi
   [----------] finished processing OSUBandwidthTest (OSU bandwidth test)

   [----------] started processing OSUAllreduceTest_2 (OSU Allreduce test)
   [ RUN      ] OSUAllreduceTest_2 on daint:gpu using gnu
   [      DEP ] OSUAllreduceTest_2 on daint:gpu using gnu
   [ RUN      ] OSUAllreduceTest_2 on daint:gpu using intel
   [      DEP ] OSUAllreduceTest_2 on daint:gpu using intel
   [ RUN      ] OSUAllreduceTest_2 on daint:gpu using pgi
   [      DEP ] OSUAllreduceTest_2 on daint:gpu using pgi
   [----------] finished processing OSUAllreduceTest_2 (OSU Allreduce test)

   [----------] started processing OSUAllreduceTest_4 (OSU Allreduce test)
   [ RUN      ] OSUAllreduceTest_4 on daint:gpu using gnu
   [      DEP ] OSUAllreduceTest_4 on daint:gpu using gnu
   [ RUN      ] OSUAllreduceTest_4 on daint:gpu using intel
   [      DEP ] OSUAllreduceTest_4 on daint:gpu using intel
   [ RUN      ] OSUAllreduceTest_4 on daint:gpu using pgi
   [      DEP ] OSUAllreduceTest_4 on daint:gpu using pgi
   [----------] finished processing OSUAllreduceTest_4 (OSU Allreduce test)

   [----------] started processing OSUAllreduceTest_8 (OSU Allreduce test)
   [ RUN      ] OSUAllreduceTest_8 on daint:gpu using gnu
   [      DEP ] OSUAllreduceTest_8 on daint:gpu using gnu
   [ RUN      ] OSUAllreduceTest_8 on daint:gpu using intel
   [      DEP ] OSUAllreduceTest_8 on daint:gpu using intel
   [ RUN      ] OSUAllreduceTest_8 on daint:gpu using pgi
   [      DEP ] OSUAllreduceTest_8 on daint:gpu using pgi
   [----------] finished processing OSUAllreduceTest_8 (OSU Allreduce test)

   [----------] started processing OSUAllreduceTest_16 (OSU Allreduce test)
   [ RUN      ] OSUAllreduceTest_16 on daint:gpu using gnu
   [      DEP ] OSUAllreduceTest_16 on daint:gpu using gnu
   [ RUN      ] OSUAllreduceTest_16 on daint:gpu using intel
   [      DEP ] OSUAllreduceTest_16 on daint:gpu using intel
   [ RUN      ] OSUAllreduceTest_16 on daint:gpu using pgi
   [      DEP ] OSUAllreduceTest_16 on daint:gpu using pgi
   [----------] finished processing OSUAllreduceTest_16 (OSU Allreduce test)

   [----------] waiting for spawned checks to finish
   [       OK ] ( 1/22) OSUDownloadTest on daint:login using gnu [compile: 0.005s run: 3.373s total: 3.408s]
   [       OK ] ( 2/22) OSUBuildTest on daint:login using gnu [compile: 22.410s run: 0.035s total: 87.728s]
   [       OK ] ( 3/22) OSUBuildTest on daint:login using pgi [compile: 27.725s run: 59.918s total: 87.691s]
   [       OK ] ( 4/22) OSUBuildTest on daint:login using intel [compile: 37.437s run: 32.771s total: 98.034s]
   [       OK ] ( 5/22) OSUAllreduceTest_2 on daint:gpu using pgi [compile: 0.007s run: 139.339s total: 144.402s]
   [       OK ] ( 6/22) OSUAllreduceTest_4 on daint:gpu using pgi [compile: 0.007s run: 140.896s total: 144.395s]
   [       OK ] ( 7/22) OSUAllreduceTest_8 on daint:gpu using pgi [compile: 0.006s run: 142.451s total: 144.372s]
   [       OK ] ( 8/22) OSUAllreduceTest_16 on daint:gpu using pgi [compile: 0.005s run: 144.342s total: 144.368s]
   [       OK ] ( 9/22) OSUAllreduceTest_4 on daint:gpu using gnu [compile: 0.007s run: 140.555s total: 144.833s]
   [       OK ] (10/22) OSUAllreduceTest_16 on daint:gpu using gnu [compile: 0.005s run: 143.642s total: 144.778s]
   [       OK ] (11/22) OSUAllreduceTest_8 on daint:gpu using gnu [compile: 0.007s run: 142.456s total: 145.151s]
   [       OK ] (12/22) OSUAllreduceTest_2 on daint:gpu using gnu [compile: 0.005s run: 139.685s total: 145.510s]
   [       OK ] (13/22) OSUBandwidthTest on daint:gpu using gnu [compile: 0.009s run: 193.440s total: 200.818s]
   [       OK ] (14/22) OSUBandwidthTest on daint:gpu using pgi [compile: 0.006s run: 194.465s total: 201.080s]
   [       OK ] (15/22) OSULatencyTest on daint:gpu using intel [compile: 0.009s run: 278.603s total: 283.389s]
   [       OK ] (16/22) OSUAllreduceTest_4 on daint:gpu using intel [compile: 0.006s run: 281.112s total: 283.365s]
   [       OK ] (17/22) OSULatencyTest on daint:gpu using pgi [compile: 0.006s run: 285.499s total: 293.712s]
   [       OK ] (18/22) OSUAllreduceTest_2 on daint:gpu using intel [compile: 0.006s run: 280.693s total: 283.756s]
   [       OK ] (19/22) OSUAllreduceTest_8 on daint:gpu using intel [compile: 0.006s run: 282.550s total: 283.971s]
   [       OK ] (20/22) OSUAllreduceTest_16 on daint:gpu using intel [compile: 0.005s run: 284.573s total: 284.596s]
   [       OK ] (21/22) OSULatencyTest on daint:gpu using gnu [compile: 0.006s run: 286.186s total: 295.202s]
   [       OK ] (22/22) OSUBandwidthTest on daint:gpu using intel [compile: 0.005s run: 340.005s total: 343.927s]
   [----------] all spawned checks have finished

   [  PASSED  ] Ran 22 test case(s) from 8 check(s) (0 failure(s))
   [==========] Finished on Tue Nov  3 13:27:54 2020
   Log file(s) saved in: '/tmp/rfm-n4lrqiqf.log'

Before starting running the tests, ReFrame topologically sorts them based on their dependencies and schedules them for running using the selected execution policy.
With the serial execution policy, ReFrame simply executes the tests to completion as they "arrive", since the tests are already topologically sorted.
In the asynchronous execution policy, tests are spawned and not waited for.
If a test's dependencies have not yet completed, it will not start its execution and a ``DEP`` message will be printed to denote this.

Finally, ReFrame's runtime takes care of properly cleaning up the resources of the tests respecting dependencies.
Normally when an individual test finishes successfully, its stage directory is cleaned up.
However, if other tests are depending on this one, this would be catastrophic, since most probably the dependent tests would need the outcome of this test.
ReFrame fixes that by not cleaning up the stage directory of a test until all its dependent tests have finished successfully.
