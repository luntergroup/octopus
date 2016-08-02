This folder contains all of the tests for Octopus. These can be divded into three categories:

NOTE: Many of the tests use real data. In order to run the tests the files specified in 'test_common.h' must be present in your system.

1. Component unit tests: these tests cover functionality requirments of the major components of Octopus. They are designed to ensure expected functionality, especially at edge cases, and avoid common bugs (e.g. off-by-one errors). Note many of the tests here are run on real data.
2. Benchmarks: these tests contain benchmarks for various key components. Generally these are tests that have directed design decisions (e.g. using virtual methods).
3. Data: these are tests on real data, usually 1000G. They are designed to measure and improve calling performance.
