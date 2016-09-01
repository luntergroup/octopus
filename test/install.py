#!/usr/bin/env python3

import os
from subprocess import call
import argparse

def run_unit_tests(octopus_build_dir, use_verbose_output):
    octopus_test_dir = octopus_build_dir + "/test"
    os.chdir(octopus_test_dir)
    ctest_options = []
    if use_verbose_output:
        ctest_options.append("--verbose")
    call(["ctest"] + ctest_options)

parser = argparse.ArgumentParser()
parser.add_argument('--type', help='C++ compiler path', default="unit")
parser.add_argument('--verbose', help='Output verbose test information', action='store_true')
parser.add_argument('--compiler', help='C++ compiler path')
args = vars(parser.parse_args())

if args["type"] not in ["unit", "valgrind", "regression"]:
    print("Unknown test type " + type)
    exit()

# This file is in octopus-dir/test
octopus_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

root_cmake = octopus_dir + "/CMakeLists.txt"

if not os.path.exists(root_cmake):
    print("octopus source directory corrupted: root CMakeLists.txt is missing. Please re-download source code.")
    exit()

octopus_build_dir = octopus_dir + "/build"

if not os.path.exists(octopus_build_dir):
    print("octopus source directory corrupted: build directory is missing. Please re-download source code.")
    exit()

os.chdir(octopus_build_dir) # so cmake doesn't pollute root directory

cmake_options = []

if args["type"] == "unit":
    cmake_options.extend(["-DBUILD_TESTING=ON", octopus_dir])
elif args["type"] == "valgrind":
    cmake_options.append("-DCMAKE_BUILD_TYPE=Debug")

if args["compiler"]:
    cmake_options.append("-DCMAKE_CXX_COMPILER=" + args["compiler"])

ret = call(["cmake"] + cmake_options + [".."])

if ret == 0:
    ret = call(["make"])
    if ret == 0:
        if args["type"] == "unit":
            run_unit_tests(octopus_build_dir, args["verbose"])
        elif args["type"] == "valgrind":
            call(["make", "install"])
