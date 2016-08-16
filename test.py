#!/usr/bin/env python3

import os
from subprocess import call
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--verbose', help='Output verbose test information', action='store_true')
parser.add_argument('--compiler', help='C++ compiler path')
args = vars(parser.parse_args())

octopus_dir = os.path.dirname(os.path.realpath(__file__))

root_cmake = octopus_dir + "/CMakeLists.txt"

if not os.path.exists(root_cmake):
    print("octopus source directory corrupted: root CMakeLists.txt is missing. Please redownload source code.")
    exit()

octopus_build_dir = octopus_dir + "/build"

if not os.path.exists(octopus_build_dir):
    print("octopus source directory corrupted: build directory is missing. Please redownload source code.")
    exit()

os.chdir(octopus_build_dir) # so cmake doesn't pollute root directory

cmake_options = ["-DBUILD_TESTING=ON", octopus_dir]

if args["compiler"]:
    cmake_options.append("-DCMAKE_CXX_COMPILER=" + args["compiler"])

ret = call(["cmake"] + cmake_options + [".."])
if ret == 0:
    ret = call(["make"])
    if ret == 0:
        octopus_test_dir = octopus_build_dir + "/test"
        os.chdir(octopus_test_dir)
        ctest_options = []
        if (args["verbose"]):
            ctest_options.append("--verbose")
        call(["ctest"] + ctest_options)
