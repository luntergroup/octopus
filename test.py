#!/usr/bin/env python3

import os
from subprocess import call

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

ret = call(["cmake", "-DBUILD_TESTING=ON", octopus_dir])
if ret == 0:
    ret = call(["make"])
    if ret == 0:
        octopus_test_dir = octopus_build_dir + "/test"
        os.chdir(octopus_test_dir)
        call("ctest")
