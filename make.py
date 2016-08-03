#!/usr/bin/env python3

import os
from subprocess import call
import platform

def is_unix():
    system = platform.system()
    return system == "Darwin" or system == "Linux"

octopus_dir = os.path.dirname(os.path.realpath(__file__))

root_cmake = octopus_dir + "/CMakeLists.txt"

if not os.path.exists(root_cmake):
    print("Octopus source directory corrupted: root CMakeLists.txt is missing. Please redownload source code.")
    exit()

octopus_build_dir = octopus_dir + "/build"

if not os.path.exists(octopus_build_dir):
    print("Octopus source directory corrupted: build directory is missing. Please redownload source code.")
    exit()

bin_dir = octopus_dir + "/bin"

if not os.path.exists(bin_dir):
    print("No bin directory found, making one")
    os.makedirs(bin_dir)

os.chdir(octopus_build_dir) # so cmake doesn't pollute root directory

call(["cmake", octopus_dir])

if is_unix():
    call(["make", "install"])
else:
    print("TODO: make for Windows!")
