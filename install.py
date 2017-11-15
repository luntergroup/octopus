#!/usr/bin/env python3

import os
import sys
from subprocess import call
import platform
import argparse
from shutil import move, rmtree
import multiprocessing

def is_unix():
    system = platform.system()
    return system == "Darwin" or system == "Linux"

parser = argparse.ArgumentParser()
parser.add_argument('--clean', help='Do a clean install', action='store_true')
parser.add_argument('--root', help='Install into /usr/local/bin', action='store_true')
parser.add_argument('--compiler', help='C++ compiler path')
parser.add_argument('--keep_cache', help='Do not refresh CMake cache', action='store_true')
parser.add_argument('--debug', help='Builds in debug mode', action='store_true')
parser.add_argument('--sanitize', help='Builds in release mode with sanitize flags', action='store_true')
parser.add_argument('--static', help='Builds using static libraries', action='store_true')
parser.add_argument('--threads', help='The number of threads to use for building', type=int)
parser.add_argument('--boost', help='The Boost library root')
args = vars(parser.parse_args())

octopus_dir = os.path.dirname(os.path.realpath(__file__))
root_cmake = octopus_dir + "/CMakeLists.txt"

if not os.path.exists(root_cmake):
    print("octopus source directory corrupted: root CMakeLists.txt is missing. Please re-download source code.")
    exit(1)

octopus_build_dir = octopus_dir + "/build"

if not os.path.exists(octopus_build_dir):
    print("octopus source directory corrupted: build directory is missing. Please re-download source code.")
    exit(1)

bin_dir = octopus_dir + "/bin"

if not os.path.exists(bin_dir):
    print("No bin directory found, making one")
    os.makedirs(bin_dir)

if args["clean"]:
    print("Cleaning build directory")
    move(octopus_build_dir + "/cmake", octopus_dir + "/cmake")
    rmtree(octopus_build_dir)
    os.makedirs(octopus_build_dir)
    move(octopus_dir + "/cmake", octopus_build_dir + "/cmake")

cmake_cache_file = "CMakeCache.txt"
os.chdir(octopus_build_dir) # so cmake doesn't pollute root directory

if not args["keep_cache"] and os.path.exists(cmake_cache_file):
    os.remove(cmake_cache_file)

cmake_options = []
if args["root"]:
    cmake_options.extend(["-DINSTALL_ROOT=ON", octopus_dir])
if args["compiler"]:
    cmake_options.append("-DCMAKE_CXX_COMPILER=" + args["compiler"])
if args["debug"]:
    cmake_options.append("-DCMAKE_BUILD_TYPE=Debug")
elif args["sanitize"]:
    cmake_options.append("-DCMAKE_BUILD_TYPE=RelWithDebInfo")
else:
    cmake_options.append("-DCMAKE_BUILD_TYPE=Release")
if args["static"]:
    cmake_options.append("-DUSE_STATIC_BOOST=ON")
if args["boost"]:
    cmake_options.append("-DBOOST_ROOT=" + args["boost"])

ret = call(["cmake"] + cmake_options + [".."])

if ret == 0:
    make_options = []
    if args["threads"]:
        if (args["threads"] > 1):
            make_options.append("-j" + str(args["threads"]))
    else:
        make_options.append("-j" + str(multiprocessing.cpu_count()))
    
    if is_unix():
        ret = call(["make"] + make_options)
        if ret == 0 and args["root"]:
            ret = call(["make", "install"])
    else:
        print("Windows make files not supported. Build files have been written to " + octopus_build_dir)

sys.exit(ret)
