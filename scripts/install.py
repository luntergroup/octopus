#!/usr/bin/env python3

import os
import os.path
import sys
from subprocess import call
import platform
import argparse
from shutil import move, rmtree
import multiprocessing
import urllib.request

google_cloud_octopus_base = "https://storage.googleapis.com/luntergroup/octopus"
forest_url_base = os.path.join(google_cloud_octopus_base, "forests")
forests = ['germline', 'somatic']

def get_octopus_version():
    return "0.5.2-beta"

def is_unix():
    system = platform.system()
    return system == "Darwin" or system == "Linux"

def download_file(url, file_name):
    urllib.request.urlretrieve(url, file_name)

def download_forests(forest_dir, version):
    if not os.path.exists(forest_dir):
        print("No forest directory found, making one")
        os.makedirs(forest_dir)
    for forest in forests:
        forest_name = forest + '.v' + version + '.forest'
        forest_url = os.path.join(forest_url_base, forest_name)
        forest_file = os.path.join(forest_dir, forest_name)
        try:
            print("Downloading " + forest_url + " to " + forest_file)
            download_file(forest_url, forest_file)
        except:
            print("Failed to download forest " + forest_name)

def main(args):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    octopus_dir = os.path.dirname(script_dir)
    root_cmake = os.path.join(octopus_dir, "CMakeLists.txt")
    octopus_version = get_octopus_version()

    print("Installing Octopus " + octopus_version)

    if not os.path.exists(root_cmake):
        print("octopus source directory corrupted: root CMakeLists.txt is missing. Please re-download source code.")
        exit(1)

    octopus_build_dir = os.path.join(octopus_dir, "build")

    if not os.path.exists(octopus_build_dir):
        print("octopus source directory corrupted: build directory is missing. Please re-download source code.")
        exit(1)

    bin_dir = os.path.join(octopus_dir, "bin")

    if not os.path.exists(bin_dir):
        print("No bin directory found, making one")
        os.makedirs(bin_dir)

    if args["clean"]:
        print("Cleaning build directory")
        move(os.path.join(octopus_build_dir, "cmake"), os.path.join(octopus_dir, "cmake"))
        rmtree(octopus_build_dir)
        os.makedirs(octopus_build_dir)
        move(os.path.join(octopus_dir, "cmake"), os.path.join(octopus_build_dir, "cmake"))

    cmake_cache_file = "CMakeCache.txt"
    os.chdir(octopus_build_dir) # so cmake doesn't pollute root directory

    if not args["keep_cache"] and os.path.exists(cmake_cache_file):
        os.remove(cmake_cache_file)

    cmake_options = []
    if args["prefix"]:
        cmake_options.append("-DCMAKE_INSTALL_PREFIX=" + args["prefix"])
    if args["c_compiler"]:
        cmake_options.append("-DCMAKE_C_COMPILER=" + args["c_compiler"])
    if args["cxx_compiler"]:
        cmake_options.append("-DCMAKE_CXX_COMPILER=" + args["cxx_compiler"])
    if args["debug"]:
        cmake_options.append("-DCMAKE_BUILD_TYPE=Debug")
    elif args["sanitize"]:
        cmake_options.append("-DCMAKE_BUILD_TYPE=RelWithDebInfo")
    else:
        cmake_options.append("-DCMAKE_BUILD_TYPE=Release")
    if args["static"]:
        cmake_options.append("-DBUILD_SHARED_LIBS=OFF")
    if args["boost"]:
        cmake_options.append("-DBOOST_ROOT=" + args["boost"])
    if args["htslib"]:
        cmake_options.append("-DHTSLIB_ROOT=" + args["htslib"])
    if args["verbose"]:
        cmake_options.append("CMAKE_VERBOSE_MAKEFILE:BOOL=ON")

    try:
        # CMake version 3 is called cmake3 in CentOS (see https://github.com/luntergroup/octopus/issues/37).
        ret = call(["cmake3"] + cmake_options + [".."])
    except FileNotFoundError:
        ret = call(["cmake"] + cmake_options + [".."])

    if ret == 0:
        make_options = []
        if args["threads"]:
            if (args["threads"] > 1):
                make_options.append("-j" + str(args["threads"]))
        else:
            make_options.append("-j" + str(multiprocessing.cpu_count()))

        if is_unix():
            ret = call(["make", "install"] + make_options)
        else:
            print("Windows make files not supported. Build files have been written to " + octopus_build_dir)

    if args["download"]:
        if len(forests) > 0:
            forest_dir = os.path.join(octopus_dir, "resources/forests")
            download_forests(forest_dir, octopus_version)

    sys.exit(ret)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix',
                        required=False,
                        type=str,
                        help='Install into given location')
    parser.add_argument('--clean',
                        help='Do a clean install',
                        action='store_true')
    parser.add_argument('-c', '--c_compiler',
                        required=False,
                        type=str,
                        help='C compiler path to use')
    parser.add_argument('-cxx', '--cxx_compiler',
                        required=False,
                        type=str,
                        help='C++ compiler path to use')
    parser.add_argument('--keep_cache',
                        help='Do not refresh CMake cache',
                        action='store_true')
    parser.add_argument('--debug',
                        help='Builds in debug mode',
                        action='store_true')
    parser.add_argument('--sanitize',
                        help='Builds in release mode with sanitize flags',
                        action='store_true')
    parser.add_argument('--static',
                        help='Builds using static libraries',
                        action='store_true')
    parser.add_argument('--threads',
                        help='The number of threads to use for building',
                        type=int)
    parser.add_argument('--boost',
                        required=False,
                        type=str,
                        help='The Boost library root')
    parser.add_argument('--htslib',
                        required=False,
                        type=str,
                        help='The HTSlib library root')
    parser.add_argument('--download',
                        required=False,
                        help='Try to download octopus classifiers',
                        action='store_true')
    parser.add_argument('--verbose',
                        help='Ouput verbose make information',
                        action='store_true')
    args = vars(parser.parse_args())
    main(args)
