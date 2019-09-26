#!/usr/bin/env python3

import os
import os.path
import sys
from subprocess import call, check_output
import platform
import argparse
from shutil import move, rmtree
import multiprocessing
import urllib.request

google_cloud_octopus_base = "https://storage.googleapis.com/luntergroup/octopus"
forest_url_base = os.path.join(google_cloud_octopus_base, "forests")
forests = ['germline', 'somatic']

latest_llvm = 'llvm'
latest_gcc = 'gcc@9'

required_git_version = 1, 8, 0

class Version(object):
    def __init__(self):
        self.major = None
        self.minor = None
        self.patch = None
        self.release = None
        self.branch = None
        self.commit = None

    def __str__(self):
        res = '.'.join(str(x) for x in [self.major, self.minor, self.patch])
        if self.release:
            res += '-' + self.release
        if self.branch or self.commit:
            res += ' ('
            pad = ''
            if self.branch:
                res += self.branch
                pad = ' '
            if self.commit:
                res += pad + self.commit
            res += ')'
        return res

def to_short_version_str(version):
    res = '.'.join(str(x) for x in [version.major, version.minor, version.patch])
    if version.release:
        res += '-' + version.release
    return res

def get_octopus_version(octopus_build_dir):
    cmake_generated_dir = os.path.join(octopus_build_dir, 'generated')
    cmake_version_header = os.path.join(cmake_generated_dir, "version.hpp")
    header = open(cmake_version_header).read().splitlines()
    result = Version()
    for line in header:
        tokens = line.split()
        if len(tokens) == 3 and tokens[0] == "#define":
            field, value = tokens[1], tokens[2].replace('"', '')
            if field == "VERSION_MAJOR":
                result.major = int(value)
            elif field == "VERSION_MINOR":
                result.minor = int(value)
            elif field == "VERSION_PATCH":
                result.patch = int(value)
            elif field == "VERSION_RELEASE":
                result.release = value
            elif field == "GIT_BRANCH":
                result.branch = value
            elif field == "GIT_COMMIT_HASH":
                result.commit = value
    return result

def is_unix():
    system = platform.system()
    return system == "Darwin" or system == "Linux"

def is_osx():
    return platform.system() == "Darwin"

def is_centos(version=None):
    if platform.system() == "Linux":
        dist, dist_version, _ = platform.linux_distribution()
        if dist != "CentOS":
            return False
        if version is None:
            return True
        else:
            dist_version = tuple(int(v) for v in dist_version.split('.'))
            if len(version) < len(dist_version):
                dist_version = dist_version[:len(version)]
            if len(dist_version) < len(version):
                dist_version = tuple(list(dist_version) + (len(version) - len(dist_version)) * [0])
            return dist_version == version
    else:
        return False

def download_file(url, file_name):
    urllib.request.urlretrieve(url, file_name)

def git_clone(project):
    call(['git', 'clone', project])

def get_homebrew_name():
    return 'brew'

def download_homebrew():
    git_clone('https://github.com/Homebrew/brew')

def is_old_brew_config_git(brew_bin):
    brew_config = check_output([brew_bin, 'config']).decode("utf-8").split()
    if 'Git:' in brew_config:
        brew_git_version = tuple(int(v) for v in brew_config[brew_config.index('Git:') + 1].split('.'))
        return brew_git_version < required_git_version
    else:
        return False

def which(program):
    return check_output(['which', program]).decode("utf-8").strip()

def git_version(git_bin):
    return tuple([int(v) for v in check_output([git_bin, '--version']).decode("utf-8").strip().split()[-1].split('.')])

def hack_old_brewed_git():
    env_git = which('git')
    if git_version(env_git) < required_git_version:
        required_git_version_str = '.'.join([str(v) for v in required_git_version])
        raise Exception('Could not find git version >= ' + required_git_version_str)
    else:
        os.environ["HOMEBREW_GIT_PATH"] = env_git
        os.environ["HOMEBREW_NO_ENV_FILTERING"] = "1"

def gcc_version(gcc_bin):
    return tuple(int(v) for v in check_output([gcc_bin, '-dumpversion']).decode("utf-8").strip().split('.'))

def hack_centos6_brewed_gcc(brew_bin_dir):
    # See https://github.com/Homebrew/linuxbrew-core/issues/4803
    # and https://github.com/Homebrew/linuxbrew-core/issues/4077
    # and https://github.com/Linuxbrew/brew/wiki/Symlink-GCC
    try:
        os.symlink(which('gcc'), os.path.join(brew_bin_dir, 'gcc-' + '.'.join(str(v) for v in gcc_version(which('gcc'))[:2])))
        os.symlink(which('g++'), os.path.join(brew_bin_dir, 'g++-' + '.'.join(str(v) for v in gcc_version(which('g++'))[:2])))
        os.symlink(which('gfortran'), os.path.join(brew_bin_dir, 'gfortran-' + '.'.join(str(v) for v in gcc_version(which('gfortran'))[:2])))
    except FileExistsError:
        return

def init_homebrew(brew_bin_dir):
    brew_bin = os.path.join(brew_bin_dir, 'brew')
    if is_old_brew_config_git(brew_bin):
        hack_old_brewed_git()
    if is_centos((6,)):
        hack_centos6_brewed_gcc(brew_bin_dir)
    call([brew_bin, 'update'])

def install_homebrew(build_dir):
    brew_dir = os.path.join(build_dir, get_homebrew_name())
    if not os.path.exists(brew_dir):
        download_homebrew()
    brew_bin_dir = os.path.join(brew_dir, 'bin')
    os.environ['PATH']= brew_bin_dir + os.pathsep + os.path.join(brew_dir, 'sbin') + os.pathsep + os.environ['PATH']
    init_homebrew(brew_bin_dir)
    return os.path.join(brew_bin_dir, 'brew') # brew binary

def get_required_dependencies():
    result = ['cmake']
    if is_osx():
        result.append(latest_llvm)
    else:
        result.append(latest_gcc)
        result.append('binutils')
    result += ['boost', 'htslib', 'gmp']
    return result

def get_brewed_compiler_binaries(homebrew_dir):
    cellar_dir = os.path.join(homebrew_dir, 'Cellar')
    if is_osx():
        llvm_dir = os.path.join(cellar_dir, latest_llvm)
        llvm_version = os.listdir(llvm_dir)[0]
        llvm_bin_dir = os.path.join(llvm_dir, os.path.join(llvm_version, 'bin'))
        return os.path.join(llvm_bin_dir, 'clang'), os.path.join(llvm_bin_dir, 'clang++')
    else:
        gcc_dir = os.path.join(cellar_dir, latest_gcc)
        gcc_version = os.listdir(gcc_dir)[0]
        gcc_bin_dir = os.path.join(gcc_dir, os.path.join(gcc_version, 'bin'))
        gcc_bin_name = latest_gcc.replace('@', '-')
        gxx_bin_name =  gcc_bin_name.replace('cc', '++')
        return os.path.join(gcc_bin_dir, gcc_bin_name), os.path.join(gcc_bin_dir, gxx_bin_name)

def install_dependencies(build_dir):
    brew_bin = install_homebrew(build_dir)
    dependencies = get_required_dependencies()
    brew_command = [brew_bin, 'install'] + dependencies
    if call(brew_command) == 0:
        dependencies_dir = os.path.join(build_dir, get_homebrew_name())
        cc, cxx = get_brewed_compiler_binaries(dependencies_dir)
        binaries = {'cmake': os.path.join(dependencies_dir, 'bin/cmake'),
                    'c_compiler': cc, 'cxx_compiler': cxx}
        return dependencies_dir, binaries
    else:
        return None, None

def download_forests(forest_dir, version):
    if not os.path.exists(forest_dir):
        print("No forest directory found, making one")
        os.makedirs(forest_dir)
    for forest in forests:
        forest_name = forest + '.v' + to_short_version_str(version) + '.forest'
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

    dependencies_dir, dependencies_binaries = None, None
    if args["install_dependencies"]:
        dependencies_dir, dependencies_binaries = install_dependencies(octopus_build_dir)

    cmake_options = []
    if args["prefix"]:
        cmake_options.append("-DCMAKE_INSTALL_PREFIX=" + args["prefix"])
    if args["debug"]:
        cmake_options.append("-DCMAKE_BUILD_TYPE=Debug")
    elif args["sanitize"]:
        cmake_options.append("-DCMAKE_BUILD_TYPE=RelWithDebInfo")
    else:
        cmake_options.append("-DCMAKE_BUILD_TYPE=Release")
    if args["static"]:
        cmake_options.append("-DBUILD_SHARED_LIBS=OFF")
    if args["verbose"]:
        cmake_options.append("-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON")
    if dependencies_dir is not None:
        if args["c_compiler"]:
            cmake_options.append("-DCMAKE_C_COMPILER=" + args["c_compiler"])
        else:
            cmake_options.append("-DCMAKE_C_COMPILER=" + dependencies_binaries["c_compiler"])
        if args["cxx_compiler"]:
            cmake_options.append("-DCMAKE_CXX_COMPILER=" + args["cxx_compiler"])
        else:
            cmake_options.append("-DCMAKE_CXX_COMPILER=" + dependencies_binaries["cxx_compiler"])
        if args["boost"]:
            cmake_options.append("-DBOOST_ROOT=" + args["boost"])
            cmake_options.append("-DBoost_NO_BOOST_CMAKE=TRUE")
            cmake_options.append("-DBoost_NO_SYSTEM_PATHS=TRUE")
        else:
            cmake_options.append("-DBOOST_ROOT=" + dependencies_dir)
            cmake_options.append("-DBoost_NO_BOOST_CMAKE=TRUE")
            cmake_options.append("-DBoost_NO_SYSTEM_PATHS=TRUE")
        if args["htslib"]:
            cmake_options.append("-DHTSLIB_ROOT=" + args["htslib"])
            cmake_options.append("-DHTSlib_NO_SYSTEM_PATHS=TRUE")
        else:
            cmake_options.append("-DHTSLIB_ROOT=" + dependencies_dir)
            cmake_options.append("-DHTSlib_NO_SYSTEM_PATHS=TRUE")
        if args["gmp"]:
            cmake_options.append("-DGMPlib_ROOT=" + args["gmp"])
        else:
            cmake_options.append("-DGMPlib_ROOT=" + dependencies_dir)

        ret = call([dependencies_binaries['cmake']] + cmake_options + [".."])
    else:
        if args["c_compiler"]:
            cmake_options.append("-DCMAKE_C_COMPILER=" + args["c_compiler"])
        if args["cxx_compiler"]:
            cmake_options.append("-DCMAKE_CXX_COMPILER=" + args["cxx_compiler"])
        if args["boost"]:
            cmake_options.append("-DBOOST_ROOT=" + args["boost"])
            cmake_options.append("-DBoost_NO_BOOST_CMAKE=TRUE")
            cmake_options.append("-DBoost_NO_SYSTEM_PATHS=TRUE")
        if args["htslib"]:
            cmake_options.append("-DHTSLIB_ROOT=" + args["htslib"])
            cmake_options.append("-DHTSlib_NO_SYSTEM_PATHS=TRUE")
        if args["gmp"]:
            cmake_options.append("-DGMPlib_ROOT=" + args["gmp"])

        try:
            # CMake version 3 is called cmake3 in CentOS (see https://github.com/luntergroup/octopus/issues/37).
            ret = call(["cmake3"] + cmake_options + [".."])
        except FileNotFoundError:
            ret = call(["cmake"] + cmake_options + [".."])

    if ret == 0:
        octopus_version = get_octopus_version(octopus_build_dir)
        print("Installing Octopus " + str(octopus_version))

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

        if args["download_forests"]:
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
    parser.add_argument('--install-dependencies',
                        default=False,
                        help='Install all dependencies locally into build directory',
                        action='store_true')
    parser.add_argument('--clean',
                        default=False,
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
                        default=False,
                        help='Do not refresh CMake cache',
                        action='store_true')
    parser.add_argument('--debug',
                        default=False,
                        help='Builds in debug mode',
                        action='store_true')
    parser.add_argument('--sanitize',
                        default=False,
                        help='Builds in release mode with sanitize flags',
                        action='store_true')
    parser.add_argument('--static',
                        default=False,
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
    parser.add_argument('--gmp',
                        required=False,
                        type=str,
                        help='The GMP library root')
    parser.add_argument('--download-forests',
                        default=False,
                        help='Try to download pre-trained random forests for filtering',
                        action='store_true')
    parser.add_argument('--verbose',
                        default=False,
                        help='Ouput verbose make information',
                        action='store_true')
    args = vars(parser.parse_args())
    main(args)
