#!/usr/bin/env python3

import os
import sys
from pathlib import Path
from subprocess import call, check_output
import platform
import distro
import argparse
import shutil
import multiprocessing
import urllib.request
import gzip

google_cloud_octopus_base = "https://storage.googleapis.com/luntergroup/octopus"
forest_url_base = google_cloud_octopus_base + "/forests"
forests = ['germline', 'somatic']

latest_llvm = 'llvm'
latest_gcc = 'gcc'

required_curl_version = 7, 41, 0
required_git_version = 2, 7, 0

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
    cmake_generated_dir = octopus_build_dir / 'generated'
    cmake_version_header = cmake_generated_dir / "version.hpp"
    header = cmake_version_header.open().read().splitlines()
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
        dist, dist_version, _ = distro.linux_distribution()
        if "CentOS" not in dist:
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
    urllib.request.urlretrieve(url, str(file_name))

def git_clone(project):
    call(['git', 'clone', project])

def get_homebrew_name():
    return 'brew'

def download_homebrew():
    git_clone('https://github.com/Homebrew/brew')

def is_old_brew_config_curl(brew_bin):
    brew_config = check_output([str(brew_bin), 'config']).decode("utf-8").split()
    if 'Curl:' in brew_config:
        brew_curl_version = tuple(int(v) for v in brew_config[brew_config.index('Curl:') + 1].split('.'))
        return brew_curl_version < required_curl_version
    else:
        raise Exception("No Curl in brew config")

def curl_version(curl_bin):
    return tuple([int(v) for v in check_output([str(curl_bin), '--version']).decode("utf-8").split('\n')[0].strip().split()[1].split('.')])

def is_old_brew_config_git(brew_bin):
    brew_config = check_output([str(brew_bin), 'config']).decode("utf-8").split()
    if 'Git:' in brew_config:
        brew_git_version = tuple(int(v) for v in brew_config[brew_config.index('Git:') + 1].split('.'))
        return brew_git_version < required_git_version
    else:
        return False

def which(program):
    return check_output(['which', program]).decode("utf-8").strip()

def git_version(git_bin):
    return tuple([int(v) for v in check_output([str(git_bin), '--version']).decode("utf-8").strip().split()[-1].split('.')])

def hack_old_brewed_curl():
    env_curl = which('curl')
    if curl_version(env_curl) < required_curl_version:
        required_curl_version_str = '.'.join([str(v) for v in required_curl_version])
        raise Exception('Could not find cURL version >= ' + required_curl_version_str)
    else:
        os.environ["HOMEBREW_CURL_PATH"] = env_curl
        os.environ["HOMEBREW_NO_ENV_FILTERING"] = "1"

def hack_old_brewed_git():
    env_git = which('git')
    if git_version(env_git) < required_git_version:
        required_git_version_str = '.'.join([str(v) for v in required_git_version])
        raise Exception('Could not find git version >= ' + required_git_version_str)
    else:
        os.environ["HOMEBREW_GIT_PATH"] = env_git
        os.environ["HOMEBREW_NO_ENV_FILTERING"] = "1"

def install_curl_from_source(install_dir):
    curl_dir = install_dir / "curl"
    if not curl_dir.exists():
        call(["wget", "https://github.com/curl/curl/releases/download/curl-7_74_0/curl-7.74.0.tar.gz"])
        call("tar -xzf curl-7.74.0.tar.gz".split())
        Path("curl-7.74.0.tar.gz").unlink()
        shutil.move("curl-7.74.0", curl_dir)
        cwd = Path.cwd()
        os.chdir(curl_dir)
        call(["./configure"])
        call(["make"])
        os.chdir(cwd)
    os.environ['PATH'] = os.pathsep.join([str(curl_dir / "src"), os.environ['PATH']])
    os.environ["HOMEBREW_CURL_PATH"] = str(curl_dir / "src" / "curl")
    os.environ["HOMEBREW_NO_ENV_FILTERING"] = "1"

def install_git_from_source(install_dir):
    git_dir = install_dir / "git"
    if not git_dir.exists():
        call(["git", "clone", "https://github.com/git/git.git"])
        shutil.move("git", install_dir)
        cwd = Path.cwd()
        os.chdir(git_dir)
        call(["make", "configure"])
        call(["./configure"])
        call(["make"])
        os.chdir(cwd)
    os.environ['PATH'] = os.pathsep.join([str(git_dir), os.environ['PATH']])
    os.environ["HOMEBREW_GIT_PATH"] = str(git_dir / "git")
    os.environ["HOMEBREW_NO_ENV_FILTERING"] = "1"
    
def gcc_version(gcc_bin):
    return tuple(int(v) for v in check_output([str(gcc_bin), '-dumpversion']).decode("utf-8").strip().split('.'))

def hack_centos6_brewed_gcc(brew_bin_dir):
    # See https://github.com/Homebrew/linuxbrew-core/issues/4803
    # and https://github.com/Homebrew/linuxbrew-core/issues/4077
    # and https://github.com/Linuxbrew/brew/wiki/Symlink-GCC
    try:
        os.symlink(which('gcc'), brew_bin_dir / ('gcc-' + '.'.join(str(v) for v in gcc_version(which('gcc'))[:2])))
        os.symlink(which('g++'), brew_bin_dir / ('g++-' + '.'.join(str(v) for v in gcc_version(which('g++'))[:2])))
        os.symlink(which('gfortran'), brew_bin_dir / ('gfortran-' + '.'.join(str(v) for v in gcc_version(which('gfortran'))[:2])))
    except FileExistsError:
        return

def init_homebrew(brew_bin_dir):
    brew_bin = brew_bin_dir / 'brew'
    brew_dependency_dir = brew_bin_dir.parent.parent / "brew_dependencies"
    if is_old_brew_config_curl(brew_bin):
        try:
            hack_old_brewed_curl()
        except:
            print("Could not find viable cURL, attempting to install from source...")
            install_curl_from_source(brew_dependency_dir)
    if is_old_brew_config_git(brew_bin):
        try:
            hack_old_brewed_git()
        except:
            print("Could not find viable git, attempting to install from source...")
            install_git_from_source(brew_dependency_dir)
    if is_centos((6,)):
        hack_centos6_brewed_gcc(brew_bin_dir)
    call([str(brew_bin), 'update'])

def install_homebrew(build_dir):
    brew_dir = build_dir / get_homebrew_name()
    if not brew_dir.exists():
        download_homebrew()
    brew_bin_dir = brew_dir / 'bin'
    os.environ['PATH'] = os.pathsep.join([str(brew_bin_dir), str(brew_dir / 'sbin'), os.environ['PATH']])
    init_homebrew(brew_bin_dir)
    return brew_bin_dir / 'brew' # brew binary

def get_glibc_version():
    try:
        return tuple(int(v) for v in check_output(['ldd', '--version']).decode("utf-8").strip().split('\n')[0].split()[-1].split('.'))
    except:
        return True

def is_old_glibc():
    return get_glibc_version() < (2, 18)

def get_required_dependencies():
    compiled, recompile = ['cmake'], []
    if is_osx():
        compiled.append(latest_llvm)
    else:
        compiled.append(latest_gcc)
        compiled += ['binutils']
        if is_old_glibc():
            compiled.append('glibc')
    compiled += ['boost', 'gmp']
    if is_osx():
        compiled += ['htslib']
    else:
        recompile += ['htslib']
    return compiled, recompile

def get_brewed_compiler_binaries(homebrew_dir):
    cellar_dir = homebrew_dir / 'Cellar'
    if is_osx():
        llvm_dir = cellar_dir / latest_llvm
        llvm_version = os.listdir(str(llvm_dir))[0]
        llvm_bin_dir = llvm_dir / llvm_version / 'bin'
        return llvm_bin_dir / 'clang', llvm_bin_dir / 'clang++'
    else:
        gcc_dir = cellar_dir / latest_gcc
        gcc_version = os.listdir(str(gcc_dir))[0]
        gcc_bin_dir = gcc_dir / gcc_version / 'bin'
        return gcc_bin_dir / 'gcc', gcc_bin_dir / 'g++'

def patch_homebrew_centos_gcc9(homebrew_dir):
    homebrew_bin_dir = homebrew_dir / 'bin'
    patchelf_bin = homebrew_bin_dir / 'patchelf'
    homebrew_ld = homebrew_dir / 'lib' / 'ld.so'
    gcc_dir = homebrew_dir / 'Cellar' / latest_gcc
    gcc_version = os.listdir(str(gcc_dir))[0]
    gcc_dir = gcc_dir / gcc_version
    gcc_bin_dir = gcc_dir / 'bin'
    bins = [str(gcc_bin_dir / ex) for ex in os.listdir(str(gcc_bin_dir))]
    for b in bins: os.chmod(b, 0o755)
    call([str(patchelf_bin), "--set-interpreter", str(homebrew_ld)] + bins)
    for b in bins: os.chmod(b, 0o555)
    libexec_dir = gcc_dir / 'libexec' / 'gcc' / 'x86_64-pc-linux-gnu' / gcc_version
    libexec_libs = ['cc1', 'cc1obj', 'cc1objplus', 'cc1plus', 'collect2', 'f951', 'lto1', 'lto-wrapper', 'plugin/gengtype', 'install-tools/fixincl']
    libexec_libs = [str(libexec_dir / f) for f in libexec_libs]
    call([str(patchelf_bin), "--set-interpreter", str(homebrew_ld)] + libexec_libs)

def install_dependencies(build_dir):
    brew_bin = install_homebrew(build_dir)
    compiled_dependencies, recompile_dependencies = get_required_dependencies()
    if call([str(brew_bin), 'install'] + compiled_dependencies) == 0 and \
            call([str(brew_bin), 'install', '--build-from-source'] + recompile_dependencies) == 0:
        dependencies_dir = build_dir / get_homebrew_name()
        if is_centos(): patch_homebrew_centos_gcc9(dependencies_dir)
        cc, cxx = get_brewed_compiler_binaries(dependencies_dir)
        binaries = {'cmake': dependencies_dir / 'bin' / 'cmake',
                    'c_compiler': cc, 'cxx_compiler': cxx}
        return dependencies_dir, binaries
    else:
        return None, None

def unzip_gz(in_filename, out_filename=None):
    with gzip.open(str(in_filename), 'rb') as gzipped_file:
        if out_filename is None:
            out_filename = in_filename.with_suffix('')
        with out_filename.open(mode='wb') as file:
            shutil.copyfileobj(gzipped_file, file)

def download_forests(forest_dir, version):
    if not forest_dir.exists():
        print("No forest directory found, making one")
        forest_dir.mkdir(parents=True)
    for forest in forests:
        gzipped_forest_name = forest + '.v' + to_short_version_str(version) + '.forest.gz'
        gzipped_forest_url = forest_url_base + '/' + gzipped_forest_name
        gzipped_forest_file = forest_dir / gzipped_forest_name
        try:
            print("Downloading " + gzipped_forest_url + " to " + str(gzipped_forest_file))
            download_file(gzipped_forest_url, gzipped_forest_file)
            unzip_gz(gzipped_forest_file)
            gzipped_forest_file.unlink()
        except:
            print("Failed to download forest " + gzipped_forest_name)

def main(args):
    script_dir = Path(__file__).resolve().parent
    octopus_dir = script_dir.parent
    root_cmake = octopus_dir / "CMakeLists.txt"
    
    if not root_cmake.exists():
        print("octopus source directory corrupted: root CMakeLists.txt is missing. Please re-download source code.")
        exit(1)

    octopus_build_dir = octopus_dir / "build"

    if not octopus_build_dir.exists():
        print("octopus source directory corrupted: build directory is missing. Please re-download source code.")
        exit(1)

    bin_dir = octopus_dir / "bin"

    if not bin_dir.exists():
        print("No bin directory found, making one")
        bin_dir.mkdir(parents=True)

    if args["clean"]:
        print("Cleaning build directory")
        (octopus_build_dir / "cmake").rename(octopus_dir / "cmake")
        shutil.rmtree(str(octopus_build_dir))
        octopus_build_dir.mkdir(parents=True)
        (octopus_dir / "cmake").rename(octopus_build_dir / "cmake")

    cmake_cache_file = Path("CMakeCache.txt")
    os.chdir(str(octopus_build_dir)) # so cmake doesn't pollute root directory

    if not args["keep_cache"] and cmake_cache_file.exists():
        cmake_cache_file.unlink()

    dependencies_dir, dependencies_binaries = None, None
    if args["dependencies"]:
        dependencies_dir, dependencies_binaries = install_dependencies(octopus_build_dir)

    cmake_options = []
    if args["prefix"]:
        cmake_options.append("-DCMAKE_INSTALL_PREFIX=" + str(args["prefix"]))
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
    if args["architecture"]:
            cmake_options.append("-DCOMPILER_ARCHITECTURE=" + args["architecture"])
    if dependencies_dir is not None:
        if args["c_compiler"]:
            cmake_options.append("-DCMAKE_C_COMPILER=" + str(args["c_compiler"]))
        else:
            cmake_options.append("-DCMAKE_C_COMPILER=" + str(dependencies_binaries["c_compiler"]))
        if args["cxx_compiler"]:
            cmake_options.append("-DCMAKE_CXX_COMPILER=" + str(args["cxx_compiler"]))
        else:
            cmake_options.append("-DCMAKE_CXX_COMPILER=" + str(dependencies_binaries["cxx_compiler"]))
        if args["boost"]:
            cmake_options.append("-DBOOST_ROOT=" + str(args["boost"]))
            cmake_options.append("-DBoost_NO_BOOST_CMAKE=TRUE")
            cmake_options.append("-DBoost_NO_SYSTEM_PATHS=TRUE")
        else:
            cmake_options.append("-DBOOST_ROOT=" + str(dependencies_dir))
            cmake_options.append("-DBoost_NO_BOOST_CMAKE=TRUE")
            cmake_options.append("-DBoost_NO_SYSTEM_PATHS=TRUE")
        if args["htslib"]:
            cmake_options.append("-DHTSLIB_ROOT=" + str(args["htslib"]))
            cmake_options.append("-DHTSlib_NO_SYSTEM_PATHS=TRUE")
        else:
            cmake_options.append("-DHTSLIB_ROOT=" + str(dependencies_dir))
            cmake_options.append("-DHTSlib_NO_SYSTEM_PATHS=TRUE")
        if args["gmp"]:
            cmake_options.append("-DGMP_ROOT=" + str(args["gmp"]))
        else:
            cmake_options.append("-DGMP_ROOT=" + str(dependencies_dir))

        ret = call([str(dependencies_binaries['cmake'])] + cmake_options + [".."])
    else:
        if args["c_compiler"]:
            cmake_options.append("-DCMAKE_C_COMPILER=" + str(args["c_compiler"]))
        if args["cxx_compiler"]:
            cmake_options.append("-DCMAKE_CXX_COMPILER=" + str(args["cxx_compiler"]))
        if args["boost"]:
            cmake_options.append("-DBOOST_ROOT=" + str(args["boost"]))
            cmake_options.append("-DBoost_NO_BOOST_CMAKE=TRUE")
            cmake_options.append("-DBoost_NO_SYSTEM_PATHS=TRUE")
        if args["htslib"]:
            cmake_options.append("-DHTSLIB_ROOT=" + str(args["htslib"]))
            cmake_options.append("-DHTSlib_NO_SYSTEM_PATHS=TRUE")
        if args["gmp"]:
            cmake_options.append("-DGMP_ROOT=" + str(args["gmp"]))

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

        if args["forests"]:
            if len(forests) > 0:
                forest_dir = octopus_dir / "resources" / "forests"
                download_forests(forest_dir, octopus_version)

    sys.exit(ret)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix',
                        required=False,
                        type=Path,
                        help='Install into given location')
    parser.add_argument('-D', '--dependencies',
                        default=False,
                        help='Install all dependencies locally into build directory',
                        action='store_true')
    parser.add_argument('--clean',
                        default=False,
                        help='Do a clean install',
                        action='store_true')
    parser.add_argument('-c', '--c_compiler',
                        required=False,
                        type=Path,
                        help='C compiler path to use')
    parser.add_argument('-cxx', '--cxx_compiler',
                        required=False,
                        type=Path,
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
                        type=Path,
                        help='The Boost library root')
    parser.add_argument('--htslib',
                        required=False,
                        type=Path,
                        help='The HTSlib library root')
    parser.add_argument('--gmp',
                        required=False,
                        type=Path,
                        help='The GMP library root')
    parser.add_argument('-F', '--forests',
                        default=False,
                        help='Try to download pre-trained random forests for filtering',
                        action='store_true')
    parser.add_argument('--verbose',
                        default=False,
                        help='Ouput verbose make information',
                        action='store_true')
    parser.add_argument('--architecture',
                        required=False,
                        type=str,
                        help='The architecture to compile for')
    args = vars(parser.parse_args())
    main(args)
