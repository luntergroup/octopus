#!/bin/bash

set -e
set -x

if [[ "$(uname -s)" == 'Darwin' ]]; then
    brew update || brew update
    brew unlink cmake
    brew upgrade cmake
    brew install gcc || brew link --overwrite gcc
fi
