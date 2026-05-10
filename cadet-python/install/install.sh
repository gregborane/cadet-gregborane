#!/bin/bash
set -euo pipefail

FLAG_INFO="${1:---source_build}" # --source_build | --conda | --both

if [[ "$FLAG_INFO" == "--conda" || "$FLAG_INFO" == "--both" ]]; then
    CONDA_BASE=$(conda info --base)
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda env create -f ./environment.yml || true
fi

if [[ "$FLAG_INFO" == "--source_build" || "$FLAG_INFO" == "--both" ]]; then
    # Best use lastest version of cadet-core by manually compiling it from github source code
    # Extracted from the git code-core build-linux.sh

    # refer to your distro for adapted installation

    ## Dependencies
    sudo pacman -S --noconfirm cmake hdf5 superlu superlu_dist eigen3 git lapack blas openmp

    ## Install Cadet-Core
    INSTALL_PATH="$HOME/Work/cadet"

    if [[ -d "$INSTALL_PATH/CADET-Core" ]]; then
        rm -rf "$INSTALL_PATH/CADET-Core"
    fi
    git clone https://github.com/cadet/cadet-core.git "$INSTALL_PATH/CADET-Core"

    cd "$INSTALL_PATH/CADET-Core" || exit 1
    mkdir -p build install

    cd build || exit 1
    cmake \
        -DCMAKE_INSTALL_PREFIX="../install" ..

    make -j$(nproc)
    make install
fi
