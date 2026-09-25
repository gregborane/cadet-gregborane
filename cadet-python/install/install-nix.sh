INSTALL_PATH="$HOME/Work/cadet"
git clone https://github.com/cadet/cadet-core.git "$INSTALL_PATH/CADET-Core"

cp shell.nix "$INSTALL_PATH/CADET-Core"

cd /home/greg/Work/cadet/CADET-Core

# 1. Clean previous CMake cache completely
rm -rf build install

# 2. Enter the isolated Nix build environment
nix-shell --run '
  mkdir build && cd build
  cmake -DCMAKE_INSTALL_PREFIX="../install" \
        -DBLA_VENDOR=Intel10_64lp \
        -DHDF5_PREFER_PARALLEL=OFF \
        ../
  make -j$(nproc)
  make install
'
