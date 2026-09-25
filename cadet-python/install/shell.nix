{ pkgs ? import <nixpkgs> { config.allowUnfree = true; } }:

pkgs.mkShell {
  nativeBuildInputs = with pkgs; [
    cmake
    pkg-config
    gcc
  ];

  buildInputs = with pkgs; [
    mkl
    sundials
    hdf5
    suitesparse # Provides UMFPACK
    superlu
    eigen
  ];

  shellHook = ''
    export MKLROOT="${pkgs.mkl}"
    export CMAKE_PREFIX_PATH="${pkgs.mkl}:${pkgs.sundials}:${pkgs.hdf5}:${pkgs.suitesparse}:${pkgs.superlu}:${pkgs.eigen}:$CMAKE_PREFIX_PATH"
    export PKG_CONFIG_PATH="${pkgs.hdf5}/lib/pkgconfig:${pkgs.suitesparse}/lib/pkgconfig:$PKG_CONFIG_PATH"
  '';
}
