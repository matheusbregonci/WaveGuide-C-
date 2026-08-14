#!/usr/bin/env bash
# Build the Waveguide/Cavity physics as a WebAssembly module.
#
#   source /path/to/emsdk/emsdk_env.sh
#   ./web/build.sh
#
# Produces web/dist/waveguide.{js,wasm}. Only the analytic guide is compiled:
# no Eigen, no OpenMP, no FDTD, no microstrip, no OpenGL. That is what keeps
# this a handful of source files with nothing but the C++ standard library
# behind them -- and why the module needs no SharedArrayBuffer and therefore no
# COOP/COEP headers from whatever serves it.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
out="$here/dist"
mkdir -p "$out"

if ! command -v em++ >/dev/null 2>&1; then
    echo "em++ not on PATH. Run 'source <emsdk>/emsdk_env.sh' first." >&2
    exit 1
fi

em++ \
    -std=c++17 -O3 \
    -I"$root/include" \
    "$here/bindings.cpp" \
    "$root/src/TEmnModel.cpp" \
    "$root/src/CylindricalModel.cpp" \
    "$root/src/FieldViz.cpp" \
    --bind \
    -s MODULARIZE=1 \
    -s EXPORT_NAME=createWaveguide \
    -s ENVIRONMENT=web \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s INITIAL_MEMORY=32MB \
    -s DISABLE_EXCEPTION_CATCHING=1 \
    -o "$out/waveguide.js"

echo
echo "built:"
ls -lh "$out"/waveguide.js "$out"/waveguide.wasm
