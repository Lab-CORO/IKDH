#!/bin/bash
set -e

# Source emsdk from common install locations
for SDK in ~/emsdk /opt/emsdk; do
    if [ -f "$SDK/emsdk_env.sh" ]; then
        source "$SDK/emsdk_env.sh" --quiet
        break
    fi
done

if ! command -v emcc &>/dev/null; then
    echo "Error: emcc not found. Install emsdk:"
    echo "  git clone https://github.com/emscripten-core/emsdk.git ~/emsdk"
    echo "  cd ~/emsdk && ./emsdk install latest && ./emsdk activate latest"
    exit 1
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/web"

echo "Building IKDH WASM..."
emcc \
    "$ROOT/src/hupf/Calculate.cpp" \
    "$ROOT/src/hupf/ik.cpp" \
    "$ROOT/src/ikdh.cpp" \
    "$ROOT/src/ikdh_wasm.cpp" \
    -I "$ROOT/include" \
    -I "$ROOT/src" \
    -std=c++17 \
    -O2 \
    --bind \
    -s MODULARIZE=1 \
    -s EXPORT_NAME="createIKDHModule" \
    -s ALLOW_MEMORY_GROWTH=1 \
    -o "$OUT/ikdh.js"

echo "Done: web/ikdh.js  web/ikdh.wasm"
