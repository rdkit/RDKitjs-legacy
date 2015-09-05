#!/bin/bash
cmake .. \
-DCMAKE_TOOLCHAIN_FILE=//Users/mbp/Github/emsdk_portable1/emscripten/tag-1.34.6/cmake/Modules/Platform/Emscripten.cmake \
-DRDKIT_INCLUDE_DIR=/Users/mbp/Github/rdkit/Code \
-DBoost_INCLUDE_DIR=/usr/local/Cellar/boost/1.58.0/include/  \
-DRDKIT_LIB_DIR=/Users/mbp/Github/rdkit/lib \
-DEMSCRIPTEN_BIN=/Users/mbp/Github/emsdk_portable1/emscripten/tag-1.34.6

rm src/rdkit.js

make

cp src/rdkit.js ../build/src/rdkit.js


npm run build-test

#cat ../javascript/pre.js ../build/src/rdkit.js ../javascript/post.js > ../dist/rdkit.js
# already include in the build-test method

cp ../test/rdkit.js ../../node_modules/rdkit/dist/rdkit.js

cp ../test/rdkit.js ../../../Sites/rdkitjs/rdkit.js
