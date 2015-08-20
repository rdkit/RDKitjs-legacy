#!/bin/bash
cmake .. \
-DCMAKE_TOOLCHAIN_FILE=/users/mbp/github/emk1.32/emscripten-1.32.0/cmake/Modules/Platform/Emscripten.cmake \
-DRDKIT_INCLUDE_DIR=/users/mbp/github/emk1.32/rdkit-Release_2015_03_1/Code \
-DBoost_INCLUDE_DIR=/usr/local/Cellar/boost155/1.55.0_1/include/  \
-DRDKIT_LIB_DIR=/users/mbp/github/emk1.32/rdkit-Release_2015_03_1/lib \
-DEMSCRIPTEN_BIN=/users/mbp/github/emk1.32/emscripten-1.32.0/ 

rm src/rdkit.js

make

cp src/rdkit.js ../build/src/rdkit.js

cat ../javascript/pre.js ../build/src/rdkit.js ../javascript/post.js > ../dist/rdkit.js

cp ../dist/rdkit.js ../../node_modules/rdkit/dist/rdkit.js

cp ../dist/rdkit.js ../../../Sites/rdkitjs/rdkit.js
