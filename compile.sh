# cmake  .. \
# -DCMAKE_TOOLCHAIN_FILE=~/emsdk-portable/emscripten/1.37.35/cmake/Modules/Platform/Emscripten.cmake \
# -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
# -DRDK_BUILD_CPP_TESTS=OFF \
# -DRDK_BUILD_SLN_SUPPORT=OFF \
# -DBoost_INCLUDE_DIR=/home/mzasso/Downloads/boost_1_66_0 \
# -DTHREADS_PTHREAD_ARG=OFF
#-DEIGEN3_INCLUDE_DIR=/path/to/eigen/3.3.4/include/eigen3/ \

alias emcc=/Users/GVALMTGG/Github/emsdk/emscripten/1.37.33/emcc

mkdir -p out
emcc --bind -Os -s WASM=1 -s NODEJS_CATCH_EXIT=0 -s MODULARIZE=1 -s EXPORT_NAME='"'rdk'"' -o out/rdkit.js -Iincludes -Iincludes/rdkit includes/rdkit-lib/libRDKitSmilesParse.so includes/rdkit-lib/libRDKitGraphMol.so \
           includes/rdkit-lib/libRDKitRDGeneral.so includes/rdkit-lib/libRDKitFileParsers.so includes/rdkit-lib/libRDKitDistGeometry.so includes/rdkit-lib/libRDKitForceFieldHelpers.so \
           includes/rdkit-lib/libRDKitRDGeometryLib.so includes/rdkit-lib/libRDKitDistGeomHelpers.so includes/rdkit-lib/libRDKitForceField.so \
           includes/rdkit-lib/libRDKitEigenSolvers.so includes/rdkit-lib/libRDKitAlignment.so includes/rdkit-lib/libRDKitSubstructMatch.so \
           rdkit.cc
