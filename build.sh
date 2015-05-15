#!/bin/bash

PATH_BOOST_INCLUDE="/home/mzasso/boost/include"
PATH_RDKIT="/home/mzasso/rdkit"
PATH_EMSCRIPTEN="/home/mzasso/emsdk_portable/emscripten/master"

$PATH_EMSCRIPTEN/em++  --bind -o rdkit.js rdmol.cpp \
    -I$PATH_RDKIT/Code -I$PATH_BOOST_INCLUDE \
    $PATH_RDKIT/lib/libGraphMol.so $PATH_RDKIT/lib/libDescriptors.so $PATH_RDKIT/lib/libRDGeneral.so $PATH_RDKIT/lib/libRDGeometryLib.so $PATH_RDKIT/lib/libSmilesParse.so $PATH_RDKIT/lib/libDataStructs.so $PATH_RDKIT/lib/libFingerprints.so $PATH_RDKIT/lib/libSubgraphs.so $PATH_RDKIT/lib/libDistGeomHelpers.so $PATH_RDKIT/lib/libForceField.so $PATH_RDKIT/lib/libDepictor.so $PATH_RDKIT/lib/libDistGeometry.so $PATH_RDKIT/lib/libEigenSolvers.so $PATH_RDKIT/lib/libAlignment.so $PATH_RDKIT/lib/libForceFieldHelpers.so $PATH_RDKIT/lib/libFileParsers.so $PATH_RDKIT/lib/libSubstructMatch.so $PATH_RDKIT/lib/libPartialCharges.so  $PATH_RDKIT/lib/libMolDraw2D.so $PATH_RDKIT/lib/libMolTransforms.so $PATH_RDKIT/lib/libChemTransforms.so \
    -O3 --closure 1 --memory-init-file 0
