set(RDKITJS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/rdmol.cpp) 
set(RDKITJS_H ${CMAKE_CURRENT_SOURCE_DIR}/rdmol.h)


if(Boost_FOUND)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/rdkit.js
        COMMAND ${EMSCRIPTEN_BIN}/em++ --bind -o rdkit.js -s MODULARIZE=1  -s EXPORT_NAME='"'rdk'"'  ${RDKITJS_SRC}  ${CMAKE_CXX_FLAGS}  -I${RDKIT_INCLUDE_DIR} -I${Boost_INCLUDE_DIR} ${Boost_LIBRARIES} ${RDKIT_LIB_DIR}/libGraphMol.so ${RDKIT_LIB_DIR}/libDescriptors.so ${RDKIT_LIB_DIR}/libRDGeneral.so ${RDKIT_LIB_DIR}/libRDGeometryLib.so ${RDKIT_LIB_DIR}/libSmilesParse.so ${RDKIT_LIB_DIR}/libDataStructs.so ${RDKIT_LIB_DIR}/libFingerprints.so ${RDKIT_LIB_DIR}/libSubgraphs.so ${RDKIT_LIB_DIR}/libDistGeomHelpers.so ${RDKIT_LIB_DIR}/libForceField.so ${RDKIT_LIB_DIR}/libDepictor.so ${RDKIT_LIB_DIR}/libDistGeometry.so ${RDKIT_LIB_DIR}/libEigenSolvers.so ${RDKIT_LIB_DIR}/libAlignment.so ${RDKIT_LIB_DIR}/libForceFieldHelpers.so ${RDKIT_LIB_DIR}/libFileParsers.so ${RDKIT_LIB_DIR}/libSubstructMatch.so ${RDKIT_LIB_DIR}/libPartialCharges.so  ${RDKIT_LIB_DIR}/libMolDraw2D.so ${RDKIT_LIB_DIR}/libChemTransforms.so  ${RDKIT_LIB_DIR}/libAlignment.so ${RDKIT_LIB_DIR}/libMolTransforms.so  ${RDKIT_LIB_DIR}/libMolAlign.so --memory-init-file 0 -s DISABLE_EXCEPTION_CATCHING=0 -s TOTAL_MEMORY=16777216 -s ALLOW_MEMORY_GROWTH=1  -s ASSERTIONS=1 
        )
    add_custom_target(rdkitjs ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/rdkit.js)
endif()

