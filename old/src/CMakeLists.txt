set(RDKITJS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/rdmol.cpp)
set(RDKITJS_H ${CMAKE_CURRENT_SOURCE_DIR}/rdmol.h)


# boost found

if(Boost_FOUND)
	# RDkit Original lib format
	if(NOT RDKitNEWLIB)
	    add_custom_command(
	        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/rdkit.js
	        COMMAND ${EMSCRIPTEN_BIN}/em++ --bind -o rdkit.js -s MODULARIZE=1  -s EXPORT_NAME='"'rdk'"'  ${RDKITJS_SRC}  ${CMAKE_CXX_FLAGS}  -I${RDKIT_INCLUDE_DIR} -I${Boost_INCLUDE_DIR} ${Boost_LIBRARIES} ${RDKIT_LIB_DIR}/libGraphMol.so ${RDKIT_LIB_DIR}/libDescriptors.so ${RDKIT_LIB_DIR}/libRDGeneral.so ${RDKIT_LIB_DIR}/libRDGeometryLib.so ${RDKIT_LIB_DIR}/libSmilesParse.so ${RDKIT_LIB_DIR}/libDataStructs.so ${RDKIT_LIB_DIR}/libFingerprints.so ${RDKIT_LIB_DIR}/libSubgraphs.so ${RDKIT_LIB_DIR}/libDistGeomHelpers.so ${RDKIT_LIB_DIR}/libForceField.so ${RDKIT_LIB_DIR}/libDepictor.so ${RDKIT_LIB_DIR}/libDistGeometry.so ${RDKIT_LIB_DIR}/libEigenSolvers.so ${RDKIT_LIB_DIR}/libAlignment.so ${RDKIT_LIB_DIR}/libForceFieldHelpers.so ${RDKIT_LIB_DIR}/libFileParsers.so ${RDKIT_LIB_DIR}/libSubstructMatch.so ${RDKIT_LIB_DIR}/libPartialCharges.so  ${RDKIT_LIB_DIR}/libMolDraw2D.so ${RDKIT_LIB_DIR}/libChemTransforms.so  ${RDKIT_LIB_DIR}/libAlignment.so ${RDKIT_LIB_DIR}/libMolTransforms.so  ${RDKIT_LIB_DIR}/libMolAlign.so --memory-init-file 0 -s DISABLE_EXCEPTION_CATCHING=0 -s TOTAL_MEMORY=16777216 -s ALLOW_MEMORY_GROWTH=1)
	    add_custom_target(rdkitjs ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/rdkit.js)
	endif()

	# RDkit new lib format since 2017.03.1
	if(RDKitNEWLIB)
    	# ASM compiler
		# RDkit new lib format since 2017.03.1
		if(ASM)
		    add_custom_command(
		        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/rdkit.js
		        COMMAND ${EMSCRIPTEN_BIN}/em++ --bind -o rdkit.js -s MODULARIZE=1  -s EXPORT_NAME='"'rdk'"'  ${RDKITJS_SRC}  ${CMAKE_CXX_FLAGS}  -I${RDKIT_INCLUDE_DIR} -I${Boost_INCLUDE_DIR} ${Boost_LIBRARIES} ${RDKIT_LIB_DIR}/libRDKitGraphMol.so ${RDKIT_LIB_DIR}/libRDKitDescriptors.so ${RDKIT_LIB_DIR}/libRDKitRDGeneral.so ${RDKIT_LIB_DIR}/libRDKitRDGeometryLib.so ${RDKIT_LIB_DIR}/libRDKitSmilesParse.so ${RDKIT_LIB_DIR}/libRDKitDataStructs.so ${RDKIT_LIB_DIR}/libRDKitFingerprints.so ${RDKIT_LIB_DIR}/libRDKitSubgraphs.so ${RDKIT_LIB_DIR}/libRDKitDistGeomHelpers.so ${RDKIT_LIB_DIR}/libRDKitForceField.so ${RDKIT_LIB_DIR}/libRDKitDepictor.so ${RDKIT_LIB_DIR}/libRDKitDistGeometry.so ${RDKIT_LIB_DIR}/libRDKitEigenSolvers.so ${RDKIT_LIB_DIR}/libRDKitAlignment.so ${RDKIT_LIB_DIR}/libRDKitForceFieldHelpers.so ${RDKIT_LIB_DIR}/libRDKitFileParsers.so ${RDKIT_LIB_DIR}/libRDKitSubstructMatch.so ${RDKIT_LIB_DIR}/libRDKitPartialCharges.so  ${RDKIT_LIB_DIR}/libRDKitMolDraw2D.so ${RDKIT_LIB_DIR}/libRDKitChemTransforms.so ${RDKIT_LIB_DIR}/libRDKitMolTransforms.so  ${RDKIT_LIB_DIR}/libRDKitMolAlign.so --memory-init-file 0 -s DISABLE_EXCEPTION_CATCHING=0 -s TOTAL_MEMORY=16777216 -s ALLOW_MEMORY_GROWTH=1)
		endif()

    	# WASM compiler

		if(WASM)
		    add_custom_command(
		        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/rdkit.js
		        COMMAND ${EMSCRIPTEN_BIN}/em++ --bind -o rdkit.js -s WASM=1 -s MODULARIZE=1  -s EXPORT_NAME='"'rdk'"'  ${RDKITJS_SRC}  ${CMAKE_CXX_FLAGS}  -I${RDKIT_INCLUDE_DIR} -I${Boost_INCLUDE_DIR} ${Boost_LIBRARIES} ${RDKIT_LIB_DIR}/libRDKitGraphMol.so ${RDKIT_LIB_DIR}/libRDKitDescriptors.so ${RDKIT_LIB_DIR}/libRDKitRDGeneral.so ${RDKIT_LIB_DIR}/libRDKitRDGeometryLib.so ${RDKIT_LIB_DIR}/libRDKitSmilesParse.so ${RDKIT_LIB_DIR}/libRDKitDataStructs.so ${RDKIT_LIB_DIR}/libRDKitFingerprints.so ${RDKIT_LIB_DIR}/libRDKitSubgraphs.so ${RDKIT_LIB_DIR}/libRDKitDistGeomHelpers.so ${RDKIT_LIB_DIR}/libRDKitForceField.so ${RDKIT_LIB_DIR}/libRDKitDepictor.so ${RDKIT_LIB_DIR}/libRDKitDistGeometry.so ${RDKIT_LIB_DIR}/libRDKitEigenSolvers.so ${RDKIT_LIB_DIR}/libRDKitAlignment.so ${RDKIT_LIB_DIR}/libRDKitForceFieldHelpers.so ${RDKIT_LIB_DIR}/libRDKitFileParsers.so ${RDKIT_LIB_DIR}/libRDKitSubstructMatch.so ${RDKIT_LIB_DIR}/libRDKitPartialCharges.so  ${RDKIT_LIB_DIR}/libRDKitMolDraw2D.so ${RDKIT_LIB_DIR}/libRDKitChemTransforms.so ${RDKIT_LIB_DIR}/libRDKitMolTransforms.so  ${RDKIT_LIB_DIR}/libRDKitMolAlign.so --memory-init-file 0 -s DISABLE_EXCEPTION_CATCHING=0 -s TOTAL_MEMORY=16777216 -s ALLOW_MEMORY_GROWTH=1)
		endif()

	    add_custom_target(rdkitjs ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/rdkit.js)
    endif()
endif()
