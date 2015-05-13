
#include <emscripten/bind.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/BitOps.h>

using namespace emscripten;
using RDKit::ROMol;
using RDKit::RWMol;

class Molecule {
    
public:

  Molecule(RWMol* mol): rdmol(mol) {};
  
  unsigned int getNumAtoms() {
      return rdmol->getNumAtoms();
  };
  
  
   std::string getFP()
    {
        ExplicitBitVect *finger =  RDKit::RDKFingerprintMol(*rdmol);
        return BitVectToText(*finger);      
   };
  
  
  
  
  int getMW() {
      return RDKit::Descriptors::calcAMW(*rdmol);
  };
  
  static Molecule *fromSmiles(std::string smiles) {
      rdErrorLog->df_enabled = false;
      return new Molecule(RDKit::SmilesToMol(smiles));
  };

private:
  RWMol* rdmol;
  
};

// Binding code
EMSCRIPTEN_BINDINGS(rdmol) {
  class_<Molecule>("Molecule")
    .function("getNumAtoms", &Molecule::getNumAtoms, allow_raw_pointers())
    .function("getMW", &Molecule::getMW, allow_raw_pointers())
    .function("getFP", &Molecule::getFP, allow_raw_pointers())
    .class_function("fromSmiles", &Molecule::fromSmiles, allow_raw_pointers());
}



// /emscripten/emscripten/em++  --bind -o rdmol.js ../rdmol.cpp -Icode -Iinclude lib/libGraphMol.so lib/libDescriptors.so lib/libRDGeneral.so lib/libRDGeometryLib.so lib/libSmilesParse.so lib/libDataStructs.so lib/libFingerprints.so lib/libSubgraphs.so  -O2 


