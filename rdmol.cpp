#include <emscripten/bind.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>

#include <DataStructs/BitOps.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
// comments thegodone & Paolo => MMFF.h and Builder.h need to be patch to avoid class issues! 13_05_2015
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <utility>      // std::pair, std::get


#include <GraphMol/FileParsers/MolWriters.h>

#include <ForceField/MMFF/Params.h>

#include <boost/cstdint.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>

#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Substruct/SubstructUtils.h>




using namespace emscripten;
using RDKit::ROMol;
using RDKit::RWMol;

class Molecule {
    
public:
    
    
    Molecule(RWMol* mol): rdmol(mol), rdquery(mol) {};
    
    unsigned int getNumAtoms() {
        return rdmol->getNumAtoms();
    };
    
    std::string getFP()
    {
        ExplicitBitVect* finger =  RDKit::RDKFingerprintMol(*rdmol);
        return BitVectToText(*finger);
    };
    
    
    std::string getMorganFP2()
    {
        ExplicitBitVect* finger =  RDKit::MorganFingerprints::getFingerprintAsBitVect(*rdmol,2,2048);
        return BitVectToText(*finger);
    };
    
    
    
    std::string getMorganFP3()
    {
        ExplicitBitVect* finger =  RDKit::MorganFingerprints::getFingerprintAsBitVect(*rdmol,3,2048);
        return BitVectToText(*finger);
    };

    
    /*
    std::pair<int, double> MMFFoptimizeMolecule()
    {
        std::pair<int, double> p= RDKit::MMFF::MMFFOptimizeMolecule(*rdmol);
        
       // std::vector<double> g;
       // g.push_back((double)std::get<0>(p));
       // g.push_back(std::get<0>(p));
        
        //return std::get<0>(p);
        return p;
    }
    */
    
    
    std::string MMFFoptimizeMolecule()
    {
        std::string res="";
        std::pair<int, double> p = RDKit::MMFF::MMFFOptimizeMolecule(*rdmol);
        res += boost::lexical_cast<std::string>(p.first);
        res += boost::lexical_cast<std::string>(",");
        res += boost::lexical_cast<std::string>(p.second);
        return res;
    }
    
    
    
    std::vector<std::string> getproplist()
    
    {
        
        return rdmol->getPropList();
    
    }
    
    
    std::string smilewrite()
    {
        std::stringstream ss;
        RDKit::SmilesWriter *writer = new RDKit::SmilesWriter(&ss, " ","Name",false);
        writer->write(*rdmol);
        writer->flush();
        return ss.str();
    }
   
    
    
     std::string sdwrite()
    {
        std::stringstream ss;
        RDKit::SDWriter *writer = new RDKit::SDWriter(&ss,false);
        writer->write(*rdmol);
        writer->flush();
        return ss.str();
    }
   
    
    
    int Embedmolecule3D()
    {
        return RDKit::DGeomHelpers::EmbedMolecule(*rdmol);
        
    }
    
    void addHs()
    {
        return RDKit::MolOps::addHs(*rdmol);
    }
    
    
    int getMW() {
        return RDKit::Descriptors::calcAMW(*rdmol);
    };
    
   
    std::string GetSubstructMatches()
    {
        RDKit::MatchVectType matchV;
        std::vector< RDKit::MatchVectType > matches;
       int matched = RDKit::SubstructMatch(*rdmol,*rdquery,matches,true);
       std::string res = "";
        for(int idx=0;idx<matched;idx++){
            res +=".";
        }
        return res;
    
    }
    
    
    static Molecule *fromSmiles(std::string smiles) {
        rdErrorLog->df_enabled = false;
        return new Molecule(RDKit::SmilesToMol(smiles));
    };
    
    
    static Molecule *fromSmarts(std::string smarts) {
        rdErrorLog->df_enabled = false;
        return new Molecule(RDKit::SmartsToMol(smarts));
    };
    
    
    
private:
    RWMol* rdmol;
    RWMol* rdquery;
    
};



// Binding code
EMSCRIPTEN_BINDINGS(rdmol) {
    class_<Molecule>("Molecule")
    .function("getNumAtoms", &Molecule::getNumAtoms, allow_raw_pointers())
    .function("getMW", &Molecule::getMW, allow_raw_pointers())
    .function("getFP", &Molecule::getFP, allow_raw_pointers())
    .function("getMorganFP2", &Molecule::getMorganFP2, allow_raw_pointers())
    .function("getMorganFP3", &Molecule::getMorganFP3, allow_raw_pointers())
    .function("addHs", &Molecule::addHs, allow_raw_pointers())
    .function("Embedmolecule3D", &Molecule::Embedmolecule3D, allow_raw_pointers())
    .function("MMFFoptimizeMolecule", &Molecule::MMFFoptimizeMolecule, allow_raw_pointers())
    .function("sdwrite", &Molecule::sdwrite, allow_raw_pointers())
    .function("smilewrite", &Molecule::smilewrite, allow_raw_pointers())
    .function("getproplist", &Molecule::getproplist, allow_raw_pointers())
    .function("GetSubstructMatches", &Molecule::GetSubstructMatches, allow_raw_pointers())
    .class_function("fromSmiles", &Molecule::fromSmiles, allow_raw_pointers())
    .class_function("fromSmarts", &Molecule::fromSmarts, allow_raw_pointers());
    register_vector<std::string>("VectorString");

}







// /emscripten/emscripten/em++  --bind -o rdmol.js ../rdmol.cpp -Icode -Iinclude lib/libGraphMol.so lib/libDescriptors.so lib/libRDGeneral.so lib/libRDGeometryLib.so lib/libSmilesParse.so lib/libDataStructs.so lib/libFingerprints.so lib/libSubgraphs.so  -O2