#include <emscripten/bind.h>

#include <string>
#include <vector>

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>


using namespace std;
using namespace emscripten;
using RDKit::ROMol;
using RDKit::RWMol;

class Molecule
{
    public:
      Molecule(RWMol *mol){
        if(rdmol != 0)
            delete rdmol;
        rdmol = mol;
        rdquery = 0;
      }
      ~Molecule(){
            if(rdmol != 0)
              delete rdmol;

            if(rdquery != 0)
              delete rdquery;
       };

      unsigned int getNumAtoms();
      void MolToBinary();
      std::string getFP();
      std::string getMorganFP2();
      std::string getMorganFP3();
      std::vector<double> MMFFoptimizeMolecule();
      std::vector<double> UFFOptimizeMolecule();
      void Murcko();
      std::vector<std::string> getproplist();
      std::string smilewrite();
      std::string sdwrite();
      unsigned int compute2DCoords();
      std::string Drawing2D();
      int Embedmolecule3D();
      void addHs();
      void removeHs();
      void sanitizeMol();
      void cleanUp();
      void Kekulize();
      int getMW();
      double ExactMW();
      std::string Formula();
      double Chi0v();
      double Chi1v();
      double Chi2v();
      double Chi3v();
      double Chi4v();
      double Chi0n();
      double Chi1n();
      double Chi2n();
      double Chi3n();
      double Chi4n();
      double HallKierAlpha();
      double Kappa1();
      double Kappa2();
      double Kappa3();
      void logp_mr();
      unsigned int LipinskiHBA();
      unsigned int LipinskiHBD();
      unsigned int NumRotatableBonds();
      unsigned int NumHBD();
      unsigned int NumHBA();
      unsigned int NumHeteroatoms();
      unsigned int NumAmideBonds();
      double FractionCSP3();
      unsigned int NumRings();
      unsigned int NumAromaticRings();
      unsigned int NumAliphaticRings();
      unsigned int NumSaturatedRings();
      unsigned int NumHeterocycles();
      unsigned int NumAromaticHeterocycles();
      unsigned int NumAromaticCarbocycles ();
      unsigned int NumSaturatedHeterocycles();
      unsigned int NumSaturatedCarbocycles();
      unsigned int NumAliphaticHeterocycles();
      unsigned int NumAliphaticCarbocycles();
      double LabuteASA();

      double TPSA();
      std::vector< double > SlogP_VSA();
      std::vector< double > SMR_VSA();
      std::vector< double > PEO_VSA();
      std::vector< unsigned int > MQNs();
      std::string GetSubstructMatches(std::string smilesref);
      bool HasSubstructMatchStr(std::string smilesref);
      /// get & set & has properties
      std::string getProp(std::string key);
      int setProp(std::string key, std::string value);
      bool hasProp(std::string key);

      static Molecule* fromSmiles(std::string smiles) {
          rdErrorLog->df_enabled = false;
          return new Molecule(RDKit::SmilesToMol(smiles));
      }


      static Molecule* Mol2BlockToMol(std::string molBlock) {
        rdErrorLog->df_enabled = false;
        return new Molecule(RDKit::Mol2BlockToMol(molBlock,true,true));
      }

      static Molecule* MolBlockToMol(std::string molBlock)
      {
        rdErrorLog->df_enabled = false;
        return new Molecule(RDKit::MolBlockToMol(molBlock, true, true));
      }

      static Molecule* fromSmarts(std::string smarts) {
        rdErrorLog->df_enabled = false;
        return new Molecule(RDKit::SmartsToMol(smarts));
      }

    /*
      static Molecule *Mol2FileToMol(std::string fname)
      {
          return new Molecule(Mol2FileToMol(fname));
      }
      */

    private:
        RWMol* rdmol;
        RWMol* rdquery;
};

Molecule* passThrough(Molecule* ptr) { return ptr; }

// Binding code
EMSCRIPTEN_BINDINGS(rdmol) {
    class_<Molecule>("Molecule")
    //
    .function("getNumAtoms", &Molecule::getNumAtoms, allow_raw_pointers())
    
    // fingerprints
    .function("getFP", &Molecule::getFP, allow_raw_pointers())
    .function("getMorganFP2", &Molecule::getMorganFP2, allow_raw_pointers())
    .function("getMorganFP3", &Molecule::getMorganFP3, allow_raw_pointers())
    
    // molops basic functions
    .function("addHs", &Molecule::addHs, allow_raw_pointers())
    .function("removeHs", &Molecule::removeHs, allow_raw_pointers())
    .function("sanitizeMol", &Molecule::sanitizeMol, allow_raw_pointers())
    .function("cleanUp", &Molecule::cleanUp, allow_raw_pointers())
    .function("Kekulize", &Molecule::Kekulize, allow_raw_pointers())

    // 2D & 3D molecules
    .function("compute2DCoords", &Molecule::compute2DCoords, allow_raw_pointers())
    .function("Embedmolecule3D", &Molecule::Embedmolecule3D, allow_raw_pointers())
    .function("MMFFoptimizeMolecule", &Molecule::MMFFoptimizeMolecule, allow_raw_pointers())
    .function("UFFOptimizeMolecule", &Molecule::UFFOptimizeMolecule, allow_raw_pointers())
    
    // drawing molecules
    .function("Drawing2D", &Molecule::Drawing2D, allow_raw_pointers())

    // murcko
    .function("Murcko", &Molecule::Murcko, allow_raw_pointers())
    .function("MolToBinary", &Molecule::MolToBinary, allow_raw_pointers())

    
    
    // writer basic functions
    .function("sdwrite", &Molecule::sdwrite, allow_raw_pointers())
    .function("smilewrite", &Molecule::smilewrite, allow_raw_pointers())
    
    // properties
    .function("getproplist", &Molecule::getproplist, allow_raw_pointers())
    .function("getProp", &Molecule::getProp, allow_raw_pointers())
    .function("setProp", &Molecule::setProp, allow_raw_pointers())
    .function("hasProp", &Molecule::hasProp, allow_raw_pointers())
    
    // susbtructure
    .function("GetSubstructMatches", &Molecule::GetSubstructMatches, allow_raw_pointers())
    .function("HasSubstructMatchStr", &Molecule::HasSubstructMatchStr, allow_raw_pointers())
    
    // descriptors
    .function("getMW", &Molecule::getMW, allow_raw_pointers())
    .function("ExactMW",&Molecule::ExactMW ,allow_raw_pointers())
    .function("Formula",&Molecule::Formula ,allow_raw_pointers())
    .function("Chi0v",&Molecule::Chi0v ,allow_raw_pointers())
    .function("Chi1v",&Molecule::Chi1v ,allow_raw_pointers())
    .function("Chi2v",&Molecule::Chi2v ,allow_raw_pointers())
    .function("Chi3v",&Molecule::Chi3v ,allow_raw_pointers())
    .function("Chi4v",&Molecule::Chi4v ,allow_raw_pointers())
    .function("Chi0n",&Molecule::Chi0n ,allow_raw_pointers())
    .function("Chi1n",&Molecule::Chi1n ,allow_raw_pointers())
    .function("Chi2n",&Molecule::Chi2n ,allow_raw_pointers())
    .function("Chi3n",&Molecule::Chi3n ,allow_raw_pointers())
    .function("Chi4n",&Molecule::Chi4n ,allow_raw_pointers())
    .function("HallKierAlpha",&Molecule::HallKierAlpha ,allow_raw_pointers())
    .function("Kappa1",&Molecule::Kappa1 ,allow_raw_pointers())
    .function("Kappa2",&Molecule::Kappa2 ,allow_raw_pointers())
    .function("Kappa3",&Molecule::Kappa3 ,allow_raw_pointers())
    .function("logp_mr",&Molecule::logp_mr ,allow_raw_pointers())
    .function("LipinskiHBA",&Molecule::LipinskiHBA ,allow_raw_pointers())
    .function("LipinskiHBD",&Molecule::LipinskiHBD ,allow_raw_pointers())
    .function("NumRotatableBonds",&Molecule::NumRotatableBonds ,allow_raw_pointers())
    .function("NumHBD",&Molecule::NumHBD ,allow_raw_pointers())
    .function("NumHBA",&Molecule::NumHBA ,allow_raw_pointers())
    .function("NumHeteroatoms",&Molecule::NumHeteroatoms ,allow_raw_pointers())
    .function("NumAmideBonds",&Molecule::NumAmideBonds ,allow_raw_pointers())
    .function("FractionCSP3",&Molecule::FractionCSP3 ,allow_raw_pointers())
    .function("NumRings",&Molecule::NumRings ,allow_raw_pointers())
    .function("NumAromaticRings",&Molecule::NumAromaticRings ,allow_raw_pointers())
    .function("NumAliphaticRings",&Molecule::NumAliphaticRings ,allow_raw_pointers())
    .function("NumSaturatedRings ",&Molecule::NumSaturatedRings ,allow_raw_pointers())
    .function("NumHeterocycles",&Molecule::NumHeterocycles ,allow_raw_pointers())
    .function("NumAromaticHeterocycles",&Molecule::NumAromaticHeterocycles ,allow_raw_pointers())
    .function("NumAromaticCarbocycles",&Molecule::NumAromaticCarbocycles ,allow_raw_pointers())
    .function("NumSaturatedHeterocycles",&Molecule::NumSaturatedHeterocycles ,allow_raw_pointers())
    .function("NumSaturatedCarbocycles",&Molecule::NumSaturatedCarbocycles ,allow_raw_pointers())
    .function("NumAliphaticHeterocycles",&Molecule::NumAliphaticHeterocycles ,allow_raw_pointers())
    .function("NumAliphaticCarbocycles",&Molecule::NumAliphaticCarbocycles ,allow_raw_pointers())
    .function("LabuteASA",&Molecule::LabuteASA ,allow_raw_pointers())
    .function("TPSA",&Molecule::TPSA ,allow_raw_pointers())
    .function("SlogP_VSA",&Molecule::SlogP_VSA ,allow_raw_pointers())
    .function("SMR_VSA",&Molecule::SMR_VSA ,allow_raw_pointers())
    .function("PEO_VSA",&Molecule::PEO_VSA ,allow_raw_pointers())
    .function("MQNs",&Molecule::MQNs ,allow_raw_pointers())
    
    // create class from smiles or smarts
   // .function("Mol2FileToMol", &Molecule::Mol2FileToMol, allow_raw_pointers())
    .class_function("MolBlockToMol", &Molecule::MolBlockToMol, allow_raw_pointers())
    .class_function("Mol2BlockToMol", &Molecule::Mol2BlockToMol, allow_raw_pointers())
    .class_function("fromSmiles", &Molecule::fromSmiles, allow_raw_pointers())
    .class_function("fromSmarts", &Molecule::fromSmarts, allow_raw_pointers());
    // register the vectors
    register_vector<std::string>("VectorString");
    register_vector<double>("VectorDouble");
    register_vector<unsigned int>("VectorUint");
}
