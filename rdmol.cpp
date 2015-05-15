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

//Drawing
//#include <GraphMol/MolDrawing/MolDrawing.h>
//#include <GraphMol/MolDrawing/DrawingToSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

// 2D
#include <GraphMol/Depictor/RDDepictor.h>

// cpickle
#include <GraphMol/MolPickler.h>


// murko
#include <GraphMol/ChemTransforms/ChemTransforms.h>


#include <GraphMol/DistGeomHelpers/Embedder.h>
// comments thegodone & Paolo => MMFF.h and Builder.h need to be patch to avoid class issues! 13_05_2015
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <utility>      // std::pair, std::get


#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>



#include <GraphMol/FileParsers/MolWriters.h>

#include <ForceField/MMFF/Params.h>

#include <boost/cstdint.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <string>


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
    
    
    void MolToBinary()
    {
        std::string res;
        RDKit::MolPickler::pickleMol(*rdmol,res);
    }
    
    
    
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
    
    
    
    std::vector<double> MMFFoptimizeMolecule()
    {
        std::vector<double> res(2);
        std::pair<int, double> p = RDKit::MMFF::MMFFOptimizeMolecule(*rdmol);
        res[0] = static_cast<double>(p.first);
        res[1] = p.second;
        return res;
    }
    
   
    
    
    std::vector<double> UFFOptimizeMolecule()
    {
        std::vector<double> res(2);
        std::pair<int, double> p = RDKit::UFF::UFFOptimizeMolecule(*rdmol);
        res[0] = static_cast<double>(p.first);
        res[1] = p.second;
        return res;
    }
    
    
    void Murcko()
    {
        RDKit::MurckoDecompose(*rdmol);
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
    
    
    unsigned int compute2DCoords()
    {
        return RDDepict::compute2DCoords(*rdmol);
    }
    
    
    std::string Drawing2D()
    {
        
        
        RDDepict::compute2DCoords(*rdmol);
        WedgeMolBonds(*rdmol,&(rdmol->getConformer()));
        RDKit::MolDraw2DSVG drawer(300,300);
        drawer.drawMolecule(*rdmol);
        drawer.finishDrawing();
        return drawer.getDrawingText();

        //return  RDKit::Drawing::MolToDrawing(*mol);
    }
    //svg = RDKit::Drawing::DrawingToSVG(drawing, 4);

    
    

    int Embedmolecule3D()
    {
        return RDKit::DGeomHelpers::EmbedMolecule(*rdmol);
        
    }
    
    void addHs()
    {
        return RDKit::MolOps::addHs(*rdmol);
    }
    
    
    void removeHs()
    {
        return RDKit::MolOps::removeHs(*rdmol);
    }
    
    
    void sanitizeMol()
    {
        return RDKit::MolOps::sanitizeMol(*rdmol);
    }
    
    void cleanUp()
    {
        return RDKit::MolOps::cleanUp(*rdmol);
    }
    
    void Kekulize()
    {
        return RDKit::MolOps::Kekulize(*rdmol);
    }
    
    
    int getMW()
    {
        return RDKit::Descriptors::calcAMW(*rdmol);
    };
    
    
    double ExactMW()
    {
        return RDKit::Descriptors::calcExactMW(*rdmol);
    }
    
    
    std::string Formula()
    {
        return RDKit::Descriptors::calcMolFormula(*rdmol);
    }
    
    
    double Chi0v()
    {
        return   RDKit::Descriptors::calcChi0v(*rdmol);
    }
    
    
    double Chi1v()
    {
        return   RDKit::Descriptors::calcChi1v (*rdmol);
    }
    
    double Chi2v()
    {
        return   RDKit::Descriptors::calcChi2v (*rdmol);
    }
    
    
    double Chi3v()
    {
        return   RDKit::Descriptors::calcChi3v (*rdmol);
    }
    
    
    
    double Chi4v()
    {
        return   RDKit::Descriptors::calcChi4v (*rdmol);
    }
    
    
    double Chi0n()
    {
        return   RDKit::Descriptors::calcChi0n (*rdmol);
    }
    
    
    double Chi1n()
    {
        return   RDKit::Descriptors::calcChi1n (*rdmol);
    }
    
    
    double Chi2n()
    {
        return   RDKit::Descriptors::calcChi2n (*rdmol);
    }
    
    
    double Chi3n()
    {
        return   RDKit::Descriptors::calcChi3n (*rdmol);
    }
    
    double Chi4n()
    {
        return   RDKit::Descriptors::calcChi4n (*rdmol);
    }
    
    
    double HallKierAlpha()
    {
        return RDKit::Descriptors::calcHallKierAlpha (*rdmol);
    }
    
    double Kappa1()
    {
        return    RDKit::Descriptors::calcKappa1 (*rdmol);
    }
    
    double Kappa2()
    {
        return    RDKit::Descriptors::calcKappa2 (*rdmol);
    }
    
    double Kappa3()
    {
        return    RDKit::Descriptors::calcKappa3 (*rdmol);
    }
    
    void logp_mr()
    {
        
        double logp;
        double mr;
        
        return RDKit::Descriptors::calcCrippenDescriptors (*rdmol,logp,mr);
    }
    
    
    unsigned int LipinskiHBA()
    {
        return RDKit::Descriptors::calcLipinskiHBA (*rdmol);
    }
    
    unsigned int LipinskiHBD()
    {
        return RDKit::Descriptors::calcLipinskiHBD (*rdmol);
    }
    
    unsigned int NumRotatableBonds()
    {
        return RDKit::Descriptors::calcNumRotatableBonds (*rdmol);
    }
    
    unsigned int NumHBD()
    {
        return RDKit::Descriptors::calcNumHBD (*rdmol);
    }
    
    unsigned int NumHBA()
    {
        return RDKit::Descriptors::calcNumHBA (*rdmol);
    }
    
    unsigned int NumHeteroatoms()
    {
        return RDKit::Descriptors::calcNumHeteroatoms (*rdmol);
    }
    
    unsigned int NumAmideBonds()
    {
        return RDKit::Descriptors::calcNumAmideBonds (*rdmol);
    }
    
    double FractionCSP3()
    {
        return RDKit::Descriptors::calcFractionCSP3 (*rdmol);
    }
    
    unsigned int NumRings()
    {
        return RDKit::Descriptors::calcNumRings (*rdmol);
    }
    
    unsigned int NumAromaticRings()
    {
        return   RDKit::Descriptors::calcNumAromaticRings (*rdmol);
    }
    
    unsigned int NumAliphaticRings()
    {
        return   RDKit::Descriptors::calcNumAliphaticRings (*rdmol);
    }
    
    unsigned int NumSaturatedRings ()
    {
        return  RDKit::Descriptors::calcNumSaturatedRings (*rdmol);
    }
    
    unsigned int NumHeterocycles()
    {
        return   RDKit::Descriptors::calcNumHeterocycles (*rdmol);
    }
    
    unsigned int NumAromaticHeterocycles()
    {
        return    RDKit::Descriptors::calcNumAromaticHeterocycles (*rdmol);
    }
    
    unsigned int NumAromaticCarbocycles ()
    {
        return RDKit::Descriptors::calcNumAromaticCarbocycles (*rdmol);
    }
    
    
    unsigned int NumSaturatedHeterocycles()
    {
        return   RDKit::Descriptors::calcNumSaturatedHeterocycles (*rdmol);
    }
    
    unsigned int NumSaturatedCarbocycles()
    {
        return RDKit::Descriptors::calcNumSaturatedCarbocycles (*rdmol);
    }
    
    unsigned int NumAliphaticHeterocycles()
    {
        return RDKit::Descriptors::calcNumAliphaticHeterocycles (*rdmol);
    }
    unsigned int NumAliphaticCarbocycles()
    {
        return RDKit::Descriptors::calcNumAliphaticCarbocycles (*rdmol);
    }
    double LabuteASA()
    {
        return RDKit::Descriptors::calcLabuteASA (*rdmol);
    }
    
    double TPSA()
    {
        return RDKit::Descriptors::calcTPSA (*rdmol);
    }
    
    std::vector< double > SlogP_VSA()
    {
        return RDKit::Descriptors::calcSlogP_VSA (*rdmol);
        
    }
    
    
    std::vector< double > SMR_VSA() {
        return RDKit::Descriptors::calcSMR_VSA (*rdmol);
    }
    
    
    std::vector< double > PEO_VSA()
    {
        return RDKit::Descriptors::calcPEOE_VSA (*rdmol);
        
    }
    
    
    std::vector< unsigned int > MQNs()
    {
        return RDKit::Descriptors::calcMQNs (*rdmol);
    }
    
    
    
    
    
    std::string GetSubstructMatches(std::string smilesref)
    {
        RDKit::MatchVectType matchV;
        std::vector< RDKit::MatchVectType > matches;
        //RWMol* rdquery = fromSmarts(smilesref);
        
        rdErrorLog->df_enabled = false;
        rdquery = RDKit::SmartsToMol(smilesref);
        
        
        
        int matched = RDKit::SubstructMatch(*rdmol,*rdquery,matches,true);
        std::string res = "";
        for(int idx=0;idx<matched;idx++){
            res +=".";
        }
        return res;
        
    }
    
    
    bool HasSubstructMatchStr(std::string smilesref)
    {
        rdErrorLog->df_enabled = false;
        rdquery = RDKit::SmartsToMol(smilesref);
       
        RDKit::MatchVectType res;
        
        return RDKit::SubstructMatch(*rdmol,*rdquery,res);
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
    .class_function("fromSmiles", &Molecule::fromSmiles, allow_raw_pointers())
    .class_function("fromSmarts", &Molecule::fromSmarts, allow_raw_pointers());
    // register the vectors
    register_vector<std::string>("VectorString");
    register_vector<double>("VectorDouble");
    register_vector<unsigned int>("VectorUint");
    
}
