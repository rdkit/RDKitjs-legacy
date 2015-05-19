#include "rdmol.h"

//standard libraries
#include <string>
#include <vector>
#include <utility>  // std::pair, std::get

// Boost include for rdkit dependence
#include <boost/cstdint.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>


/* RDkit libraries */
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/MolOps.h>

// stream Mol2File
#include <RDGeneral/StreamOps.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

//Drawing
//#include <GraphMol/MolDrawing/MolDrawing.h>
//#include <GraphMol/MolDrawing/DrawingToSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>

// 2D
#include <GraphMol/Depictor/RDDepictor.h>

// cpickle
#include <GraphMol/MolPickler.h>


// murko
#include <GraphMol/ChemTransforms/ChemTransforms.h>


#include <GraphMol/DistGeomHelpers/Embedder.h>
// comments thegodone & Paolo => MMFF.h and Builder.h need to be patch to avoid class issues! 13_05_2015
// #include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
// #include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>

#include <GraphMol/FileParsers/MolWriters.h>

#include <ForceField/MMFF/Params.h>

#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Substruct/SubstructUtils.h>


using namespace std;
using RDKit::ROMol;
using RDKit::RWMol;

Molecule::Molecule(RWMol* mol): rdmol(mol), rdquery(0) {
}

Molecule::~Molecule() {
  if(rdmol != 0)
    delete rdmol;

  if(rdquery != 0)
    delete rdquery;
}

Molecule* Molecule::fromSmiles(std::string smiles) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::SmilesToMol(smiles));
}

Molecule* Molecule::Mol2BlockToMol(std::string molBlock) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::Mol2BlockToMol(molBlock,true,true));
}

Molecule* Molecule::MolBlockToMol(std::string molBlock) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::MolBlockToMol(molBlock, true, true));
}

Molecule* Molecule::fromSmarts(std::string smarts) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::SmartsToMol(smarts));
}

unsigned int Molecule::getNumAtoms()
{
    return rdmol->getNumAtoms();
}

void Molecule::MolToBinary()
{
    std::string res;
    RDKit::MolPickler::pickleMol(*rdmol,res);
}


std::string Molecule::getFP()
{
    ExplicitBitVect* finger =  RDKit::RDKFingerprintMol(*rdmol);
    return BitVectToText(*finger);
}


std::string Molecule::getMorganFP2()
{
    ExplicitBitVect* finger =  RDKit::MorganFingerprints::getFingerprintAsBitVect(*rdmol,2,2048);
    return BitVectToText(*finger);
}


std::string Molecule::getMorganFP3()
{
    ExplicitBitVect* finger =  RDKit::MorganFingerprints::getFingerprintAsBitVect(*rdmol,3,2048);
    return BitVectToText(*finger);
}

std::vector<double> Molecule::MMFFoptimizeMolecule()
{
    std::vector<double> res(2);
    /*std::pair<int, double> p = RDKit::MMFF::MMFFOptimizeMolecule(*rdmol);
    res[0] = static_cast<double>(p.first);
    res[1] = p.second;*/
    return res;
}

std::vector<double> Molecule::UFFOptimizeMolecule()
{
    std::vector<double> res(2);
    std::pair<int, double> p = RDKit::UFF::UFFOptimizeMolecule(*rdmol);
    res[0] = static_cast<double>(p.first);
    res[1] = p.second;
    return res;
}

void Molecule::Murcko()
{
    RDKit::MurckoDecompose(*rdmol);
}

std::vector<std::string> Molecule::getproplist()
{
    return rdmol->getPropList();
}

std::string Molecule::smilewrite()
{
    std::stringstream ss;
    RDKit::SmilesWriter *writer = new RDKit::SmilesWriter(&ss, " ","Name",false);
    writer->write(*rdmol);
    writer->flush();
    return ss.str();
}

std::string Molecule::sdwrite()
{
    std::stringstream ss;
    RDKit::SDWriter *writer = new RDKit::SDWriter(&ss,false);
    writer->write(*rdmol);
    writer->flush();
    return ss.str();
}

unsigned int Molecule::compute2DCoords()
{
    return RDDepict::compute2DCoords(*rdmol);
}

std::string Molecule::Drawing2D()
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


int Molecule::Embedmolecule3D()
{
    return RDKit::DGeomHelpers::EmbedMolecule(*rdmol);
}

void Molecule::addHs()
{
    return RDKit::MolOps::addHs(*rdmol);
}


void Molecule::removeHs()
{
    return RDKit::MolOps::removeHs(*rdmol);
}


void Molecule::sanitizeMol()
{
    return RDKit::MolOps::sanitizeMol(*rdmol);
}

void Molecule::cleanUp()
{
    return RDKit::MolOps::cleanUp(*rdmol);
}

void Molecule::Kekulize()
{
    return RDKit::MolOps::Kekulize(*rdmol);
}


int Molecule::getMW()
{
    return RDKit::Descriptors::calcAMW(*rdmol);
}


double Molecule::Molecule::ExactMW()
{
    return RDKit::Descriptors::calcExactMW(*rdmol);
}


std::string Molecule::Formula()
{
    return RDKit::Descriptors::calcMolFormula(*rdmol);
}


double Molecule::Chi0v()
{
    return   RDKit::Descriptors::calcChi0v(*rdmol);
}


double Molecule::Chi1v()
{
    return   RDKit::Descriptors::calcChi1v (*rdmol);
}

double Molecule::Chi2v()
{
    return   RDKit::Descriptors::calcChi2v (*rdmol);
}


double Molecule::Chi3v()
{
    return   RDKit::Descriptors::calcChi3v (*rdmol);
}



double Molecule::Chi4v()
{
    return   RDKit::Descriptors::calcChi4v (*rdmol);
}


double Molecule::Chi0n()
{
    return   RDKit::Descriptors::calcChi0n (*rdmol);
}


double Molecule::Chi1n()
{
    return   RDKit::Descriptors::calcChi1n (*rdmol);
}


double Molecule::Chi2n()
{
    return   RDKit::Descriptors::calcChi2n (*rdmol);
}


double Molecule::Chi3n()
{
    return   RDKit::Descriptors::calcChi3n (*rdmol);
}

double Molecule::Chi4n()
{
    return   RDKit::Descriptors::calcChi4n (*rdmol);
}


double Molecule::HallKierAlpha()
{
    return RDKit::Descriptors::calcHallKierAlpha (*rdmol);
}

double Molecule::Kappa1()
{
    return    RDKit::Descriptors::calcKappa1 (*rdmol);
}

double Molecule::Kappa2()
{
    return    RDKit::Descriptors::calcKappa2 (*rdmol);
}

double Molecule::Kappa3()
{
    return    RDKit::Descriptors::calcKappa3 (*rdmol);
}

void Molecule::logp_mr()
{
    
    double logp;
    double mr;
    
    return RDKit::Descriptors::calcCrippenDescriptors (*rdmol,logp,mr);
}


unsigned int Molecule::LipinskiHBA()
{
    return RDKit::Descriptors::calcLipinskiHBA (*rdmol);
}

unsigned int Molecule::LipinskiHBD()
{
    return RDKit::Descriptors::calcLipinskiHBD (*rdmol);
}

unsigned int Molecule::NumRotatableBonds()
{
    return RDKit::Descriptors::calcNumRotatableBonds (*rdmol);
}

unsigned int Molecule::NumHBD()
{
    return RDKit::Descriptors::calcNumHBD (*rdmol);
}

unsigned int Molecule::NumHBA()
{
    return RDKit::Descriptors::calcNumHBA (*rdmol);
}

unsigned int Molecule::NumHeteroatoms()
{
    return RDKit::Descriptors::calcNumHeteroatoms (*rdmol);
}

unsigned int Molecule::NumAmideBonds()
{
    return RDKit::Descriptors::calcNumAmideBonds (*rdmol);
}

double Molecule::FractionCSP3()
{
    return RDKit::Descriptors::calcFractionCSP3 (*rdmol);
}

unsigned int Molecule::NumRings()
{
    return RDKit::Descriptors::calcNumRings (*rdmol);
}

unsigned int Molecule::NumAromaticRings()
{
    return   RDKit::Descriptors::calcNumAromaticRings (*rdmol);
}

unsigned int Molecule::NumAliphaticRings()
{
    return   RDKit::Descriptors::calcNumAliphaticRings (*rdmol);
}

unsigned int Molecule::NumSaturatedRings ()
{
    return  RDKit::Descriptors::calcNumSaturatedRings (*rdmol);
}

unsigned int Molecule::NumHeterocycles()
{
    return   RDKit::Descriptors::calcNumHeterocycles (*rdmol);
}

unsigned int Molecule::NumAromaticHeterocycles()
{
    return    RDKit::Descriptors::calcNumAromaticHeterocycles (*rdmol);
}

unsigned int Molecule::NumAromaticCarbocycles ()
{
    return RDKit::Descriptors::calcNumAromaticCarbocycles (*rdmol);
}


unsigned int Molecule::NumSaturatedHeterocycles()
{
    return   RDKit::Descriptors::calcNumSaturatedHeterocycles (*rdmol);
}

unsigned int Molecule::NumSaturatedCarbocycles()
{
    return RDKit::Descriptors::calcNumSaturatedCarbocycles (*rdmol);
}

unsigned int Molecule::NumAliphaticHeterocycles()
{
    return RDKit::Descriptors::calcNumAliphaticHeterocycles (*rdmol);
}

unsigned int Molecule::NumAliphaticCarbocycles()
{
    return RDKit::Descriptors::calcNumAliphaticCarbocycles (*rdmol);
}

double Molecule::LabuteASA()
{
    return RDKit::Descriptors::calcLabuteASA (*rdmol);
}

double Molecule::TPSA()
{
    return RDKit::Descriptors::calcTPSA (*rdmol);
}

std::vector< double > Molecule::SlogP_VSA()
{
    return RDKit::Descriptors::calcSlogP_VSA (*rdmol);
}


std::vector< double > Molecule::SMR_VSA() {
    return RDKit::Descriptors::calcSMR_VSA (*rdmol);
}


std::vector< double > Molecule::PEO_VSA()
{
    return RDKit::Descriptors::calcPEOE_VSA (*rdmol);
}


std::vector< unsigned int > Molecule::MQNs()
{
    return RDKit::Descriptors::calcMQNs (*rdmol);
}

std::string Molecule::GetSubstructMatches(std::string smilesref)
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


bool Molecule::HasSubstructMatchStr(std::string smilesref)
{
    rdErrorLog->df_enabled = false;
    rdquery = RDKit::SmartsToMol(smilesref);
    
    RDKit::MatchVectType res;
    return RDKit::SubstructMatch(*rdmol,*rdquery,res);
}



/// get & set & has properties


std::string Molecule::getProp(std::string key) {
    std::string res;
    rdmol->getProp(key,res);
    return res;
}


int Molecule::setProp(std::string key, std::string value) {
    rdmol->setProp(key,value);
    return 0;
}

bool Molecule::hasProp(std::string key) {
    return rdmol->hasProp(key);
}
