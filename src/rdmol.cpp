#include "rdmol.h"

//standard libraries
#include <string>
#include <vector>
#include <utility>  // pair, get

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
#include <GraphMol/Conformer.h>


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
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
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

Molecule::Molecule(RWMol* mol): rdmol(mol), rdquery(mol) {


}

Molecule::~Molecule() {
  if(rdmol != 0)
    delete rdmol;

  if(rdquery != 0)
    delete rdquery;
}

Molecule* Molecule::fromSmiles(string smiles) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::SmilesToMol(smiles));
}

Molecule* Molecule::Mol2BlockToMol(string molBlock) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::Mol2BlockToMol(molBlock,true,true));
}

Molecule* Molecule::MolBlockToMol(string molBlock) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::MolBlockToMol(molBlock, true, true));
}

Molecule* Molecule::fromSmarts(string smarts) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::SmartsToMol(smarts));
}

Molecule* Molecule::molFromPickle(string pickle) {
  rdErrorLog->df_enabled = false;

  ROMol* mol =  new ROMol;
  RDKit::MolPickler::molFromPickle(pickle, mol);

  //RWMol* rdquery =  new RWMol();
  //RWMol* rdquery = RWMol(*RDKit::ROMol::ROMol(pickle));
  //RDKit::MolPickler::molFromPickle(pickle,*rdquery);
  return new Molecule(dynamic_cast<RWMol *>(mol));
}


Molecule* Molecule::MurckofromSmiles(string smi)
{

//     ROMol *mol=static_cast<RWMol &>(RDKit::SmilesToMol(smi));
     ROMol *mol=RDKit::SmilesToMol(smi);
     return new Molecule(dynamic_cast<RWMol *>(RDKit::MurckoDecompose(*mol)));

    //return RDKit::MolToSmiles(*nmol);


}

Molecule* Molecule::newmolecule()
{

     RWMol *mol= new RWMol();
     return new Molecule(mol);



}




unsigned int Molecule::addAtom (int atomid)
{
   RDKit::Atom* atom= new RDKit::Atom(atomid);  
   return rdmol->addAtom(atom);

}       


unsigned int Molecule::addBond (unsigned int beginAtomIdx, unsigned int endAtomIdx,RDKit::Bond::BondType bondtype)
{
   return rdmol->addBond(beginAtomIdx,endAtomIdx,bondtype);

}       


void Molecule::setBondDir (int Bondid, RDKit::Bond::BondDir bonddir)
{
   rdmol->getBondWithIdx(Bondid)->setBondDir(bonddir);

}       



unsigned int Molecule::getNumAtoms()
{
    return rdmol->getNumAtoms();
}

string Molecule::MolToBinary()
{
    string res;
    RDKit::MolPickler::pickleMol(*rdmol,res);
    return res;
}







string Molecule::getFP()
{
    ExplicitBitVect* finger =  RDKit::RDKFingerprintMol(*rdmol);
    return BitVectToText(*finger);
}


string Molecule::getMorganFP2()
{
    ExplicitBitVect* finger =  RDKit::MorganFingerprints::getFingerprintAsBitVect(*rdmol,2,2048);
    return BitVectToText(*finger);
}


string Molecule::getMorganFP3()
{
    ExplicitBitVect* finger =  RDKit::MorganFingerprints::getFingerprintAsBitVect(*rdmol,3,2048);
    return BitVectToText(*finger);
}

vector<double> Molecule::MMFFoptimizeMolecule()
{
    vector<double> res(2);
    pair<int, double> p = RDKit::MMFF::MMFFOptimizeMolecule(*rdmol);
    res[0] = static_cast<double>(p.first);
    res[1] = p.second;
    return res;
}


vector<double> Molecule::MMFFoptimizeMolecule(int maxIters, string mmffVariant)
{
    pair<int, double> p = RDKit::MMFF::MMFFOptimizeMolecule(*rdmol,maxIters,mmffVariant);
    vector<double> res(2);
    res[0] = static_cast<double>(p.first);
    res[1] = p.second;
    return res;
}



vector<double> Molecule::MMFFOptimizeMoleculeConfs (unsigned int numThreads,int  maxIters,string mmffVariant)
{
   vector<pair< int, double>> p;
   RDKit::MMFF::MMFFOptimizeMoleculeConfs(*rdmol,p,numThreads,maxIters,mmffVariant);
   vector<double> res(2*p.size());
   for (int i = 0; i < 2*p.size(); i++) {
        if ( i % 2== 0 )
            res[i] = static_cast<double>(p[i].first);
        else
            res[p.size()+i] = p[i].second;

    }

    return res;
}






vector<double> Molecule::UFFOptimizeMolecule()
{
    vector<double> res(2);
    pair<int, double> p = RDKit::UFF::UFFOptimizeMolecule(*rdmol);
    res[0] = static_cast<double>(p.first);
    res[1] = p.second;
    return res;
}

void Molecule::Murcko()
{
    RDKit::MurckoDecompose(*rdmol);
}

vector<string> Molecule::getproplist()
{
    return rdmol->getPropList();
}

string Molecule::smilewrite()
{
    stringstream ss;
    RDKit::SmilesWriter *writer = new RDKit::SmilesWriter(&ss, " ","Name",false);
    writer->write(*rdmol);
    writer->flush();
    return ss.str();
}

string Molecule::sdwrite()
{
    stringstream ss;
    RDKit::SDWriter *writer = new RDKit::SDWriter(&ss,false);
    writer->write(*rdmol);
    writer->flush();
    return ss.str();
}




string Molecule::sdwriteConfs()
{
    stringstream ss;
    RDKit::SDWriter *writer = new RDKit::SDWriter(&ss,false);
    
    for(int i=0; i<rdmol->getNumConformers(); ++i){
         writer->write(*rdmol,i);
     }

   /* for (int i = 0; i<confsize ;++i)
    {
        writer->write(*rdmol,i);
    }
*/
    //writer->flush();
    return ss.str();
}


unsigned int Molecule::compute2DCoords()
{
    return RDDepict::compute2DCoords(*rdmol);
}

string Molecule::Drawing2D()
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


int Molecule::EmbedMolecule(unsigned int maxIterations, int seed)

{   
 
    return RDKit::DGeomHelpers::EmbedMolecule(*rdmol,maxIterations,seed);
}


int Molecule::EmbedMolecule()

{   
 
    return RDKit::DGeomHelpers::EmbedMolecule(*rdmol);
}


vector<int> Molecule::EmbedMultipleConfs(unsigned int numConfs, unsigned int maxIterations, int seed)
{   
 
    return RDKit::DGeomHelpers::EmbedMultipleConfs(*rdmol,numConfs,maxIterations,seed);
}


vector<int> Molecule::EmbedMultipleConfs()
{   
 
    return RDKit::DGeomHelpers::EmbedMultipleConfs(*rdmol);
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


string Molecule::Formula()
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

vector< double > Molecule::SlogP_VSA()
{
    return RDKit::Descriptors::calcSlogP_VSA (*rdmol);
}


vector< double > Molecule::SMR_VSA() {
    return RDKit::Descriptors::calcSMR_VSA (*rdmol);
}


vector< double > Molecule::PEO_VSA()
{
    return RDKit::Descriptors::calcPEOE_VSA (*rdmol);
}


vector< unsigned int > Molecule::MQNs()
{
    return RDKit::Descriptors::calcMQNs (*rdmol);
}

string Molecule::GetSubstructMatches(string smilesref)
{
    RDKit::MatchVectType matchV;
    vector< RDKit::MatchVectType > matches;

    rdErrorLog->df_enabled = false;
    rdquery = RDKit::SmartsToMol(smilesref);

    int matched = RDKit::SubstructMatch(*rdmol,*rdquery,matches,true);
    string res = "";
    for(int idx=0;idx<matched;idx++){
        res +=".";
    }
    return res;
}


bool Molecule::HasSubstructMatchStr(string smilesref)
{
    rdErrorLog->df_enabled = false;
    rdquery = RDKit::SmartsToMol(smilesref);
    
    RDKit::MatchVectType res;
    return RDKit::SubstructMatch(*rdmol,*rdquery,res);
}



/// get & set & has properties


string Molecule::getProp(string key) {
    string res;
    rdmol->getProp(key,res);
    return res;
}

int Molecule::getNumConformers() {
    return rdmol->getNumConformers();
}

/*
int Molecule::getConformer() {
    return rdmol->getConformer();
}
*/


int Molecule::setProp(string key, string value) {
    rdmol->setProp(key,value);
    return 0;
}

bool Molecule::hasProp(string key) {
    return rdmol->hasProp(key);
}
