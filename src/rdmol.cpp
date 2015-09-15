#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>


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
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/FileParsers.h>

// fkingreprints
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/MACCS.h>

#include <DataStructs/BitOps.h>
#include <DataStructs/SparseIntVect.h>

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

// bond
#include <GraphMol/Bond.h>

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

#include <boost/tokenizer.hpp>


#include <emscripten/emscripten.h>

#include <GraphMol/molAlign/AlignMolecules.h>

#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/PartialCharges/GasteigerParams.h>
#include <map>
#include <boost/cstdint.hpp>
#include <typeinfo>



typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
//template<typename IndexType>
//typedef std::map<IndexType,int> StorageType;



using namespace std;
using RDKit::ROMol;
using RDKit::RWMol;

/* adding another class... ???
Similarity::Similarity(RWMol* mol1, RWMol* mol2): rdmol1(mol1),rdmol2(mol2) {}


Similarity::~Similarity() {
  if(mol1 != 0)
    delete mol1;
  if(mol2 != 0)
    delete mol2;
}


double  Similarity::Tanimoto ()       
{

   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol1,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*rdmol2,2);

   return TanimotoSimilarity (*v1,*v2);
}
*/


Molecule::Molecule(RWMol* mol): rdmol(mol) {}


Molecule::~Molecule() {
  if(rdmol != 0) {
    //printf("Destroy called\n");
    delete rdmol;
    rdmol =0;
  }
}


Molecule* Molecule::fromSmiles(string smiles) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::SmilesToMol(smiles));
}

Molecule* Molecule::Mol2BlockToMol(string molBlock) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::Mol2BlockToMol(molBlock,true,false));
}

Molecule* Molecule::MolBlockToMol(string molBlock) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::MolBlockToMol(molBlock, true, false));
}


Molecule* Molecule::fromSmarts(string smarts) {
  rdErrorLog->df_enabled = false;
  return new Molecule(RDKit::SmartsToMol(smarts));
}
  Molecule* Molecule::molFromPickle(string pickle) {
  RWMol *mol = new RWMol();
  RDKit::MolPickler::molFromPickle(pickle, *mol);
  return new Molecule(mol);
}


Molecule* Molecule::MurckofromSmiles(string smi)
{
     ROMol *mol=RDKit::SmilesToMol(smi);
     return new Molecule(dynamic_cast<RWMol *>(RDKit::MurckoDecompose(*mol)));
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


//// need to enumerate the BondType & BondDir lists ... fro emscripten ???
unsigned int Molecule::addBond (unsigned int beginAtomIdx, unsigned int endAtomIdx, int bondtypeid)
{  
 
   // RDKit::Bond::BondType realbondtype = (RDKit::Bond::BondType)enum.Parse(typeof(RDKit::Bond::BondType), bondtype);
 RDKit::Bond::BondType castEnum = (RDKit::Bond::BondType)bondtypeid;
   

   return rdmol->addBond(beginAtomIdx,endAtomIdx,castEnum);

}       


void Molecule::setBondDir (int Bondid, int bonddirid)
{

 RDKit::Bond::BondDir castEnum = (RDKit::Bond::BondDir)bonddirid;

   rdmol->getBondWithIdx(Bondid)->setBondDir(castEnum);

}       



double  Molecule::TanimotoSimilarityfromSmile ( string smilesref)       
{
   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return TanimotoSimilarity (*v1,*v2);
}

double  Molecule::DiceSimilarityfromSmile ( string smilesref)       
{

   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return DiceSimilarity (*v1,*v2);
}

double  Molecule::TverskySimilarityfromSmile ( string smilesref, double a, double b)       
{
   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return TverskySimilarity (*v1,*v2,a,b);
}






unsigned int Molecule::getNumAtoms()
{
    return rdmol->getNumAtoms();
}


boost::uint32_t Molecule::getAtomCode(int atomid) {
  return RDKit::AtomPairs::getAtomCode(rdmol->getAtomWithIdx(atomid));
}



vector<int> Molecule::getAtomPairFingerprint() {
   RDKit::SparseIntVect<int> * finger =  RDKit::AtomPairs::getAtomPairFingerprint(*rdmol);
   RDKit::SparseIntVect<int>::StorageType gnze = finger->getNonzeroElements();

    int elementsize=gnze.size();
    vector<int>  result(2*elementsize);
    
    map<int, int>::iterator it;
    
    int idx = 0;
    for (it = gnze.begin(); it != gnze.end(); it++ )
    {
        result[idx]=it->first;
        result[elementsize+idx]=it->second;
        idx=idx+1;
    }
    return result;
}

vector<int> Molecule::getHashedAtomPairFingerprint(int size, int atomid1, int atomid2) {
  RDKit::SparseIntVect<int> * finger =  RDKit::AtomPairs::getHashedAtomPairFingerprint(*rdmol,size,atomid1, atomid2);
  RDKit::SparseIntVect<int>::StorageType gnze = finger->getNonzeroElements();
  
  int elementsize=gnze.size();
    vector<int>  result(2*elementsize);
    
    map<int, int>::iterator it;
    
    int idx = 0;
    for (it = gnze.begin(); it != gnze.end(); it++ )
    {
        result[idx]=it->first;
        result[elementsize+idx]=it->second;
        idx=idx+1;
    }
    return result;
}


string Molecule::getHashedAtomPairFingerprintAsBitVect(int size, int atomid1, int atomid2) {
  ExplicitBitVect * finger = RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect(*rdmol,size,atomid1, atomid2);
  return BitVectToText(*finger);
}



string Molecule::MolToBinary()
{
    string res;
    RDKit::MolPickler::pickleMol(*rdmol,res);
    return res;
}


string Molecule::getRDKFP()
{
    ExplicitBitVect* finger =  RDKit::RDKFingerprintMol(*rdmol);
    return BitVectToText(*finger);
}


string Molecule::getLayeredFP(unsigned int layer,unsigned int sizes,unsigned int lengths)
{
    ExplicitBitVect* finger =  RDKit::LayeredFingerprintMol(*rdmol,layer, sizes,lengths);
    return BitVectToText(*finger);
}






string Molecule::getMACCSFP()
{
    ExplicitBitVect* finger =  RDKit::MACCSFingerprints::getFingerprintAsBitVect(*rdmol);
    return BitVectToText(*finger);
}

string Molecule::getPatternFP()
{
    ExplicitBitVect* finger =  RDKit::PatternFingerprintMol(*rdmol);
    return BitVectToText(*finger);
}



string Molecule::getMorganFP(unsigned int sizes,unsigned int lengths)
{
    ExplicitBitVect* finger =  RDKit::MorganFingerprints::getFingerprintAsBitVect(*rdmol,sizes, lengths);
    return BitVectToText(*finger);
}


vector<int> Molecule::getMorganFP_GetOnBits(unsigned int sizes,unsigned int lengths)
{
    ExplicitBitVect* finger =  RDKit::MorganFingerprints::getFingerprintAsBitVect(*rdmol,sizes, lengths);
    IntVect onBits;
    finger->getOnBits(onBits);
    return onBits;
}



map<boost::uint32_t, int> Molecule::getMorganFP_getNonzeroElements(unsigned int sizes)
{
    RDKit::SparseIntVect<boost::uint32_t> * fp =  RDKit::MorganFingerprints::getFingerprint(*rdmol,sizes);
    RDKit::SparseIntVect<boost::uint32_t>::StorageType gnze = fp->getNonzeroElements();
    //map<boost::uint32_t,int> gnze ;
    
    /*for(RDKit::SparseIntVect<boost::uint32_t>::StorageType::const_iterator iter=nze.begin();
        iter!=nze.end();++iter){
        boost::uint32_t v = iter->first;
        int v2 = iter->second;
    }*/
    //gnze.insert(nze.begin(), nze.end());
    //std::cout << gnze << std::endl;

    /*
    map<boost::uint32_t, int>::iterator it;

    for ( it = gnze.begin(); it != gnze.end(); it++ )
    {
        std::cout << it->first  // string (key)
                  << ':'
                  << it->second   // string's value 
                  << std::endl ;
    }
    */
    return gnze;
}




vector<boost::uint32_t> Molecule::getMorganFPlist(unsigned int sizes)
{
    RDKit::SparseIntVect<boost::uint32_t> * fp =  RDKit::MorganFingerprints::getFingerprint(*rdmol,sizes);
    RDKit::SparseIntVect<boost::uint32_t>::StorageType gnze = fp->getNonzeroElements();

    int elementsize=gnze.size();
    vector<boost::uint32_t>  result(2*elementsize);
    
    map<boost::uint32_t, int>::iterator it;
    
    int idx = 0;
    for (it = gnze.begin(); it != gnze.end(); it++ )
    {
        result[idx]=it->first;
        result[elementsize+idx]=it->second;
        idx=idx+1;
    }
    return result;
}


vector<string> Molecule::FindMolChiralCenters(bool includeUnassigned)
{
  RDKit::MolOps::assignStereochemistry(*rdmol,false,true,true);
  int p = 0;
  const RDKit::Atom *atom;
  vector<string> centers;
    for (int i =0;i<rdmol->getNumAtoms();i++) { 
      atom = rdmol->getAtomWithIdx(i);
      std::cout << i  << std::endl; 
      vector<string>  plist=atom->getPropList();
      for (int x=0;x<plist.size();x++)
      {
      std::cout << plist[x]  << std::endl; 
      }
      if (atom->hasProp(RDKit::common_properties::_CIPCode)) {
         string cip;
         atom->getProp(RDKit::common_properties::_CIPCode, cip);
         std::cout << cip << std::endl; 

         centers[p] = i;
         centers[p+1] = cip;
         p=p+2;
      }
      else if (includeUnassigned and atom->hasProp(RDKit::common_properties::_ChiralityPossible)){
         centers[p] = i;
         centers[p+1] = '?';
         std::cout << "includeUnassigned true" << std::endl; 

      p=p+2;
    }
  }
  return centers;
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



vector<string> Molecule::getproplist()
{
    return rdmol->getPropList();
}

string Molecule::toSmiles()
{
    string smile =  RDKit::MolToSmiles(*rdmol);  
    return smile;

}

string Molecule::toMolfile()
{
    stringstream ss;
    RDKit::SDWriter *writer = new RDKit::SDWriter(&ss);
    writer->write(*rdmol);
    writer->flush();
    delete writer;
    return ss.str();
}


string Molecule::getPath()

{
    char cCurrentPath[FILENAME_MAX];
    int size =sizeof(cCurrentPath);
    getcwd(cCurrentPath, size);
    string res = string(cCurrentPath);
    return res;

}
    

string Molecule::sdwriteConfs()
{
    stringstream ss;
    RDKit::SDWriter *writer = new RDKit::SDWriter(&ss); // remove false argument!
    
    for(int i=0; i<rdmol->getNumConformers(); ++i){
         writer->write(*rdmol,i);
     }

    writer->flush();
    delete writer;
    return ss.str();
}


unsigned int Molecule::compute2DCoords()
{
    return RDDepict::compute2DCoords(*rdmol);
}



vector<float> Molecule::getAtomsPos2D()
{
    double sigma = 0;
    int atomnumber = rdmol->getNumAtoms();
    vector<float>  res(2*atomnumber+1);

    RDDepict::compute2DCoords(*rdmol);
    WedgeMolBonds(*rdmol,&(rdmol->getConformer()));
    RDKit::MolDraw2DSVG drawer(300,300);
    drawer.drawMolecule(*rdmol);
    drawer.finishDrawing();


    for (int i =0;i<2*atomnumber;i=i+2) { 
      RDGeom::Point2D atomcoords = drawer.getDrawCoords(i/2); // replace the getAtomsCoords by Draw
      res[i]=atomcoords[0]/300;  // rescale to 0-1 range!
      res[i+1]=atomcoords[1]/300;
    }
    

    if (rdmol->getNumBonds() >0){
      RDKit::Bond * bond = rdmol->getBondWithIdx(0);
      int idx1 = bond->getBeginAtomIdx();
      int idx2 = bond->getEndAtomIdx();
      sigma = 0.3*sqrt(pow(res[2*idx1]-res[2*idx2],2)+pow(res[2*idx1+1]-res[2*idx2+1],2));
    }

    else {
      sigma = 0.3*sqrt(pow(res[0]-res[2],2)+pow(res[1]-res[3],2));
    }
     // store sigma in the last position => odd list!
      res[2*atomnumber+1]=sigma;
    return res;


}



double Molecule::get2DScale(vector<float> atcds, double width, double height) {
    double x_min =10;
    double y_min =10;
    double x_max = -10;
    double y_max = -10;
    for( int i = 0; i < atcds.size() -1 ; i=i+2 ) {
      if (atcds[i]<x_min ) x_min = atcds[i];
      if (atcds[i]>x_max ) x_max = atcds[i];
      if (atcds[i+1]<y_min ) y_min = atcds[i+1];
      if (atcds[i+1]>y_max ) y_max = atcds[i+1];
    }

    double x_range = x_max - x_min;
    double y_range = y_max - y_min;
    //return scale = min( double( width ) / x_range , double( height ) / y_range );
    double scale;
    if (width/x_range > height/y_range) scale =  height/y_range;
    else scale = width/x_range;
    return scale ;
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


void Molecule::computeGasteigerCharges()
{
    vector<double> charges(rdmol->getNumAtoms(),0);
    RDKit::computeGasteigerCharges(*rdmol, charges,12,false);
}

// the results are not equivalent to the NPscore.py function need to check another method!!!
vector<int> Molecule::getSpiroBridgeMacrocycles()
{
    RDKit::VECT_INT_VECT rings = rdmol->getRingInfo()->atomRings();

    int macroCycleCount = 0;
    int bridgeAtomCount = 0;
    int spiroAtomCount = 0;
    int i=0;
    for (RDKit::VECT_INT_VECT_CI ringIter = rings.begin();ringIter != rings.end();ringIter++) {
        if (ringIter->size()>8) {
            macroCycleCount++;
        }
        std::sort(rings[i].begin(),rings[i].end());
        int x=1;
        for (RDKit::VECT_INT_VECT_CI ringIter2 = ringIter+1; ringIter2 != rings.end(); ringIter2++) {
            vector<int> v;
            RDKit::Intersect(rings[i],rings[i+x],v);
            if (v.size() == 1) {
                spiroAtomCount++;
            } else if (v.size() > 2) {
                bridgeAtomCount++;
            }
            x++;
        }
        i++;
    }

    vector<int> result(3);
    result[0]=spiroAtomCount;
    result[1]=bridgeAtomCount;
    result[2]=macroCycleCount;

    return result;
}



/*
string Molecule::sdwritefile(string filename)
{


string smiString = "CN([CH](Cc3ccc(OS(c2cccc1c2ccnc1)(=O)=O)cc3)C(=O)N5CCN(c4ccccc4)CC5)S(c7cccc6c7ccnc6)(=O)=O \
                           C[n](n1)ccc1NC(=O)C(/C#N)=C/c2[s]cc(c2)Br O=C2C1CC3CC(C1)CC2C3 \
                           c8cc9cc3c5c(c(=O)cc4n1nc6c2c(ccc(c12)c(cc3)c45)cc7c6cc(=O)cc7)c9cc8 \
                           c1cccc(c1)CCCCNC(=O)C(/C#N)=C/c2ccc(O)c(c2)O \
                           COc(cc1)ccc1C(=C2)C=C(NC2=C)c3ccc(OC)cc3O \
                           CC(OC1C(CCCC3)C3C(CCCC2)C2C1OC(C)=O)=O \
                           Cl.C[N+](C)(C)CCO O=C2C1CC3CC(C1)CC2C3 \
                           CN3CCC25[CH]4CCC(=O)C5(C)Oc1c(O)ccc(c12)C[CH]34 \
                           c1ccccc1\\C=C/C=C\\Cl.c1ccccc1C=C(Cl)C1SCCC1.Cl\\C=C1/SCCC1 Cl\\C=C/Br \
                           c1ccccc1\\C=C/C=C\\C=C/Cl c1ccccc1C=C(Cl)C1SCCC1 Cl\\C=C1/SCCC1 Cl\\C=C/Br \
                           CN2C3CC(OC(=O)C(CO)c1ccccc1)CC2CC3 \
                           N2C3CC(OC(=O)C(CO)c1ccccc1)CC2CC3 \
                           C2C3CC(OC(=O)C(CO)c1ccccc1)CC2CC3 \
                           ClC=C1SCCC1 C/C=C/C(C1CCCCCC1)=O c1ccccc1\\C=C/C=C\\C=C/Cl \
                           C[n](n1)ccc1NC(=O)C(/C#N)=C/c2[s]cc(c2)Br \
                           C1CC(C(Cl)(Br)F)CC1 \
                           C1(C2)CCC2C=C1 \
                           OC2C1(C)CCCCC1CCC2 \
                           CC12CCCC3OC1CCCC23 \
                           CC(OC1C(CCCC3)C3C(CCCC2)C2C1OC(C)=O)=O \
                           ON=C(CC1=CC=CC=C1)[CH](C#N)C2=CC=CC=C2 \
                           COc(cc1)ccc1C(=C2/C#N)\\C=C(NC2=C(C#N)C#N)\\c3ccc(OC)cc3O \
                           COc(cc1)ccc1C(=C2)C=C(NC2=C)c3ccc(OC)cc3O";



    string ofile = filename;
    RDKit::SDWriter *sdfWriter = new RDKit::SDWriter(ofile);
    boost::char_separator<char> spaceSep(" ");
    tokenizer tokens(smiString,spaceSep);
    int j=0;
    for(tokenizer::iterator token=tokens.begin();token!=tokens.end(); ++token){
        ++j;
        std::string smi=*token;
        RWMol *m = RDKit::SmilesToMol(smi, 0, 1); 
        sdfWriter->write(*m);
        delete m;
    }

    sdfWriter->flush();
    delete sdfWriter;

    return to_string(j);
}



string Molecule::sdreadfile(string filename)

 {
    string ofile = filename;

  RDKit::SDMolSupplier reader(ofile);
  int i = 0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    std::string mname;
    delete mol;
    i++;
  }

    return to_string(i);
}



void Molecule::save(string path, string data)
{
 std::string asm_code;
    asm_code += "localStorage.setItem( '";
    asm_code += path;
    asm_code += "', '";
    asm_code += data;
    asm_code += "' );";
    emscripten_run_script( asm_code.data() );

}


string Molecule::load(string path) 
{
    std::string asm_code;
    asm_code += "localStorage.getItem( '";
    asm_code += path;
    asm_code += "' );";
    std::string buffer = emscripten_run_script_string( asm_code.data() );
    return buffer;

}



int Molecule::readfile(string filename)

 {
  string line;

  ifstream myfile;
  myfile.open(filename);
  int j=0;
  if (myfile.is_open())
   {
     while ( getline (myfile,line)  && j<100)
     { 
       ++j;
       cout << line << '\n' ;
    }
    myfile.close();
  }

  else 
    std::cout << "Unable to open file" << std::endl; 

  return 1;


}




int Molecule::nodereadwritewithdata(string path, string data) {
    

    string asm_code;
    asm_code += "var fs = require('fs');";
    asm_code += "fs.writeFileSync('";
    asm_code += path;
    asm_code += "','";
    asm_code += data;
    asm_code += "');";
    emscripten_run_script( asm_code.data() );
    return 1;
}
*/

/*
int Molecule::nodereadwrite() {

  FILE *file;
  int res;
  char buffer[512];

  // write something locally with node
  EM_ASM(
    var fs = require('fs');
    fs.writeFileSync('foobar.txt', 'yeehaw');
  );

  // mount the current folder as a NODEFS instance
  // inside of emscripten
  EM_ASM(
    FS.mkdir('/working');
    FS.mount(NODEFS, { root: '.' }, '/working');
  );

  // read and validate the contents of the file
  file = fopen("/working/foobar.txt", "r");
  res = fread(buffer, sizeof(char), 6, file);
  fclose(file);

  // write out something new
  file = fopen("/working/foobar.txt", "w");
  res = fwrite("cheez", sizeof(char), 5, file);
  fclose(file);

  // validate the changes were persisted to the underlying fs
  EM_ASM(
    var fs = require('fs');
    var contents = fs.readFileSync('foobar.txt', { encoding: 'utf8' });

  );

  puts("success");

  return 0;
}









int Molecule::writefile(string filename, string data)

 {
    std::string asm_code;
    asm_code += "FS.writeFile('";
    asm_code += filename;
    asm_code += "', '";
    asm_code += data;
    asm_code += "' );";
    emscripten_run_script_string( asm_code.data() );


    std::ifstream file(filename);

    while(!file.eof() && !file.fail())
    {
        std::string line;
        getline(file, line);
        std::string name;
    
        std::cout << "read " << line << std::endl;

    }



  ofstream myfile;
  myfile.open (filename);
  myfile << data;
  myfile.close();
  

   return 1;
}



*/

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


int Molecule::findSSSR (std::vector< std::vector< int >> res)
{
    return RDKit::MolOps::findSSSR(*rdmol,res);
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

//cannot access the data of this vector in js! need to split it! 
vector<double> Molecule::logp_mr()
{    
    vector<double> res(2); 
    double logp, mr;
    RDKit::Descriptors::calcCrippenDescriptors (*rdmol,logp,mr);
    res[0] = logp;
    res[1] = mr;
    return res;
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

unsigned int Molecule::NumSaturatedRings()
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

unsigned int Molecule::NumAromaticCarbocycles()
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



vector<double> Molecule::getCrippenAtomContribs() {
    vector<double> logp(rdmol->getNumAtoms());
    vector<double> mr(rdmol->getNumAtoms());
    RDKit::Descriptors::getCrippenAtomContribs(*rdmol,logp,mr,true);
    vector<double> results;
    results.reserve(logp.size() + mr.size());
    results.insert(results.end(), logp.begin(), logp.end());
    results.insert(results.end(), mr.begin(), mr.end());
    return results;
}


vector<double>  Molecule::getTPSAAtomContribs() {
    vector<double> contribs(rdmol->getNumAtoms());
    RDKit::Descriptors::getTPSAAtomContribs(*rdmol,contribs,false);
    return contribs;
}



vector<double> Molecule::getASAContribs() {
    double hContrib=0.0;
    bool includeHs=true;
    vector<double> contribs(rdmol->getNumAtoms());
    RDKit::Descriptors::getLabuteAtomContribs(*rdmol,contribs,hContrib,includeHs,false);
    return contribs;
}

int Molecule::GetSubstructMatches(string smilesref)
{
    RDKit::MatchVectType matchV;
    vector< RDKit::MatchVectType > matches;

    rdErrorLog->df_enabled = false;
    RWMol* rdquery = RDKit::SmartsToMol(smilesref);

    int matched = RDKit::SubstructMatch(*rdmol,*rdquery,matches,true);
  /*  vector< int > res;
    for(int idx=0;idx<matched;idx++){
        res +=".";
    }*/
    return matched;
}




/*
double Molecule::GetConformersRMS(int confId1, int confId2){
    RDKit::Conformer  conf1 = rdmol->getConformer(confId1);
    RDKit::Conformer  conf2 = rdmol->getConformer(confId2);
    //RDKit::MolAlign::alignMolConformers(*rdmol,{confId1,confId2});
    double ssr = 0 ;
    double d;
    for (int i =0;i<rdmol->getNumAtoms();i++) { 
      d = RDGeom::Point3D::Distance(conf1.getAtomPos(i),conf2.getAtomPos(i)) ;
      ssr += d*d ;
    }
    ssr = ssr/rdmol->getNumAtoms() ;
    return ssr;
}
*/


double Molecule::AlignMol(string smilesref){
    rdErrorLog->df_enabled = false;
    RWMol* rdquery = RDKit::SmilesToMol(smilesref);
    printf("smiles ok");

    RDKit::MolOps::addHs(*rdquery);
    printf("addH ok");

    RDKit::DGeomHelpers::EmbedMolecule(*rdquery);
    printf("add 3D ok");

    pair<int, double> p = RDKit::MMFF::MMFFOptimizeMolecule(*rdquery);
    printf("optimze 3D ok");

    double res = RDKit::MolAlign::alignMol(*rdmol,*rdquery);
    printf("Align ok");
   // string sres = (string)res;
   // printf("%s\n", sres.c_str());
    // can return the rmds values need to look at that closely... in python not there
    return res;
}



void Molecule::AlignMolConformers(){
    // can return the rmds values need to look at that closely... in python not there
    RDKit::MolAlign::alignMolConformers(*rdmol);

}


/*
vector<double> Molecule::AlignMolConformersRMSlist(){
    // can return the rmds values need to look at that closely... in python not there
    std::vector<double> *RMSlist;
    RDKit::MolAlign::alignMolConformers(*rdmol,0,0,0,false,50,RMSlist);
    return *RMSlist;

}
*/

bool Molecule::HasSubstructMatchStr(string smilesref)
{
    rdErrorLog->df_enabled = false;
    RWMol* rdquery = RDKit::SmartsToMol(smilesref);
    
    RDKit::MatchVectType res;
    return RDKit::SubstructMatch(*rdmol,*rdquery,res);
}



/// get & set & has properties


string Molecule::getProp(string key) {
    string res;
    rdmol->getProp(key,res);
    return res;
}


string Molecule::getAtomProp(string key,int atomid) {
    string res;
    rdmol->getAtomWithIdx(atomid)->getProp(key,res);
    return res;
}


int Molecule::getNumConformers() {
    return rdmol->getNumConformers();
}


RDKit::Conformer Molecule::getConformer(int id) {
    return rdmol->getConformer(id);
}



int Molecule::setProp(string key, string value) {
    rdmol->setProp(key,value);
    return 0;
}

bool Molecule::hasProp(string key) {
    return rdmol->hasProp(key);
}
