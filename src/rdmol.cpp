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


typedef boost::tokenizer<boost::char_separator<char> > tokenizer;




using namespace std;
using RDKit::ROMol;
using RDKit::RWMol;

Molecule::Molecule(RWMol* mol): rdmol(mol) {}


Molecule::~Molecule() {
  if(rdmol != 0)
    delete rdmol;
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

/*
double  Molecule::RusselSimilarityfromSmile ( string smilesref)       
{


   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return RusselSimilarityfromSmile (*v1,*v2);
}


double  Molecule::CosineSimilarityfromSmile ( string smilesref)       
{


   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return CosineSimilarity (*v1,*v2);
}


double  Molecule::KulczynskiSimilarityfromSmile ( string smilesref)       
{


   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return KulczynskiSimilarity (*v1,*v2);
}


double  Molecule::McConnaugheySimilarityfromSmile ( string smilesref)       
{


   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return McConnaugheySimilarity (*v1,*v2);
}


double  Molecule::SokalSimilarityfromSmile ( string smilesref)       
{


   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return SokalSimilarity (*v1,*v2);
}



double  Molecule::AsymmetricSimilarityfromSmile ( string smilesref)       
{


   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return AsymmetricSimilarity (*v1,*v2);
}



double  Molecule::BraunBlanquetSimilarityfromSmile ( string smilesref)       
{


   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return BraunBlanquetSimilarity (*v1,*v2);
}
*/


/*
double  Molecule::RogotGoldbergSimilarityfromSmile ( string smilesref)       
{


   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return RogotGoldbergSimilarity (*v1,*v2);
}
*/
/*
double  Molecule::OnBitSimilarityfromSmile ( string smilesref)       
{


   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return OnBitSimilarity (*v1,*v2);
}


int  Molecule::NumBitsInCommonfromSmile ( string smilesref)       
{


   ROMol *mol=RDKit::SmilesToMol(smilesref);
   RDKit::SparseIntVect< boost::uint32_t > * v1 =  RDKit::MorganFingerprints::getFingerprint(*rdmol,2);
   RDKit::SparseIntVect< boost::uint32_t > * v2 =  RDKit::MorganFingerprints::getFingerprint(*mol,2);

   return NumBitsInCommon (*v1,*v2);
}

*/





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


string Molecule::getRDKFP()
{
    ExplicitBitVect* finger =  RDKit::RDKFingerprintMol(*rdmol);
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


/*
SparseIntVect<boost::int32_t> *fp1;fp1 = AtomPairs::getAtomPairFingerprint(*m1);
fp1 = AtomPairs::getTopologicalTorsionFingerprint(*m1);
SparseIntVect<boost::int64_t> *fp1; fp1 = AtomPairs::getHashedTopologicalTorsionFingerprint(*m3,4096);
SparseIntVect<boost::uint32_t> *iv; iv = MorganFingerprints::getHashedFingerprint(*mol,2,2048,0,0,false,true,false,&bitInfo2);
  */




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


string Molecule::getPath()

{

    char cCurrentPath[FILENAME_MAX];
    int size =sizeof(cCurrentPath);

    getcwd(cCurrentPath, size);
    string res = string(cCurrentPath);

    return res;

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





string Molecule::sdwriteConfs()
{
    stringstream ss;
    RDKit::SDWriter *writer = new RDKit::SDWriter(&ss,false);
    
    for(int i=0; i<rdmol->getNumConformers(); ++i){
         writer->write(*rdmol,i);
     }

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

vector< double > Molecule::logp_mr()
{
    vector< double > v = new vector< double >[2]; 
    double logp;
    double mr;
    RDKit::Descriptors::calcCrippenDescriptors (*rdmol,logp,mr)

    v[0]=logp;
    v[1]=mr
    return v;
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
    RWMol* rdquery = RDKit::SmartsToMol(smilesref);

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
