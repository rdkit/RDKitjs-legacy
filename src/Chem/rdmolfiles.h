#include <string>

#include <GraphMol/RWMol.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/SequenceWriters.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using std::string;

using RDKit::RWMol;

RWMol *MolBlockToMol(string molBlock,
                     bool sanitize,
                     bool removeHs,
                     bool strictParsing)
{
  return RDKit::MolBlockToMol(molBlock,
                              sanitize,
                              removeHs,
                              strictParsing);
}

RWMol *SmartsToMol(string smarts,
                   int debugParse,
                   bool mergeHs,
                   std::map<string, string> *replacements)
{
  return RDKit::SmartsToMol(smarts,
                            debugParse,
                            mergeHs,
                            replacements);
}

RWMol *SmilesToMol(string smiles,
                   int debugParse,
                   bool sanitize,
                   std::map<string, string> *replacements,
                   bool allowCXSMILES,
                   bool parseName,
                   bool removeHs)
{
  RDKit::SmilesParserParams params;
  params.debugParse = debugParse;
  params.sanitize = sanitize;
  params.replacements = replacements;
  params.allowCXSMILES = allowCXSMILES;
  params.parseName = parseName;
  params.removeHs = removeHs;
  return RDKit::SmilesToMol(smiles, params);
}

string MolToFASTA(RWMol *mol)
{
  return RDKit::MolToFASTA(*mol);
}

string MolToHELM(RWMol *mol)
{
  return RDKit::MolToHELM(*mol);
}

string MolToMolBlock(RWMol *mol,
                     bool includeStereo,
                     int confId,
                     bool kekulize,
                     bool forceV3000)
{
  return RDKit::MolToMolBlock(*mol,
                              includeStereo,
                              confId,
                              kekulize,
                              forceV3000);
}

string MolToSmarts(RWMol *mol,
                   bool isomericSmarts)
{
  return RDKit::MolToSmarts(*mol,
                            isomericSmarts);
}

string MolToSmiles(RWMol *mol,
                   bool isomericSmiles,
                   bool kekuleSmiles,
                   int rootedAtAtom,
                   bool canonical,
                   bool allBondsExplicit,
                   bool allHsExplicit)
{
  return RDKit::MolToSmiles(*mol,
                            isomericSmiles,
                            kekuleSmiles,
                            rootedAtAtom,
                            canonical,
                            allBondsExplicit,
                            allHsExplicit);
}

#define BIND_Chem_rdmolfiles()                                       \
  {                                                                  \
    function("MolBlockToMol", &MolBlockToMol, allow_raw_pointers()); \
    function("SmartsToMol", &SmartsToMol, allow_raw_pointers());     \
    function("SmilesToMol", &SmilesToMol, allow_raw_pointers());     \
    function("MolToFASTA", &MolToFASTA, allow_raw_pointers());       \
    function("MolToHELM", &MolToHELM, allow_raw_pointers());         \
    function("MolToMolBlock", &MolToMolBlock, allow_raw_pointers()); \
    function("MolToSmarts", &MolToSmarts, allow_raw_pointers());     \
    function("MolToSmiles", &MolToSmiles, allow_raw_pointers());     \
  }
