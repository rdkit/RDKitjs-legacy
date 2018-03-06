#include <string>

#include <GraphMol/MolPickler.h>
#include <GraphMol/RWMol.h>

using RDKit::RWMol;

RWMol *molFromPickle(std::string pickle)
{
  RWMol *result = new RWMol();
  RDKit::MolPickler::molFromPickle(pickle, result);
  return result;
}

std::string pickleMol(RWMol *mol)
{
  std::string res;
  RDKit::MolPickler::pickleMol(mol, res);
  return res;
}

unsigned int getNumAtoms(RWMol *mol, bool onlyExplicit)
{
  return mol->getNumAtoms(onlyExplicit);
}

unsigned int getNumBonds(RWMol *mol, bool onlyHeavy)
{
  return mol->getNumBonds(onlyHeavy);
}

unsigned int getNumConformers(RWMol *mol)
{
  return mol->getNumConformers();
}

unsigned int getNumHeavyAtoms(RWMol *mol)
{
  return mol->getNumHeavyAtoms();
}

#define BIND_Chem_rdchem()                                                 \
  {                                                                        \
    function("molFromPickle", &molFromPickle, allow_raw_pointers());       \
    function("pickleMol", &pickleMol, allow_raw_pointers());               \
    function("getNumAtoms", &getNumAtoms, allow_raw_pointers());           \
    function("getNumBonds", &getNumBonds, allow_raw_pointers());           \
    function("getNumConformers", &getNumConformers, allow_raw_pointers()); \
    function("getNumHeavyAtoms", &getNumHeavyAtoms, allow_raw_pointers()); \
  }
