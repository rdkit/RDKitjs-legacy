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

#define BIND_Chem_rdchem()                                           \
  {                                                                  \
    function("molFromPickle", &molFromPickle, allow_raw_pointers()); \
    function("pickleMol", &pickleMol, allow_raw_pointers());         \
  }
