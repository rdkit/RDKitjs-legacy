#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>

using RDKit::RWMol;

void addHs(RWMol *mol, bool explicitOnly, bool addCoords)
{
  return RDKit::MolOps::addHs(*mol, explicitOnly, addCoords);
}

#define BIND_Chem_rdmolops()                         \
  {                                                  \
    function("addHs", &addHs, allow_raw_pointers()); \
  }
