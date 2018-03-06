#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>

using RDKit::RWMol;

void addHs(RWMol *mol,
           bool explicitOnly,
           bool addCoords)
{
  return RDKit::MolOps::addHs(*mol,
                              explicitOnly,
                              addCoords);
}

void removeHs(RWMol *mol,
              bool implicitOnly,
              bool updateExplicitCount,
              bool sanitize)
{
  return RDKit::MolOps::removeHs(*mol,
                                 implicitOnly,
                                 updateExplicitCount,
                                 sanitize);
}

#define BIND_Chem_rdmolops()                               \
  {                                                        \
    function("addHs", &addHs, allow_raw_pointers());       \
    function("removeHs", &removeHs, allow_raw_pointers()); \
  }
