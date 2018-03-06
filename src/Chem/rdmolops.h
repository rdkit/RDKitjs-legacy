#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

using RDKit::ROMol;
using RDKit::RWMol;

void addHs(RWMol *mol,
           bool explicitOnly,
           bool addCoords)
{
  return RDKit::MolOps::addHs(*mol,
                              explicitOnly,
                              addCoords);
}

RWMol *deleteSubstructs(RWMol *mol,
                        RWMol *query,
                        bool onlyFrags,
                        bool useChirality)
{
  ROMol *result = RDKit::deleteSubstructs(*mol,
                                          *query,
                                          onlyFrags,
                                          useChirality);
  return new RWMol(*result);
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

#define BIND_Chem_rdmolops()                                               \
  {                                                                        \
    function("addHs", &addHs, allow_raw_pointers());                       \
    function("deleteSubstructs", &deleteSubstructs, allow_raw_pointers()); \
    function("removeHs", &removeHs, allow_raw_pointers());                 \
  }
