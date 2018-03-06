#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/RWMol.h>

using RDKit::RWMol;

void alignMolConformers(RWMol *mol)
{
  RDKit::MolAlign::alignMolConformers(*mol);
}

#define BIND_Chem_rdMolAlign()                                                 \
  {                                                                            \
    function("alignMolConformers", &alignMolConformers, allow_raw_pointers()); \
  }
