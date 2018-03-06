#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/RWMol.h>

using RDKit::RWMol;

int EmbedMolecule(RWMol *mol,
                  int method,
                  unsigned int maxIterations,
                  int seed,
                  bool clearConfs,
                  bool useRandomCoords,
                  double boxSizeMult)
{
  if (method == 0)
  {
    return RDKit::DGeomHelpers::EmbedMolecule(*mol,
                                              maxIterations,
                                              seed,
                                              clearConfs,
                                              useRandomCoords,
                                              boxSizeMult);
  }
  else if (method == 1)
  {
    return RDKit::DGeomHelpers::EmbedMolecule(*mol, RDKit::DGeomHelpers::ETDG);
  }
  else if (method == 2)
  {
    return RDKit::DGeomHelpers::EmbedMolecule(*mol, RDKit::DGeomHelpers::ETKDG);
  }
  else if (method == 3)
  {
    return RDKit::DGeomHelpers::EmbedMolecule(*mol, RDKit::DGeomHelpers::KDG);
  }

  return 0;
}

#define BIND_Chem_rdDistGeom()                                       \
  {                                                                  \
    function("EmbedMolecule", &EmbedMolecule, allow_raw_pointers()); \
  }
