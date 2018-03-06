#include <vector>

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/RWMol.h>

using RDKit::RWMol;

RDKit::DGeomHelpers::EmbedParameters getEmbedParameters(int method)
{
  if (method == 1)
  {
    return RDKit::DGeomHelpers::ETDG;
  }
  else if (method == 2)
  {
    return RDKit::DGeomHelpers::ETKDG;
  }
  else
  {
    return RDKit::DGeomHelpers::KDG;
  }
}

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
  else
  {
    return RDKit::DGeomHelpers::EmbedMolecule(*mol,
                                              getEmbedParameters(method));
  }
}

std::vector<int> EmbedMultipleConfs(RWMol *mol,
                                    int method,
                                    unsigned int numConfs,
                                    unsigned int maxIterations,
                                    int seed,
                                    bool clearConfs,
                                    bool useRandomCoords,
                                    double boxSizeMult)
{
  if (method == 0)
  {
    return RDKit::DGeomHelpers::EmbedMultipleConfs(*mol,
                                                   numConfs,
                                                   maxIterations,
                                                   seed,
                                                   clearConfs,
                                                   useRandomCoords,
                                                   boxSizeMult);
  }
  else
  {
    return RDKit::DGeomHelpers::EmbedMultipleConfs(*mol,
                                                   numConfs,
                                                   getEmbedParameters(method));
  }
}

#define BIND_Chem_rdDistGeom()                                                 \
  {                                                                            \
    function("EmbedMolecule", &EmbedMolecule, allow_raw_pointers());           \
    function("EmbedMultipleConfs", &EmbedMultipleConfs, allow_raw_pointers()); \
  }
