#include <string>
#include <sstream>

#include <emscripten/bind.h>

#include <ForceField/MMFF/Params.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>

#include "Chem/rdmolfiles.h"

using std::string;
using std::vector;

using namespace emscripten;

using RDKit::RWMol;

void addHs(RWMol *mol, bool explicitOnly, bool addCoords)
{
  return RDKit::MolOps::addHs(*mol, explicitOnly, addCoords);
}

int EmbedMolecule(RWMol *mol, unsigned int maxIterations, int seed, bool clearConfs)
{
  return RDKit::DGeomHelpers::EmbedMolecule(*mol, maxIterations, seed, clearConfs);
}

vector<double> MMFFoptimizeMolecule(RWMol *mol, int maxIters, string mmffVariant, double nonBondedThresh)
{
  vector<double> res(2);
  std::pair<int, double> p = RDKit::MMFF::MMFFOptimizeMolecule(*mol, maxIters, mmffVariant, nonBondedThresh);
  res[0] = static_cast<double>(p.first);
  res[1] = p.second;
  return res;
}

EMSCRIPTEN_BINDINGS(module)
{
  class_<RWMol>("RWMol");
  register_vector<double>("vector_double");
  register_map<string, string>("map_string_string");
  function("addHs", &addHs, allow_raw_pointers());
  function("EmbedMolecule", &EmbedMolecule, allow_raw_pointers());
  function("MMFFoptimizeMolecule", &MMFFoptimizeMolecule, allow_raw_pointers());
  BIND_Chem_rdmolfiles();
}
