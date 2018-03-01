#include <string>
#include <sstream>

#include <emscripten/bind.h>

#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
// RDkit 3D
#include <GraphMol/DistGeomHelpers/Embedder.h>
// RDkit ForceField
// comments thegodone & Paolo => MMFF.h and Builder.h need to be patch to avoid class issues! 13_05_2015
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <ForceField/MMFF/Params.h>

using emscripten::allow_raw_pointers;
using emscripten::class_;
using emscripten::function;
using emscripten::register_vector;

using RDKit::RWMol;

RWMol *smilesToMol(std::string smiles)
{
  return RDKit::SmilesToMol(smiles);
}

std::string molToMolfile2D(RWMol *mol)
{
  std::stringstream ss;
  RDKit::SDWriter writer(&ss);
  writer.write(*mol);
  writer.flush();
  return ss.str();
}

void addHs(RWMol *mol)
{
  return RDKit::MolOps::addHs(mol);
}

int EmbedMolecule(RWMol *mol, unsigned int maxIterations, int seed)
{
  return RDKit::DGeomHelpers::EmbedMolecule(mol, maxIterations, seed);
}

std::vector<double> MMFFoptimizeMolecule(RWMol *mol)
{
  std::vector<double> res(2);
  std::pair<int, double> p = RDKit::MMFF::MMFFOptimizeMolecule(*mol);
  res[0] = static_cast<double>(p.first);
  res[1] = p.second;
  return res;
}

EMSCRIPTEN_BINDINGS(module)
{
  class_<RWMol>("RWMol");

  register_vector<double>("VectorDouble");

  function("smilesToMol", &smilesToMol, allow_raw_pointers());
  function("molToMolfile2D", &molToMolfile2D, allow_raw_pointers());
  function("addHs", &addHs, allow_raw_pointers());
  function("EmbedMolecule", &EmbedMolecule, allow_raw_pointers());
  function("MMFFoptimizeMolecule", &MMFFoptimizeMolecule, allow_raw_pointers());
}
