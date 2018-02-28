#include <string>
#include <sstream>

#include <emscripten/bind.h>

#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

using emscripten::allow_raw_pointers;
using emscripten::class_;
using emscripten::function;

using RDKit::ROMol;
using RDKit::RWMol;
using RDKit::SDWriter;

RWMol *smilesToMol(std::string smiles)
{
  return RDKit::SmilesToMol(smiles);
}

std::string molToMolfile2D(RWMol *mol)
{
  std::stringstream ss;
  SDWriter writer(&ss);
  writer.write(*mol);
  writer.flush();
  return ss.str();
}

EMSCRIPTEN_BINDINGS(module)
{
  class_<RWMol>("RWMol");
  function("smilesToMol", &smilesToMol, allow_raw_pointers());
  function("molToMolfile2D", &molToMolfile2D, allow_raw_pointers());
}
