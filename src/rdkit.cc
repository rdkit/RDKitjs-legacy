#include <string>

#include <emscripten/bind.h>

#include "Chem/rdchem.h"
#include "Chem/rdDistGeom.h"
#include "Chem/rdForceFieldHelpers.h"
#include "Chem/rdMolAlign.h"
#include "Chem/rdmolfiles.h"
#include "Chem/rdmolops.h"

using namespace emscripten;

EMSCRIPTEN_BINDINGS(module)
{
  class_<RDKit::RWMol>("RWMol");
  register_vector<double>("vector_double");
  register_vector<int>("vector_int");
  register_map<std::string, std::string>("map_string_string");
  BIND_Chem_rdchem();
  BIND_Chem_rdDistGeom();
  BIND_Chem_rdForceFieldHelpers();
  BIND_Chem_rdMolAlign();
  BIND_Chem_rdmolfiles();
  BIND_Chem_rdmolops();
}
