#include <string>

#include <emscripten/bind.h>

#include "Chem/rdDistGeom.h"
#include "Chem/rdForceFieldHelpers.h"
#include "Chem/rdmolfiles.h"
#include "Chem/rdmolops.h"

using namespace emscripten;

EMSCRIPTEN_BINDINGS(module)
{
  class_<RDKit::RWMol>("RWMol");
  register_vector<double>("vector_double");
  register_map<std::string, std::string>("map_string_string");
  BIND_Chem_rdDistGeom();
  BIND_Chem_rdForceFieldHelpers();
  BIND_Chem_rdmolfiles();
  BIND_Chem_rdmolops();
}
