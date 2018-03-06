#include <vector>

#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/RWMol.h>

using std::string;
using std::vector;

using RDKit::RWMol;

vector<double> MMFFOptimizeMolecule(RWMol *mol,
                                    int maxIters,
                                    string mmffVariant,
                                    double nonBondedThresh,
                                    int confId,
                                    bool ignoreInterfragInteractions)
{
  vector<double> res(2);
  std::pair<int, double> p = RDKit::MMFF::MMFFOptimizeMolecule(*mol,
                                                               maxIters,
                                                               mmffVariant,
                                                               nonBondedThresh,
                                                               confId,
                                                               ignoreInterfragInteractions);
  res[0] = static_cast<double>(p.first);
  res[1] = p.second;
  return res;
}

#define BIND_Chem_rdForceFieldHelpers()                                            \
  {                                                                                \
    function("MMFFOptimizeMolecule", &MMFFOptimizeMolecule, allow_raw_pointers()); \
  }
