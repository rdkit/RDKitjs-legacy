'use strict';

const bindings = require('../dist/rdkit.js');

const mPtr = Symbol('Mol pointer');

class Mol {
  constructor(ptr) {
    this[mPtr] = ptr;
  }

  delete() {
    this[mPtr].delete();
  }

  /* Chem/rdDistGeom */

  embedMolecule(options = {}) {
    const {
      method = 'manual',
      maxIterations = 0,
      seed = -1,
      clearConfs = true,
      useRandomCoords = false,
      boxSizeMult = 2
    } = options;

    let methodVal;
    switch (method) {
      case 'manual':
        methodVal = 0;
        break;
      case 'ETDG':
        methodVal = 1;
        break;
      case 'ETKDG':
        methodVal = 2;
        break;
      case 'KDG':
        methodVal = 3;
        break;
      default:
        throw new TypeError(`unknown embed method: ${method}`);
    }

    return bindings.EmbedMolecule(
      this[mPtr],
      methodVal,
      maxIterations,
      seed,
      clearConfs,
      useRandomCoords,
      boxSizeMult
    );
  }

  /* Chem/rdForceFieldHelpers */

  MMFFoptimizeMolecule(options = {}) {
    const {
      maxIters = 1000,
      mmffVariant = 'MMFF94',
      nonBondedThresh = 10,
      confId = -1,
      ignoreInterfragInteractions = true
    } = options;
    const result = bindings.MMFFoptimizeMolecule(
      this[mPtr],
      maxIters,
      mmffVariant,
      nonBondedThresh,
      confId,
      ignoreInterfragInteractions
    );
    const toReturn = [result.get(0), result.get(1)];
    result.delete();
    return toReturn;
  }

  /* Chem/rdmolfiles */

  static fromMolBlock(molBlock, options = {}) {
    const { sanitize = true, removeHs = true, strictParsing = true } = options;
    return new Mol(
      bindings.MolBlockToMol(molBlock, sanitize, removeHs, strictParsing)
    );
  }

  static fromSmiles(smiles, options = {}) {
    const {
      debugParse = 0,
      sanitize = true,
      replacements = null,
      allowCXSMILES = false,
      parseName = false,
      removeHs = true
    } = options;

    let replacementsMap = null;
    if (replacements !== null) {
      replacementsMap = new bindings.map_string_string();
      for (const key in replacements) {
        replacementsMap.set(key, `${replacements[key]}`);
      }
    }

    try {
      const mol = bindings.SmilesToMol(
        smiles,
        debugParse,
        sanitize,
        replacementsMap,
        allowCXSMILES,
        parseName,
        removeHs
      );

      return new Mol(mol);
    } finally {
      if (replacementsMap !== null) {
        replacementsMap.delete();
      }
    }
  }

  toFASTA() {
    return bindings.MolToHELM(this[mPtr]);
  }

  toHELM() {
    return bindings.MolToHELM(this[mPtr]);
  }

  toMolBlock(options = {}) {
    const {
      includeStereo = false,
      confId = -1,
      kekulize = true,
      forceV3000 = false
    } = options;
    return bindings.MolToMolBlock(
      this[mPtr],
      includeStereo,
      confId,
      kekulize,
      forceV3000
    );
  }

  toSmarts(options = {}) {
    const { isomericSmarts = false } = options;
    return bindings.MolToSmarts(this[mPtr], isomericSmarts);
  }

  toSmiles(options = {}) {
    const {
      isomericSmiles = false,
      kekuleSmiles = false,
      rootedAtAtom = -1,
      canonical = true,
      allBondsExplicit = false,
      allHsExplicit = false
    } = options;
    return bindings.MolToSmiles(
      this[mPtr],
      isomericSmiles,
      kekuleSmiles,
      rootedAtAtom,
      canonical,
      allBondsExplicit,
      allHsExplicit
    );
  }

  /* Chem/rdmolops */

  addHs(options = {}) {
    const {
      explicitOnly = false,
      addCoords = false
      // onlyOnAtoms = null
    } = options;
    return bindings.addHs(this[mPtr], explicitOnly, addCoords);
  }

  removeHs(options = {}) {
    const {
      implicitOnly = false,
      updateExplicitCount = false,
      sanitize = true
    } = options;
    return bindings.removeHs(
      this[mPtr],
      implicitOnly,
      updateExplicitCount,
      sanitize
    );
  }
}

module.exports = Mol;
