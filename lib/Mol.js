'use strict';

const bindings = require('../dist/rdkit.js');

const mPtr = Symbol('Mol pointer');

function getEmbedMethodVal(method) {
  switch (method) {
    case 'manual':
      return 0;
    case 'ETDG':
      return 1;
    case 'ETKDG':
      return 2;
    case 'KDG':
      return 3;
    default:
      throw new TypeError(`unknown embed method: ${method}`);
  }
}

const embedDefaults = {
  method: 'manual',
  maxIterations: 0,
  seed: -1,
  clearConfs: true,
  useRandomCoords: false,
  boxSizeMult: 2
};

function getReplacementsMap(replacements) {
  if (replacements === null) return null;
  const replacementsMap = new bindings.map_string_string();
  for (const key in replacements) {
    replacementsMap.set(key, `${replacements[key]}`);
  }
  return replacementsMap;
}

class Mol {
  constructor(ptr) {
    this[mPtr] = ptr;
  }

  delete() {
    this[mPtr].delete();
  }

  /* Chem/rdchem */

  static fromPickle(pickle) {
    return new Mol(bindings.molFromPickle(pickle));
  }

  toPickle() {
    return bindings.pickleMol(this[mPtr]);
  }

  getNumAtoms(options = {}) {
    const { onlyExplicit = true } = options;
    return bindings.getNumAtoms(this[mPtr], onlyExplicit);
  }

  getNumBonds(options = {}) {
    const { onlyHeavy = true } = options;
    return bindings.getNumBonds(this[mPtr], onlyHeavy);
  }

  getNumConformers() {
    return bindings.getNumConformers(this[mPtr]);
  }

  getNumHeavyAtoms() {
    return bindings.getNumHeavyAtoms(this[mPtr]);
  }

  /* Chem/rdDistGeom */

  embedMolecule(options) {
    options = Object.assign({}, embedDefaults, options);
    const methodVal = getEmbedMethodVal(options.method);

    return bindings.EmbedMolecule(
      this[mPtr],
      methodVal,
      options.maxIterations,
      options.seed,
      options.clearConfs,
      options.useRandomCoords,
      options.boxSizeMult
    );
  }

  embedMultipleConfs(options) {
    options = Object.assign({ numConfs: 10 }, embedDefaults, options);
    const methodVal = getEmbedMethodVal(options.method);

    const result = bindings.EmbedMultipleConfs(
      this[mPtr],
      methodVal,
      options.numConfs,
      options.maxIterations,
      options.seed,
      options.clearConfs,
      options.useRandomCoords,
      options.boxSizeMult
    );

    const toReturn = [];
    const size = result.size();
    for (let i = 0; i < size; i++) {
      toReturn.push(result.get(i));
    }
    result.delete();
    return toReturn;
  }

  /* Chem/rdForceFieldHelpers */

  MMFFOptimizeMolecule(options = {}) {
    const {
      maxIters = 1000,
      mmffVariant = 'MMFF94',
      nonBondedThresh = 10,
      confId = -1,
      ignoreInterfragInteractions = true
    } = options;
    const result = bindings.MMFFOptimizeMolecule(
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

  /* rdMolAlign */
  alignMolConformers() {
    bindings.alignMolConformers(this[mPtr]);
  }

  /* Chem/rdmolfiles */

  static fromMolBlock(molBlock, options = {}) {
    const { sanitize = true, removeHs = true, strictParsing = true } = options;
    return new Mol(
      bindings.MolBlockToMol(molBlock, sanitize, removeHs, strictParsing)
    );
  }

  static fromSmarts(smarts, options = {}) {
    const { debugParse = 0, mergeHs = false, replacements = null } = options;

    const replacementsMap = getReplacementsMap(replacements);
    try {
      return new Mol(
        bindings.SmartsToMol(smarts, debugParse, mergeHs, replacementsMap)
      );
    } finally {
      if (replacementsMap !== null) {
        replacementsMap.delete();
      }
    }
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

    const replacementsMap = getReplacementsMap(replacements);
    try {
      return new Mol(
        bindings.SmilesToMol(
          smiles,
          debugParse,
          sanitize,
          replacementsMap,
          allowCXSMILES,
          parseName,
          removeHs
        )
      );
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

  deleteSubstructs(query, options = {}) {
    const { onlyFrags = false, useChirality = false } = options;
    return new Mol(
      bindings.deleteSubstructs(
        this[mPtr],
        query[mPtr],
        onlyFrags,
        useChirality
      )
    );
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
