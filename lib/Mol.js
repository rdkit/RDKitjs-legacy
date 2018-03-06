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

  addHs(options = {}) {
    const { explicitOnly = false, addCoords = false } = options;
    return bindings.addHs(this[mPtr], explicitOnly, addCoords);
  }

  EmbedMolecule(options = {}) {
    const { maxIterations = 0, seed = -1, clearConfs = true } = options;
    return bindings.EmbedMolecule(this[mPtr], maxIterations, seed, clearConfs);
  }

  MMFFoptimizeMolecule(options = {}) {
    const {
      maxIters = 1000,
      mmffVariant = 'MMFF94',
      nonBondedThresh = 10
    } = options;
    const result = bindings.MMFFoptimizeMolecule(
      this[mPtr],
      maxIters,
      mmffVariant,
      nonBondedThresh
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
}

module.exports = Mol;
