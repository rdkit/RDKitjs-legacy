'use strict';

const bindings = require('../dist/rdkit.js');

const mPtr = Symbol('Molecule pointer');

class Molecule {
  constructor(ptr) {
    this[mPtr] = ptr;
  }

  static fromSmiles(smiles) {
    return new Molecule(bindings.smilesToMol(smiles));
  }

  toMolfile() {
    return bindings.molToMolfile2D(this[mPtr]);
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

  delete() {
    this[mPtr].delete();
  }
}

function smilesToMolfile(smiles) {
  const mol = bindings.smilesToMol(smiles);
  try {
    return bindings.molToMolfile2D(mol);
  } finally {
    mol.delete();
  }
}

function smilesTo3D(smiles) {
  var mol = Molecule.fromSmiles(smiles);
  try {
    mol.addHs();
    mol.EmbedMolecule();
    mol.MMFFoptimizeMolecule();
    return mol.toMolfile();
  } finally {
    mol.delete();
  }
}

const api = {
  Molecule,
  smilesToMolfile,
  smilesTo3D
};

module.exports = bindings.load().then(() => api);
