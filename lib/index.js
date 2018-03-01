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

  addHs() {
    return bindings.addHs(this[mPtr]);
  }

  EmbedMolecule() {
    return bindings.EmbedMolecule(this[mPtr], 50, 0);
  }

  MMFFoptimizeMolecule() {
    return bindings.MMFFoptimizeMolecule(this[mPtr]);
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
