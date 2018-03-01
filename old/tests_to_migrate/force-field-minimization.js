'use strict';

var RDKit = require('..');

beforeEach(() => RDKit.load());

test('MMFFoptimizeMolecule', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  mol.addHs();
  mol.EmbedMolecule();
  mol.MMFFoptimizeMolecule();
  mol.removeHs();
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});

test('MMFFoptimizeMolecule with paramaters', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  mol.addHs();
  mol.EmbedMolecule();
  mol.MMFFoptimizeMoleculearg(1000, 'MMFF94');
  mol.removeHs();
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});

test('MMFFOptimizeMoleculeConfs', function() {
  var smi = 'CCCCCOC(CO)';

  var mol = RDKit.Molecule.smilesToMol(smi);
  mol.addHs();
  mol.EmbedMultipleConfsarg(3, 1000, 2015);
  mol.MMFFOptimizeMoleculeConfs(8, 1000, 'MMFF94');
  mol.removeHs();
  expect(mol.sdwriteConfs()).toMatchSnapshot();
  mol.delete();
});

test('UFFOptimizeMolecule', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  mol.addHs();
  mol.EmbedMolecule();
  mol.UFFOptimizeMolecule();
  mol.removeHs();
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});

test('EmbedMolecule', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  mol.addHs();
  mol.EmbedMoleculearg(1000, 1);
  mol.removeHs();
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});
