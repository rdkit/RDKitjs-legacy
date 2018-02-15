'use strict';

var RDKit = require('..');

beforeEach(() => RDKit.load());

test('addHs', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  mol.addHs();
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});

test('removeHs', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  mol.removeHs();
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});

test('sanitizeMol', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  mol.sanitizeMol();
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});

test('cleanUp', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  mol.cleanUp();
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});

test('Kekulize', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  mol.Kekulize();
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});
