'use strict';

var RDKit = require('..');

beforeEach(() => RDKit.load());

test('addHs', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.addHs();
  expect(mol.toMolfile()).toMatchSnapshot();
  mol.delete();
});

test('removeHs', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.removeHs();
  expect(mol.toMolfile()).toMatchSnapshot();
  mol.delete();
});

test('sanitizeMol', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.sanitizeMol();
  expect(mol.toMolfile()).toMatchSnapshot();
  mol.delete();
});

test('cleanUp', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.cleanUp();
  expect(mol.toMolfile()).toMatchSnapshot();
  mol.delete();
});

test('Kekulize', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.Kekulize();
  expect(mol.toMolfile()).toMatchSnapshot();
  mol.delete();
});
