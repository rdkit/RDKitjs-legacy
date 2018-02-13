'use strict';

const RDKit = require('..');

beforeEach(() => RDKit.load());

test.skip('getConformer', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.addHs();
  mol.EmbedMultipleConfsarg(3, 100, 2015);
  expect(mol.getConformer(1)).toBe();
  mol.delete();
});

test('getNumAtoms', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  expect(mol.getNumAtoms()).toBe(9);
  mol.delete();
});

test('getAtomNeighbors', function() {
  var smi = 'C(C)C(CCC)CCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  var p = mol.getAtomNeighbors(3);
  var e = [];
  for (var j = 0; j < p.size(); j++) {
    e.push(p.get(j));
  }
  expect(e).toEqual([2, 4]);
  mol.delete();
});

test('getBondNeighbors', function() {
  var smi = 'C(C)C(CCC)CCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  var p = mol.getBondNeighbors(3);
  var e = [];
  for (var j = 0; j < p.size(); j++) {
    e.push(p.get(j));
  }

  expect(e).toEqual([1, 1]);
  mol.delete();
});

test('getNumConformers', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.addHs();
  mol.EmbedMultipleConfsarg(3, 100, 2015);
  expect(mol.getNumConformers()).toBe(3);
  mol.delete();
});
