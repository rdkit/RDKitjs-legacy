'use strict';

var RDKit = require('..');

beforeEach(() => RDKit.load());

test('in a list of conformers', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.addHs();
  mol.EmbedMultipleConfsarg(3, 100, 2015);

  mol.AlignMolConformers();

  mol.delete();
  // to do return the aligned score !!!!
});

test('Mol vs Ref', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  expect(mol.AlignMol('CCCCC')).toBe(0.3789185721593953);
  // to do return the 2 molecules aligned!
  mol.delete();
});

test('RMS value of conformers', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.addHs();
  mol.EmbedMultipleConfsarg(10, 100, 2015);
  mol.AlignMolConformers();
  var p = mol.getConformersRMS(1, 2, 50);

  var plen = p.size();
  var e = [];
  for (var j = 0; j < plen; j++) {
    e.push(p.get(j));
  }
  // to do return the 2 molecules aligned!
  mol.delete();
});
