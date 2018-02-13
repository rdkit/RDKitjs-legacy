'use strict';

var RDKit = require('..');

beforeEach(() => RDKit.load());

test('Tanimoto', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  var smilesref = 'CCCCCOC';
  expect(mol.TanimotoSimilarityfromSmile(smilesref)).toBe(0.4666666666666667);
  mol.delete();
});

test('Dice', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  var smilesref = 'CCCCCOC';
  expect(mol.DiceSimilarityfromSmile(smilesref)).toBe(0.6363636363636364);
  mol.delete();
});

test('Tversky', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  var smilesref = 'CCCCCOC';
  expect(mol.TverskySimilarityfromSmile(smilesref, 2, 1)).toBe(
    0.34146341463414637
  );
  mol.delete();
});
