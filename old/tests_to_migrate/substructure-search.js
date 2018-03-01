'use strict';

var RDKit = require('..');

beforeEach(() => RDKit.load());

test('getSubstructMatches', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  var smilesref = 'CO';
  var d = mol.getSubstructMatches(smilesref);
  var dlen = d.size();
  var e = [];
  for (var j = 0; j < dlen; j++) {
    e.push(d.get(j));
  }

  expect(e).toEqual([4, 5, 6, 5, 7, 8]);
  mol.delete;
});

test('HasSubstructMatchStr', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  var smilesref = 'CO';
  expect(mol.HasSubstructMatchStr(smilesref)).toBe(true);
  mol.delete;
});
