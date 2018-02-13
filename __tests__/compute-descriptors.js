'use strict';

const RDKit = require('..');

beforeEach(() => RDKit.load());

test('calc_all_desc', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  var c = RDKit.calc_all_desc(mol);
  c = JSON.stringify(c);
  expect(c).toMatchSnapshot();
  mol.delete();
});

test('getCrippenAtomContribs', function() {
  var smi = 'CO';
  var mol = RDKit.Molecule.fromSmiles(smi);
  var cp = mol.getCrippenAtomContribs();
  var LP = [];
  for (i = 0; i < cp.size() / 2; i++) {
    LP.push(cp.get(i));
  }
  var MR = [];
  for (i = cp.size() / 2; i < cp.size(); i++) {
    MR.push(cp.get(i));
  }
  expect(LP).toEqual([-0.2035, -0.2893]);
  expect(MR).toEqual([2.753, 0.8238]);

  mol.delete();
});

test('getTPSAAtomContribs', function() {
  var smi = 'CCCCCOC';
  var mol = RDKit.Molecule.fromSmiles(smi);
  var cp = mol.getTPSAAtomContribs();
  var TPSA = [];
  for (i = 0; i < cp.size(); i++) {
    TPSA.push(cp.get(i));
  }
  expect(TPSA).toEqual([0, 0, 0, 0, 0, 9.23, 0]);
  mol.delete();
});

test('getASAContribs', function() {
  var smi = 'CCC';
  var mol = RDKit.Molecule.fromSmiles(smi);
  var cp = mol.getASAContribs();
  var ASA = [];
  for (i = 0; i < cp.size(); i++) {
    ASA.push(cp.get(i));
  }
  expect(ASA).toMatchSnapshot();
  mol.delete();
});
