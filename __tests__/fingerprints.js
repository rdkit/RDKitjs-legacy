'use strict';

const RDKit = require('..');

beforeEach(() => RDKit.load());

test('getRDKFP', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  expect(mol.getRDKFP()).toMatchSnapshot();
  mol.delete();
});

test('getMorganFP', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  expect(mol.getMorganFP(2, 2048)).toMatchSnapshot();
  mol.delete();
});

test('getMorganFP_GetOnBits', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  var mpf2 = mol.getMorganFP_GetOnBits(2, 2048);
  var mpf = [];
  var mpflen = mpf2.size();
  for (var j = 0; j < mpflen; j++) {
    mpf.push(mpf2.get(j));
  }
  expect(mpf).toMatchSnapshot();
  mol.delete();
});

test('getMorganFPlist', function() {
  var mol = RDKit.Molecule.fromSmiles('CCCC(CO)CCO');
  var u = mol.getMorganFPlist(2);
  var fplen = u.size();
  var v = [];
  for (var i = 0; i < fplen; i++) {
    v[i] = u.get(i);
  }
  expect(v).toMatchSnapshot();
  mol.delete();
});

test.skip('getMorganFP_getNonzeroElements', function() {
  var mol = RDKit.Molecule.fromSmiles('CCCC(CO)CCO');
  var f = mol.getMorganFP_getNonzeroElements(2);
  var p = [];
  for (i = 0; i < f.size(); i++) {
    p.push(f.get(i));
  }
  expect(p).toEqual();
  mol.delete();
});

test('getHashedAtomPairFingerprintAsBitVect', function() {
  var mol = RDKit.Molecule.fromSmiles(
    'COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21'
  );
  var e = mol.getHashedAtomPairFingerprintAsBitVect(2048, 0, 1);
  expect(e).toMatchSnapshot();
  mol.delete();
});

test('getAtomCode', function() {
  var mol = RDKit.Molecule.fromSmiles(
    'COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21'
  );
  var d = mol.getAtomCode(0);
  expect(d).toBe(33);
  mol.delete();
});

test('getAtomPairFingerprint', function() {
  var mol = RDKit.Molecule.fromSmiles('C1CCCC1CO');
  var f = mol.getAtomPairFingerprint();
  var p = [];
  for (var i = 0; i < f.size(); i++) {
    p.push(f.get(i));
  }
  expect(p).toMatchSnapshot();
  mol.delete();
});

test('getLayeredFP', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  expect(mol.getLayeredFP(2, 2, 2048)).toMatchSnapshot();
  mol.delete();
});

test('getMACCSFP', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  expect(mol.getMACCSFP()).toMatchSnapshot();
  mol.delete();
});

test('getPatternFP', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  expect(mol.getPatternFP()).toMatchSnapshot();
  mol.delete();
});
