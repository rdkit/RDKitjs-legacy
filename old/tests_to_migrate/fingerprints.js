'use strict';

const RDKit = require('..');

beforeEach(() => RDKit.load());

test('getRDKFingerprintMol', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  expect(mol.getRDKFingerprintMol()).toMatchSnapshot();
  mol.delete();
});

test('getMorganFingerprints', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  expect(mol.getMorganFingerprints(2, 2048)).toMatchSnapshot();
  mol.delete();
});

test('getMorganFingerprints_getOnBbits', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  var mpf2 = mol.getMorganFingerprints_getOnBbits(2, 2048);
  var mpf = [];
  var mpflen = mpf2.size();
  for (var j = 0; j < mpflen; j++) {
    mpf.push(mpf2.get(j));
  }
  expect(mpf).toMatchSnapshot();
  mol.delete();
});

test('getMorganFingerprintslist', function() {
  var mol = RDKit.Molecule.smilesToMol('CCCC(CO)CCO');
  var u = mol.getMorganFingerprintslist(2);
  var fplen = u.size();
  var v = [];
  for (var i = 0; i < fplen; i++) {
    v[i] = u.get(i);
  }
  expect(v).toMatchSnapshot();
  mol.delete();
});

test.skip('getMorganFingerprints_getNonzeroElements', function() {
  var mol = RDKit.Molecule.smilesToMol('CCCC(CO)CCO');
  var f = mol.getMorganFingerprints_getNonzeroElements(2);
  var p = [];
  for (i = 0; i < f.size(); i++) {
    p.push(f.get(i));
  }
  expect(p).toEqual();
  mol.delete();
});

test('getHashedAtomPairFingerprintAsBitVect', function() {
  var mol = RDKit.Molecule.smilesToMol(
    'COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21'
  );
  var e = mol.getHashedAtomPairFingerprintAsBitVect(2048, 0, 1);
  expect(e).toMatchSnapshot();
  mol.delete();
});

test('getAtomCode', function() {
  var mol = RDKit.Molecule.smilesToMol(
    'COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21'
  );
  var d = mol.getAtomCode(0);
  expect(d).toBe(33);
  mol.delete();
});

test('getAtomPairFingerprint', function() {
  var mol = RDKit.Molecule.smilesToMol('C1CCCC1CO');
  var f = mol.getAtomPairFingerprint();
  var p = [];
  for (var i = 0; i < f.size(); i++) {
    p.push(f.get(i));
  }
  expect(p).toMatchSnapshot();
  mol.delete();
});

test('getLayeredFingerprintMol', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  expect(mol.getLayeredFingerprintMol(2, 2, 2048)).toMatchSnapshot();
  mol.delete();
});

test('getMACCSFingerprints', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  expect(mol.getMACCSFingerprints()).toMatchSnapshot();
  mol.delete();
});

test('getPatternFingerprintMol', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  expect(mol.getPatternFingerprintMol()).toMatchSnapshot();
  mol.delete();
});
