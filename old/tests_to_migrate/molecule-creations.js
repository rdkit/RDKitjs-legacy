'use strict';

const RDKit = require('..');

beforeEach(() => RDKit.load());

// static constructors
//    newmolecule();
test('Murko', function() {
  var smi = 'C1CCCC1OC(CO)';
  var murcko = RDKit.Molecule.MurckosmilesToMol(smi);
  expect(murcko.molToMolfile()).toMatchSnapshot();
  murcko.delete();
});

test.skip('Mol2BlockToMol', function() {
  var molBlock = '';
  var mol = RDKit.Molecule.Mol2BlockToMol(molBlock);
  expect(mol.molToMolfile()).toBe();
  mol.delete();
});

test('MolBlockToMol', function() {
  var molBlock =
    '\n     RDKit          \n\n  9  9  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\n  5  1  1  0\nM  END\n$$$$\n';
  var mol = RDKit.Molecule.MolBlockToMol(molBlock);
  expect(mol.molToSmiles()).toBe('OCCOC1CCCC1');
  mol.delete();
});

test('smartsToMol', function() {
  var smarts = '[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]';
  var mol = RDKit.Molecule.smartsToMol(smarts);
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});

test('smilesToMol', function() {
  var smi = 'C1CCCC1OC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});

test('molFromPickle', function() {
  var pickle =
    'ï¾­Þ\u0000\u0000\u0000\u0000\u0007\u0000\u0000\u0000\u0002\u0000\u0000\u0000\u0000\u0000\u0000\u0000\t\u0000\u0000\u0000\t\u0000\u0000\u0000\u0001\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0003\u0001\b\u0000 \u0000\u0000\u0000\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\b\u0000`\u0000\u0000\u0000\u0001\u0001\u000b\u0000\u0001\u0000\u0001\u0002\u0000\u0002\u0003\u0000\u0003\u0004\u0000\u0004\u0005\u0000\u0005\u0006\u0000\u0006\u0007\u0000\u0007\b\u0000\u0004\u0000\u0000\u0014\u0001\u0005\u0000\u0001\u0002\u0003\u0004\u0017\u0000\u0000\u0000\u0000\u0016';
  var mol = RDKit.Molecule.molFromPickle(pickle);
  expect(mol.molToMolfile()).toMatchSnapshot();
  mol.delete();
});

// Pickle molecule representation
test('pickleMol', function() {
  var smi = 'C1CCCC1OC(CO)';
  var mol = RDKit.Molecule.smilesToMol(smi);
  expect(mol.pickleMol()).toMatchSnapshot();
  mol.delete();
});
