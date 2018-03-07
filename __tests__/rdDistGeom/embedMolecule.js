'use strict';

const RDKit = require('../..');

beforeEach(() => RDKit.load());

test('manual method', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  mol.addHs();
  mol.embedMolecule();
  expect(mol.toMolBlock()).toMatchSnapshot();
});

test('ETDG method', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  mol.addHs();
  mol.embedMolecule({ method: 'ETDG' });
  expect(mol.toMolBlock()).toMatchSnapshot();
});

test('ETKDG method', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  mol.addHs();
  mol.embedMolecule({ method: 'ETKDG' });
  expect(mol.toMolBlock()).toMatchSnapshot();
});

test('KDG method', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  mol.addHs();
  mol.embedMolecule({ method: 'KDG' });
  expect(mol.toMolBlock()).toMatchSnapshot();
});

test('wrong method', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  mol.addHs();
  expect(() => mol.embedMolecule({ method: 'wrong' })).toThrow(
    /unknown embed method: wrong/
  );
});
