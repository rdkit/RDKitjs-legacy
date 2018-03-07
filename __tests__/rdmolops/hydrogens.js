'use strict';

const RDKit = require('../..');

beforeEach(() => RDKit.load());

test('add and remove hydrogens', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  mol.addHs();
  expect(mol.toSmiles()).toBe('[H]OC([H])([H])OC([H])([H])[H]');
  mol.removeHs();
  expect(mol.toSmiles()).toBe('COCO');
});
