'use strict';

let RDKit;
beforeEach(async () => (RDKit = await require('../..')));

test('check some numbers', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);

  expect(mol.getNumAtoms()).toBe(4);
  expect(mol.getNumAtoms({ onlyExplicit: false })).toBe(10);

  expect(mol.getNumBonds()).toBe(3);
  expect(mol.getNumBonds({ onlyHeavy: false })).toBe(9);

  expect(mol.getNumConformers()).toBe(0);

  expect(mol.getNumHeavyAtoms()).toBe(4);
});
