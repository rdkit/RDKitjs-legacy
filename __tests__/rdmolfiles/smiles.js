'use strict';

let RDKit;
beforeEach(async () => (RDKit = await require('../..')));

test('from and to SMILES', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  expect(mol.toSmiles()).toBe('COCO');
});

test('fromSmiles with replacement', () => {
  const smiles = 'COPheCOPheX';
  const mol = RDKit.Mol.fromSmiles(smiles, {
    replacements: {
      Phe: 'c1ccccc1',
      X: 'N'
    }
  });
  expect(mol.toSmiles()).toBe('COc1ccccc1COc1ccccc1N');
});
