'use strict';

let RDKit;
beforeEach(async () => (RDKit = await require('../..')));

test('addHs', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  mol.addHs();
  expect(mol.toSmiles()).toBe('[H]OC([H])([H])OC([H])([H])[H]');
});
