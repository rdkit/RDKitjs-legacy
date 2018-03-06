'use strict';

let RDKit;
beforeEach(async () => (RDKit = await require('../..')));

test('from and to pickle', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);

  const pickle = mol.toPickle();
  expect(typeof pickle).toBe('string');

  const parsed = RDKit.Mol.fromPickle(pickle);
  expect(parsed.toSmiles()).toBe('COCO');
});
