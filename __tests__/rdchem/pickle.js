'use strict';

const RDKit = require('../..');

beforeEach(() => RDKit.load());

test('from and to pickle', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  const pickle = mol.toPickle();

  expect(typeof pickle).toBe('string');
  console.log('before new mol');

  const parsed = RDKit.Mol.fromPickle(pickle);
  console.log('after from Pickle');

  expect(parsed.toSmiles()).toBe('COCO');
});
