'use strict';

const RDKit = require('../..');

beforeEach(() => RDKit.load());

test('toMolBlock', () => {
  const smiles = 'COCO';
  const molfile = RDKit.Mol.fromSmiles(smiles).toMolBlock();
  expect(molfile).toMatchSnapshot();
  expect(molfile).toContain('4  3  0  0  0  0  0  0  0  0999 V2000');
});
