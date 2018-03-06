'use strict';

let RDKit;
beforeEach(async () => (RDKit = await require('..')));

test('smiles to molBlock', () => {
  const smiles = 'COCO';
  const molfile = RDKit.Mol.fromSmiles(smiles).toMolBlock();
  expect(molfile).toMatchSnapshot();
  expect(molfile).toContain('4  3  0  0  0  0  0  0  0  0999 V2000');
});
