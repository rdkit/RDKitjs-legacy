'use strict';

let RDKit = require('..');

beforeEach(async () => {
  RDKit = await RDKit;
});

test('smilesToMolfile', () => {
  const smiles = 'COCO';
  const molfile = RDKit.smilesToMolfile(smiles);
  expect(molfile).toMatchSnapshot();
  expect(molfile).toContain('4  3  0  0  0  0  0  0  0  0999 V2000');
});
