'use strict';

let RDKit;
beforeEach(async () => RDKit = await require('..'));

test('smilesToMolfile', () => {
  const smiles = 'COCO';
  const molfile = RDKit.smilesToMolfile(smiles);
  expect(molfile).toMatchSnapshot();
  expect(molfile).toContain('4  3  0  0  0  0  0  0  0  0999 V2000');
});
