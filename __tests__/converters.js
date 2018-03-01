'use strict';

const RDKit = require('..');

beforeEach(() => RDKit);

test('smilesToMolfile', () => {
  const smiles = 'COCO';
  const molfile = RDKit.smilesToMolfile(smiles);
  expect(molfile).toMatchSnapshot();
});
