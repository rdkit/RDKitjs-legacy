'use strict';

let RDKit;
beforeEach(async () => (RDKit = await require('../..')));

test('embedMultipleConfs', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  mol.addHs();
  const result = mol.embedMultipleConfs({ numConfs: 5 });
  expect(result).toEqual([0, 1, 2, 3, 4]);

  for (const confId of result) {
    mol.MMFFOptimizeMolecule({ confId });
  }

  mol.alignMolConformers();

  const molfiles = result.map((confId) => mol.toMolBlock({ confId }));
  expect(molfiles).toMatchSnapshot();
});
