'use strict';

const RDKit = require('../..');

beforeEach(() => RDKit.load());

test('default parameters', () => {
  const smiles = 'COCO';
  const mol = RDKit.Mol.fromSmiles(smiles);
  mol.addHs();
  mol.embedMolecule();
  const result = mol.MMFFOptimizeMolecule();
  expect(result).toHaveLength(2);
  expect(result[0]).toBe(0);
  expect(result[1]).toBeCloseTo(-19.4, 1);
  expect(mol.toMolBlock()).toMatchSnapshot();
});
