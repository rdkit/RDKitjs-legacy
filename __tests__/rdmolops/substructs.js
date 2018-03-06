'use strict';

let RDKit;
beforeEach(async () => (RDKit = await require('../..')));

test('delete substructure from smarts', () => {
  const m = RDKit.Mol.fromSmiles('CC(=O)O');
  const patt = RDKit.Mol.fromSmarts('C(=O)[OH]');
  const result = m.deleteSubstructs(patt);
  expect(result.toSmiles()).toBe('C');
});
