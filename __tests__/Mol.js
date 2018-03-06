'use strict';

let RDKit;
beforeEach(async () => (RDKit = await require('..')));

test('Mol delete', () => {
  const mol = RDKit.Mol.fromSmiles('CC');
  mol.delete();
  expect(() => mol.toSmiles()).toThrow(/Cannot pass deleted object/);
});
