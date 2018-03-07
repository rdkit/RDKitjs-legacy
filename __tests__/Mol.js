'use strict';

const RDKit = require('..');

beforeEach(() => RDKit.load());

test('Mol delete', () => {
  const mol = RDKit.Mol.fromSmiles('CC');
  mol.delete();
  expect(() => mol.toSmiles()).toThrow(/Cannot pass deleted object/);
});
