'use strict';

var RDKit = require('..');

beforeEach(() => RDKit.load());

test('set & getProp', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.setProp('prop1', 'rdkitjstest');
  expect(mol.getProp('prop1')).toBe('rdkitjstest');
});

test('hasProp', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.setProp('prop1', 'rdkitjstest');
  expect(mol.hasProp('prop1')).toBe(true);
  expect(mol.hasProp('prop2')).toBe(false);
});

test.skip('getproplist', function() {
  var smi = 'CCCCCOC(CO)';
  var mol = RDKit.Molecule.fromSmiles(smi);
  mol.setProp('prop1', 'rdkitjstest');
  expect(mol.getproplist()).toBe('');
});
