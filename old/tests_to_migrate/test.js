'use strict';

var RDKit = require('..');

beforeEach(() => RDKit.load());

// atom & bond manipulations
describe('new molecule using atom bond', function() {
  it('create molecule and add atoms and bonds', function() {
    var t = RDKit.Molecule.newmolecule();
    t.addAtom(6);
    t.addAtom(6);
    t.addAtom(6);
    t.addAtom(6);
    t.addAtom(6);
    t.addAtom(8);
    t.addBond(0, 1, 1);
    t.addBond(1, 2, 2);
    t.addBond(2, 3, 1);
    t.addBond(3, 4, 1);
    t.addBond(4, 5, 1);
    expect(t.molToSmiles()).toBe('CC=CCCO');
    t.delete();
  });
  // this is in development stage caution not working for the moment!!!!

  it.skip('setBondDir', function() {
    setBondDir(Bondid, bonddirid);
  });
});

describe.skip('getPath', function() {
  it('should work', function() {
    var p = RDKit.molecule.getPath();
    expect(p).toBe('');
    p.delete();
  });
});

describe.skip('findSSSR', function() {
  it('should work', function() {
    findSSSR(res);
  });
});
