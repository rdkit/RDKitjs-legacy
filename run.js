'use strict';

const lib = require('./lib/index');

lib.then((rdkit) => {
  console.log(rdkit);
  console.log(rdkit.smilesToMolfile('COCO'));

  const mol = rdkit.Molecule.fromSmiles('COCO');
  console.log(mol.toMolfile());
  mol.delete();
});
