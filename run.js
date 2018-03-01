'use strict';

const lib = require('./lib/index');

lib.then((rdkit) => {
  console.log(rdkit);
  console.log(rdkit.smilesToMolfile('COCO'));

  const mol = rdkit.Molecule.fromSmiles('COCO');
  console.log(mol.toMolfile());

  console.log('ici');
  console.log(rdkit.Molecule.smilesTo3D('COCO').toMolfile());

  console.log(mol.addHs());
  console.log('la');

  console.log(mol.EmbedMolecule());
  console.log('encore la');

  console.log(mol.MMFFoptimizeMolecule());
  console.log('enfin la');

  console.log(mol.toMolfile());


   

  mol.delete();

}).catch(console.error);
