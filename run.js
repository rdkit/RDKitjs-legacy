/* eslint-disable no-console */

'use strict';

const lib = require('./lib/index');

lib
  .then((rdkit) => {
    var mol = rdkit.Mol.fromSmiles('CCOCOPhe', {
      replacements: { Phe: 'c1ccccc1' }
    });
    console.log(mol.toFASTA());
    console.log(mol.toHELM());
    console.log(mol.toMolBlock());
    console.log(mol.toSmarts());
    console.log(mol.toSmiles());

    var otherMol = rdkit.Mol.fromMolBlock(mol.toMolBlock());
    console.log(otherMol.toSmiles());

    mol.delete();
    otherMol.delete();
  })
  .catch(console.error);
