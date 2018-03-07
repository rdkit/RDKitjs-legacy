/* eslint-disable no-console */

'use strict';

const rdkit = require('./lib/index');

rdkit
  .load()
  .then(() => {
    var mol = rdkit.Mol.fromSmiles('CCOCOPhe', {
      replacements: { Phe: 'c1ccccc1' }
    });

    mol.addHs();

    console.log(mol.toFASTA());
    console.log(mol.toHELM());
    console.log(mol.toMolBlock());
    console.log(mol.toSmarts());
    console.log(mol.toSmiles());

    var otherMol = rdkit.Mol.fromMolBlock(mol.toMolBlock());
    console.log(otherMol.toSmiles());

    const pickle = otherMol.toPickle();
    console.log(pickle);

    var thirdMol = rdkit.Mol.fromPickle(pickle);
    console.log(thirdMol.toSmiles());

    mol.delete();
    otherMol.delete();
    thirdMol.delete();
  })
  .catch(console.error);
