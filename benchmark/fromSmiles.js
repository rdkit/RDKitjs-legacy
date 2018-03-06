'use strict';

require('..').then((rdkit) => {
  function test() {
    for (var i = 0; i < 100; i++) {
      rdkit.Mol.fromSmiles('CCOOCCOCOCOCOCCOOCCCCOC');
    }
  }

  test();

  console.time('parse a SMILES 500 times');
  test();
  test();
  test();
  test();
  test();
  console.timeEnd('parse a SMILES 500 times');
});
