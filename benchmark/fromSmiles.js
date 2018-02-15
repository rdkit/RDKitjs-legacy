var rdkit = require('..');

rdkit.load().then(() => {

function test() {
  for (var i = 0; i < 100; i++) {
    rdkit.Molecule.smilesToMol('CCOOCCOCOCOCOCCOOCCCCOC');
  }
}

test();

console.time('test');
test();
test();
test();
test();
test();
console.timeEnd('test');

});
