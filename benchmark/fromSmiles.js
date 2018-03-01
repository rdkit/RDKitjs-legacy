var lib = require('..');

lib.then((rdkit) => {

function test() {
  for (var i = 0; i < 100; i++) {
    rdkit.Molecule.fromSmiles('CCOOCCOCOCOCOCCOOCCCCOC');
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
