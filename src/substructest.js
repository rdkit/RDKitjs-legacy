
var RDKit = require ('rdkit');
var R = RDKit;


mol = R.Molecule.fromSmiles('c1ccccc1O');
var p = mol.GetSubstructMatches('ccO');

for (i=0;i<p.size();i++){

	console.log(p.get(i));
} 

