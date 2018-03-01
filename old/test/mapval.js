
var RDKit = require('rdkit');

var mol = RDKit.Molecule.smilesToMol('CCCCCC');

//var obj = mol.getMorganFingerprints_getNonzeroElements(2048);	


var u =mol.getMorganFingerprintslist(2);
        var fplen =  u.size();
        var v= [];
        for (i =0;i<fplen;i++){   
            v[i] = u.get(i); 

        }

console.log(u);
/*
for (var prop in obj) {
  console.log("obj." + prop + " = " + obj[prop]);
}
*/