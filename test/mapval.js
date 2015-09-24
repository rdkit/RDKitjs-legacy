
var RDKit = require('rdkit');

var mol = RDKit.Molecule.fromSmiles('CCCCCC');

//var obj = mol.getMorganFP_getNonzeroElements(2048);	


var u =mol.getMorganFPlist(2);
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