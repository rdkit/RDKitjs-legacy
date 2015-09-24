
var RDKit = require('rdkit');

var mol = RDKit.Molecule.fromSmiles('CCCCCC');

var obj = mol.getMorganFP_getNonzeroElements(2);	

var u =mol.getMorganFPlist(2);
var fplen =  u.size();
var v= [];
for (i =0;i<fplen;i++){   
    v[i] = u.get(i); 
	try{
            console.log(obj.get(u.get(i)));
 	}
	catch (err){	
	}
}

console.log(JSON.stringify(v));

for (var prop in obj) {
  console.log("obj." + prop + " = " + obj[prop]);
}
