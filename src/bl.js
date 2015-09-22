var RDKit = require('RDKit');


var mol = RDKit.Molecule.fromSmiles('COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21');

var d = mol.getAtomCode(0);

var e = mol.getHashedAtomPairFingerprintAsBitVect(2048,0,1);

var f = mol.getAtomPairFingerprint();

console.log(d);   
 
console.log(e);   

for (i=0;i<f.size();i++)
 {
    console.log(f.get(i));
 }   
