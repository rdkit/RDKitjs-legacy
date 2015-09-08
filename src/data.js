//var dict = require ('./publicnp.model.json');

//console.log(Object.keys(dict).length);


//var result3 = dict['1449131080']
//console.log(result3);

var dict2 = require ('./fpscores.pkldico.json');

console.log(Object.keys(dict2).length);


var RDKit = require ('rdkit');
var R = RDKit;
m=R.Molecule.fromSmiles('Cc1ccccc1');
var u =m.getMorganFPlist(2);
var fplen =  u.size()/2;
var nf = 0;
var score1=0;
var v=0;
for (i =0;i<fplen;i++){   
	v = u.get(i+fplen);
	nf = nf + v;
	score1 = score1 + dict2[u.get(i)]*v;
}
score1=score1/nf;

console.log(score1);
