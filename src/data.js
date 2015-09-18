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

var nAtoms = m.getNumAtoms();
var p = m.FindMolChiralCenters(true);
var nChiralCenters = p.size();
var SBM = m.getSpiroBridgeMacrocycles();
var nBridgeheads = SBM.get(0);
var nSpiro = SBM.get(1);
var nMacrocycles = SBM.get(2);
var sizePenalty = Math.pow(nAtoms,1.005) - nAtoms;
var stereoPenalty = Math.log10(nChiralCenters+1);
var spiroPenalty = Math.log10(nSpiro+1);
var bridgePenalty = Math.log10(nBridgeheads+1);
var macrocyclePenalty = 0;

if (nMacrocycles > 0)
    macrocyclePenalty = math.log10(2);

var score2 = -sizePenalty-stereoPenalty-spiroPenalty-bridgePenalty-macrocyclePenalty;
console.log(score2);

var score3 = 0;
if (nAtoms > fplen)
    score3 = 0.5*Math.log(nAtoms/fplen);

var sascore = score1 + score2 + score3;
console.log(sascore);
min = -4.0;
max = 2.5;
sascore = 11. - (sascore - min + 1) / (max - min) * 9;
if (sascore > 8) 
  	sascore = 8 + math.log(sascore-8);
if (sascore > 10)
   sascore = 10.0;
else if (sascore < 1)
   sascore = 1.0;


console.log(sascore);

