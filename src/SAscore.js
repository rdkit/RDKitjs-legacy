/*
#
# calculation of synthetic accessibility score as described in:
#
# Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
# Peter Ertl and Ansgar Schuffenhauer
# Journal of Cheminformatics 1:8 (2009)
# http://www.jcheminf.com/content/1/1/8
#
# several small modifications to the original paper are included
# particularly slightly different formula for marocyclic penalty
# and taking into account also molecule symmetry (fingerprint density)
#
# for a set of 10k diverse molecules the agreement between the original method
# as implemented in PipelinePilot and this implementation is r2 = 0.97
#
# peter ertl & greg landrum, september 2013
# js adapatation by guillaume godin september 2015
*/


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

min = -4.0;
max = 2.5;

sascore = 11 - (sascore - min + 1) / (max - min) * 9;
console.log(sascore);

if (sascore > 8) 
  	sascore = 8 + math.log(sascore-8);
if (sascore > 10)
   sascore = 10.0;
else if (sascore < 1)
   sascore = 1;

console.log('SAscore:');
console.log(sascore);

/*
#
# calculation of natural product-likeness as described in:
#
# Natural Product-likeness Score and Its Application for Prioritization of Compound Libraries 
# Peter Ertl, Silvio Roggo, and Ansgar Schuffenhauer
# Journal of Chemical Information and Modeling, 48, 68-74 (2008)
# http://pubs.acs.org/doi/abs/10.1021/ci700286x
#
# for the training of this model only openly available data have been used
# ~50,000 natural products collected from various open databases
# ~1 million drug-like molecules from ZINC as a "non-NP background"
#
# peter ertl, august 2015
# javascript adaptation september 2015
*/
var dict = require ('./publicnp.model.json');

npscore = 0;
for (i =0;i<fplen;i++){   
	v = u.get(i+fplen);
	nf = nf + v;
	npscore = npscore + dict[u.get(i)];
}

npscore = npscore/nAtoms;

//# preventing score explosion for exotic molecules
if (npscore > 4) 
    npscore = 4 + Math.log10(npscore - 4 + 1);
if (npscore < -4)
    npscore = -4 - Math.log10(-4. -npscore + 1.);
console.log('NPscore:');
console.log(npscore);





