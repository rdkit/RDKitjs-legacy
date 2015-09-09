

var dict2 = require ('./fpscores.pkldico.json');


//console.log(Object.keys(dict2).length);


var dict = require ('./publicnp.model.json');


var RDKit = require ('rdkit');
var R = RDKit;



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

function SAscore (dico, nAtoms,fplen, morganFPlist, molecule){
	var nf = 0;
	var score1=0;
	var v=0;
	var f=0;
	for (i =0;i<fplen;i++){   
		v = morganFPlist.get(i+fplen);
		f = dico[morganFPlist.get(i)];
		if (f !== undefined) {  
			score1 = score1 + f*v;
		}
		// default value for undefined fragments!
		else {
			score1 = score1 -4*v;
		} 
		nf = nf + v;
	}

	score1=score1/nf;

	console.log(score1);

	var p = molecule.FindMolChiralCenters(true);
	var nChiralCenters = p.size();
	var SBM = molecule.getSpiroBridgeMacrocycles(); // there is a issue there S
	var nSpiro = SBM.get(0);
	var nBridgeheads = SBM.get(1);
	var nMacrocycles = SBM.get(2);
	var sizePenalty = Math.pow(nAtoms,1.005) - nAtoms;
	var stereoPenalty = Math.log10(nChiralCenters+1);
	var spiroPenalty = Math.log10(nSpiro+1);
	var bridgePenalty = Math.log10(nBridgeheads+1);
	var macrocyclePenalty = 0;

  	console.log("1:",nBridgeheads,"2:",nMacrocycles,"3:",nSpiro,"4:",nChiralCenters,"5:",macrocyclePenalty,"6:",sizePenalty)


	if (nMacrocycles > 0)
	    macrocyclePenalty = math.log10(2);
    
	var score2 = 0-sizePenalty-stereoPenalty-spiroPenalty-bridgePenalty-macrocyclePenalty;

	var score3 = 0;
	if (nAtoms > fplen)
	    score3 = 0.5*Math.log(nAtoms/fplen);

	var sascore = score1 + score2 + score3;

	min = -4.0;
	max = 2.5;

	sascore = 11.0 - (sascore - min + 1) / (max - min) * 9;

	if (sascore > 8.0) 
	  	sascore = 8.0 + math.log(sascore-8.0);
	if (sascore > 10.0)
	   sascore = 10.0;
	else if (sascore < 1.0)
	   sascore = 1.0;
	console.log("1:",score1,"2:",score2,"3:",score3,"final:",sascore);
	console.log("-----------------------------------------------------");

	return sascore;
}


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
function NPscore(dicts, nAtoms,fplen,morganFPlist) {
    var v=0;
    var f=0;
	var npscore = 0;
	for (i =0;i<fplen;i++){   
		v = morganFPlist.get(i+fplen);
		f =dicts[morganFPlist.get(i)];
		if (f !== undefined)
			npscore = npscore + f;
	}
	npscore = npscore/nAtoms;

	//# preventing score explosion for exotic molecules
	if (npscore > 4) 
	    npscore = 4 + Math.log10(npscore - 4 + 1);
	if (npscore < -4)
	    npscore = -4 - Math.log10(-4. -npscore + 1.);
	return npscore;
}


//loop over file or list of smiles!

var fs = require('fs');
var np = 0;
var sa = 0;
fs.readFileSync('test-sardkit.txt').toString().split(/\r?\n/).forEach(function(line){
    var linecomponent = line.split('\t'); 
    try {

		var m=R.Molecule.fromSmiles(linecomponent[0]);
		if (m !==undefined)
		{
		var morganFPlist =m.getMorganFPlist(2);
		var nAtoms = m.getNumAtoms();
		var fplen =  morganFPlist.size()/2;
		sa = SAscore (dict2,nAtoms,fplen, morganFPlist,m);
		np = NPscore(dict, nAtoms,fplen,morganFPlist);

		console.log(linecomponent[1],";sa:",sa.toPrecision(4),";np:",np.toPrecision(4),";",linecomponent[0]);

		m.delete();
	}
	}
	catch (err) {
            console.log(linecomponent);
            console.log(err);

	}

});


/*
var m=R.Molecule.fromSmiles('CC(=O)OC1C[NH+]2CCC1CC2');
var morganFPlist =m.getMorganFPlist(2);
var nAtoms = m.getNumAtoms();
var fplen =  morganFPlist.size()/2;
sa = SAscore (dict2,nAtoms,fplen, morganFPlist,m);
console.log("sa:",sa.toPrecision(4));
*/







