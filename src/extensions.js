

function calc_all_desc( smi ){
    var rdkit = require('rdkit');
	var mol = rdkit.Molecule.fromSmiles( smi );
	var mw = mol.getMW();
	var fr_sp3 = mol.FractionCSP3();
    
    var ExactMW= mol.ExactMW();
    var Formula= mol.Formula();
    var Chi0v= mol.Chi0v();
    var Chi1v= mol.Chi1v();
    var Chi2v= mol.Chi2v();
    var Chi3v= mol.Chi3v();
    var Chi4v= mol.Chi4v();
    var Chi0n= mol.Chi0n();
    var Chi1n= mol.Chi1n();
    var Chi2n= mol.Chi2n();
    var Chi3n= mol.Chi3n();
    var Chi4n= mol.Chi4n();
    var HallKierAlpha= mol.HallKierAlpha();
    var Kappa1= mol.Kappa1();
    var Kappa2= mol.Kappa2();
    var Kappa3= mol.Kappa3();
    /*
    var p ={};
    p=mol.logp_mr()[0]; 
    console.log(p);
    var logp=p[0];
    var mr=p[1];
    console.log(mr);
    console.log(logp);
    */
    var LipinskiHBA=mol.LipinskiHBA();
    var LipinskiHBD=mol.LipinskiHBD();
    var NumRotatableBonds=mol.NumRotatableBonds();
    var NumHBD=mol.NumHBD();
    var NumHBA=mol.NumHBA();
    var NumHeteroatoms=mol.NumHeteroatoms();
    var NumAmideBonds=mol.NumAmideBonds();
    var NumRings=mol.NumRings();
    var NumAromaticRings=mol.NumAromaticRings();
    var NumAliphaticRings=mol.NumAliphaticRings();
    /*
    var NumSaturatedRings=mol.NumSaturatedRings(); // not working ???
    */
    var NumHeterocycles=mol.NumHeterocycles();
    var NumAromaticHeterocycles=mol.NumAromaticHeterocycles();
    var NumAromaticCarbocycles=mol.NumAromaticCarbocycles();
    var NumSaturatedHeterocycles=mol.NumSaturatedHeterocycles();
    var NumSaturatedCarbocycles=mol.NumSaturatedCarbocycles();
    var NumAliphaticHeterocycles=mol.NumAliphaticHeterocycles();
    var NumAliphaticCarbocycles=mol.NumAliphaticCarbocycles();
    var LabuteASA= mol.LabuteASA();
    var TPSA= mol.TPSA();
    var SlogP_VSA= mol.SlogP_VSA();
    var SMR_VSA= mol.SMR_VSA();
    var PEO_VSA= mol.PEO_VSA();
    var MQN =mol.MQNs();	
   	
   	return { mw:mw, fr_sp3: fr_sp3,smi: smi };
};


var c = calc_all_desc('CCCCCOC');
console.log(c)


