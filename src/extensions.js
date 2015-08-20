var rdkit = require('rdkit');

function calc_all_desc(smi){
	var mol = rdkit.Molecule.fromSmiles(smi);
	var mw = mol.getMW();
	var fr_sp3 = mol.FractionCSP3();
    var ExactMW = mol.ExactMW();
    var Formula = mol.Formula();
    var Chi0v = mol.Chi0v();
    var Chi1v = mol.Chi1v();
    var Chi2v = mol.Chi2v();
    var Chi3v = mol.Chi3v();
    var Chi4v = mol.Chi4v();
    var Chi0n = mol.Chi0n();
    var Chi1n = mol.Chi1n();
    var Chi2n = mol.Chi2n();
    var Chi3n = mol.Chi3n();
    var Chi4n = mol.Chi4n();
    var HallKierAlpha = mol.HallKierAlpha();
    var Kappa1 = mol.Kappa1();
    var Kappa2 = mol.Kappa2();
    var Kappa3 = mol.Kappa3();
    
    var t = new rdkit.VectorDouble();
    t = mol.logp_mr(); 
    var logp = t.get(0);
    var mr = t.get(1);

    var LipinskiHBA = mol.LipinskiHBA();
    var LipinskiHBD = mol.LipinskiHBD();
    var NumRotatableBonds = mol.NumRotatableBonds();
    var NumHBD = mol.NumHBD();
    var NumHBA = mol.NumHBA();
    var NumHeteroatoms = mol.NumHeteroatoms();
    var NumAmideBonds = mol.NumAmideBonds();
    var NumRings = mol.NumRings();
    var NumAromaticRings = mol.NumAromaticRings();
    var NumAliphaticRings = mol.NumAliphaticRings();
    //var NumSaturatedRings = mol.NumSaturatedRings(); // not working ???
    var NumHeterocycles = mol.NumHeterocycles();
    var NumAromaticHeterocycles = mol.NumAromaticHeterocycles();
    var NumAromaticCarbocycles = mol.NumAromaticCarbocycles();
    var NumSaturatedHeterocycles = mol.NumSaturatedHeterocycles();
    var NumSaturatedCarbocycles = mol.NumSaturatedCarbocycles();
    var NumAliphaticHeterocycles = mol.NumAliphaticHeterocycles();
    var NumAliphaticCarbocycles = mol.NumAliphaticCarbocycles();
    var LabuteASA = mol.LabuteASA();
    var TPSA = mol.TPSA();

    // get array values
   var slogp = new rdkit.VectorDouble();
    slogp = mol.SlogP_VSA();	
    SlogP_VSA = [];
    for (i =0;i<slogp.size();i++)
    {
    	SlogP_VSA.push(slogp.get(i));
    }


    var smr = new rdkit.VectorDouble();
    smr = mol.SMR_VSA();	
    SMR_VSA = [];
    for (i =0;i<smr.size();i++)
    {
    	SMR_VSA.push(smr.get(i));
    }


    var peo = new rdkit.VectorDouble();
    peo = mol.PEO_VSA();	
    PEO_VSA = [];
    for (i =0;i<peo.size();i++)
    {
    	PEO_VSA.push(peo.get(i));
    }


    var mqn = new rdkit.VectorUint();
    mqn = mol.MQNs();	
    MQNs = [];
    for (i =0;i<mqn.size();i++)
    {
    	MQNs.push(mqn.get(i));
    }
   
   	return {smi: smi, mw:mw, ExactMW:ExactMW,Formula:Formula, fr_sp3: fr_sp3,MQNs:MQNs,PEO_VSA:PEO_VSA,SMR_VSA:SMR_VSA, 
   		TPSA:TPSA,SlogP_VSA:SlogP_VSA,logp:logp,mr:mr, LabuteASA:LabuteASA, LipinskiHBD:LipinskiHBD, LipinskiHBA:LipinskiHBA,
   		NumHeterocycles:NumHeterocycles,NumRings:NumRings,NumHeteroatoms:NumHeteroatoms,NumHBA:NumHBA,
   		NumHBD:NumHBD,NumAliphaticCarbocycles:NumAliphaticCarbocycles,NumAliphaticHeterocycles:NumAliphaticHeterocycles,
   		NumSaturatedCarbocycles:NumSaturatedCarbocycles,NumSaturatedHeterocycles:NumSaturatedHeterocycles,
   		NumAliphaticRings:NumAliphaticRings, NumAromaticRings:NumAromaticRings, NumAmideBonds:NumAmideBonds,
   		NumRotatableBonds:NumRotatableBonds,NumAromaticHeterocycles:NumAromaticHeterocycles,
   		Chi0n:Chi0n,Chi1n:Chi1n,Chi2n:Chi2n,Chi3n:Chi3n,Chi4n:Chi4n,
   		Chi0v:Chi0v, Chi1v:Chi1v, Chi2v:Chi2v, Chi3v:Chi3v, Chi4v:Chi4v, Kappa1:Kappa1,Kappa2:Kappa2, Kappa3:Kappa3,
   		HallKierAlpha:HallKierAlpha
   }; 
};


var c = calc_all_desc('CCCCCOC');
console.log(c)


