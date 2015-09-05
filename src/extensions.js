var RDKit = require('rdkit');

function calc_all_desc(mol){
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
    
    var t = new RDKit.VectorDouble();
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
    var NumSaturatedRings = mol.NumSaturatedRings();
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
   var slogp = new RDKit.VectorDouble();
    slogp = mol.SlogP_VSA();	
    SlogP_VSA = [];
    for (i =0;i<slogp.size();i++)
    {
    	SlogP_VSA.push(slogp.get(i));
    }


    var smr = new RDKit.VectorDouble();
    smr = mol.SMR_VSA();	
    SMR_VSA = [];
    for (i =0;i<smr.size();i++)
    {
    	SMR_VSA.push(smr.get(i));
    }


    var peo = new RDKit.VectorDouble();
    peo = mol.PEO_VSA();	
    PEO_VSA = [];
    for (i =0;i<peo.size();i++)
    {
    	PEO_VSA.push(peo.get(i));
    }


    var mqn = new RDKit.VectorUint();
    mqn = mol.MQNs();	
    MQNs = [];
    for (i =0;i<mqn.size();i++)
    {
    	MQNs.push(mqn.get(i));
    }
   
   	return {mw:mw,exactMW:ExactMW,formula:Formula,frsp3:fr_sp3,mqn:MQNs,peovsa:PEO_VSA,smrvsa:SMR_VSA, 
   		tpsa:TPSA,slogpvsa:SlogP_VSA,logp:logp,mr:mr,labuteASA:LabuteASA,lipinskiHBD:LipinskiHBD,lipinskiHBA:LipinskiHBA,
   		numHeterocycles:NumHeterocycles,numRings:NumRings,numHeteroatoms:NumHeteroatoms,numHBA:NumHBA,
   		numHBD:NumHBD,numAliphaticCarbocycles:NumAliphaticCarbocycles,numAliphaticHeterocycles:NumAliphaticHeterocycles,
   		numSaturatedCarbocycles:NumSaturatedCarbocycles,numSaturatedHeterocycles:NumSaturatedHeterocycles,
   		numAliphaticRings:NumAliphaticRings,numAromaticRings:NumAromaticRings,numSaturatedRings:NumSaturatedRings,
   		numAmideBonds:NumAmideBonds,numRotatableBonds:NumRotatableBonds,numAromaticHeterocycles:NumAromaticHeterocycles,
   		chi0n:Chi0n,chi1n:Chi1n,chi2n:Chi2n,chi3n:Chi3n,chi4n:Chi4n,chi0v:Chi0v,chi1v:Chi1v,chi2v:Chi2v,chi3v:Chi3v,chi4v:Chi4v,
        kappa1:Kappa1,kappa2:Kappa2,kappa3:Kappa3,hallKierAlpha:HallKierAlpha
   }; 
};

/*
var smi = 'CCCCCOC';
var mol = RDKit.Molecule.fromSmiles(smi);
var c = calc_all_desc(mol);
console.log(c);
mol.delete()
*/


/*
var m = RDKit.Molecule.fromSmiles('CC(C)=C(C)C(=O)c1ccccc1')
m.sdwrite()
m.addHs(); 
//m.sdwrite();
m.EmbedMultipleConfs(); //=> generate Conformers
//m.EmbedMultipleConfsarg(5,100,1)

console.log(m.sdwriteConfs())
var d =m.MMFFOptimizeMoleculeConfs(1,500,'MMFF94'); //=> gave the Conformers optimization result in d!
console.log('AlignMolConformers')

//m.getConformerRMS(1,2);

m.AlignMolConformers();
//var t = m.AlignMolConformersRMSlist();
//console.log(t.size());

for (i=0;i<d.size();i++)
{
    console.log(d.get(i));

}
*/

/*
debugger;

console.log('AlignMol 2')

var m1 = RDKit.Molecule.fromSmiles('CC(C)=C(C)C(=O)C')
m1.addHs();
m1.EmbedMolecule();
m1.MMFFoptimizeMolecule();
var t = m1.AlignMol('CC(C)=C(C)C(=O)CC');
*/


//m.sdwriteConfs()
//m.sdwrite()
//m.AlignMolConformers();
//m.removeHs();
//var p = m.sdwriteConfs();
//m.sdwriteConfs() //=> gave the 3D conformers list in a SDF format

m.delete()


/*

var m = RDKit.Molecule.fromSmiles('CC(C)=C(C)C(=O)c1ccccc1')
m.sdwrite()
m.addHs(); 
//m.sdwrite();
var c = m.EmbedMultipleConfs(); //=> generate Conformers
//console.log(m.sdwriteConfs())
console.log(c)

var p= m.getConformer(0);

console.log(p)

//var p = m.sdwriteConfs();
var d = m.MMFFOptimizeMoleculeConfs(1,100,'MMFF94'); //=> gave the Conformers optimization result in d!
console.log(d)

//console.log(m.sdwriteConfs())
//m.sdwrite()
var e = m.AlignMolConformers();
console.log(e)




for (i=0;i<d.length;i++)
{
    console.log(d.get(i));

}
//m.removeHs();
//var p = m.sdwriteConfs();
m.sdwriteConfs() //=> gave the 3D conformers list in a SDF format

m.delete()


*/








