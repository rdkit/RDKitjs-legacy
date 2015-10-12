// Prediction of response factors for gas chromatography with flame ionization detection
var RDKit = require('rdkit');

function fragnum(smarts, functions, mol) {
  var desc = [];
  for (i=0; i<smarts.length;i++) {
     var v = mol.GetSubstructMatchesNumber(smarts[i]);
  	 //console.log(smarts[i],functions[i],v);

     desc[i]=v;
 }
  return desc
}

// from Jean-Yves article: http://onlinelibrary.wiley.com/store/10.1002/jssc.201500106/asset/jssc4443.pdf?v=1&t=if6srj7k&s=3454006c4a3e61decd6e4d9c46dd5b651f0962a2
function rrf(smi) {
    var mol = RDKit.Molecule.fromSmiles(smi);  
    mol.addHs(); // need to add the hydrogens for the count of H!
    var Benz_count= mol.GetSubstructMatchesNumber('c1ccccc1'); // number of benzene
    var functionHoC = ['C_count','H_count','O_count', 'N_count','S_count','F_count','Cl_count','Br_count','I_count','Si_count']; // function of HoC
    var smartHoC = ['[C,c]','[H]','[O,o]','[N,n]','[S,s]','F','Cl','Br','I','[Si]']; // smarts codes
    var HoCcoef =  [103.57, 21.85, -48.18, 7.46, 74.67, -23.57, -27.43, -10.9, -2.04, 46.5]; // HoC coefs
	var HoCfinger =fragnum(smartHoC, functionHoC, mol); // get the HoC fingerprint
	//console.log('nb_benz:',Benz_count);
    
	var Heat_of_combustion = 11.06; // initialize the HoC 

	// compute the HoC
	for (j=0;j<HoCfinger.length;j++) 
	{ 
  		if (HoCfinger[j]>0) { Heat_of_combustion+= HoCcoef[j]*HoCfinger[j];}
	}

	// get the rfm
	var rfmolecular = - 0.0708  + 8.57 * 0.0001 * Heat_of_combustion + 1.27 * 0.1 * Benz_count + 6.18 * 0.01 * HoCfinger[7]; // brcount

	 // mw = molecular_weigth;  not sure if exactMW or MW ???
	var mw =mol.ExactMW();	

	 // get the rrf
    var rrf = mw*(1/rfmolecular)/158.24;


return rrf;

}
 
var smi = 'BrC(CCl)C(CCC)CCOC(CO)c1ccc(F)cc1';
console.log('rrf',rrf(smi));