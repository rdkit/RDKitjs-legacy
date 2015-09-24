
var RDKit = require('rdkit');

var mol = RDKit.Molecule.fromSmiles('C1(=C)CC#CC1');

function readvector(vect) {
	var vectolength =  vect.size();
	var vectorval= [];
	for (i =0;i<vectolength;i++){   
	    vectorval[i] = vect.get(i); 
	}

	return vectorval;
}


var bondSymbols = { "=" : 200000, "%" : 300000, "*" :  100000};
var atomicSymbols  = {"C" :9000, "O" :8900, "N" : 8800, "S" : 8700, "P" : 8600, "Si" : 8500, "B" : 8400, "F" : 8300, "Cl" : 8200, "Br" : 8100, ";" : 8000, "I" : 7900, "#" : 2100, "&" : 1100, ",": 1000};

var numSphere = 4;

/// only scans the carbons!

var atomneib = mol.getAtomNeighbors(0);	
var arrayatomid =readvector(atomneib);
var bondnum = mol.getBondNeighbors(0);	
var arraybondnum =readvector(bondnum);


var g = mol.getAdjacencyMatrix(true);
var garray =readvector(g);

console.log("grapharray",garray);
var C = mol.getNumAtoms();
console.log("C:",C);
console.log("atomids:",JSON.stringify(arrayatomid));
console.log("bondnum:",JSON.stringify(arraybondnum));

/*
for (int s;s<numSphere;s++) {
	for (int i;i<arrayatomid.length;i++) {
		var next = mol.getAtomNeighbors(arrayatomid[i]);	
		var arrayatomid =readvector(next);
	}

  }

}
*/


mol.delete();


