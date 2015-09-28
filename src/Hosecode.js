
var RDKit = require('rdkit');

var smi = 'ClC1(=C(C=CC(=C1)C(=O)C2(C(N)CCCC2))C)';
console.log(smi);

var mol = RDKit.Molecule.fromSmiles(smi);

function readvector(vect) {
	var vectolength =  vect.size();
	var vectorval= [];
	for (i =0;i<vectolength;i++){   
	    vectorval[i] = vect.get(i); 
	}

	return vectorval;
}



//var bondSymbols = { "=" : 200, "%" : 300, "*" :  100};
//var atomicSymbols  = {"C" :90, "O" :89, "N" : 88, "S" : 87, "P" : 86, "Si" : 85, "B" : 84, "F" : 83, "Cl" : 82, "Br" : 81, ";" : 80, "I" : 79, "#" : 21, "&" : 11, ",": 10};


function HosecodeSphere(sphere, atomlistvect, bondtype, ringinfo) {

 for (var j=0;j<sphere.length;j++) {
    var currentsphere = sphere[j];
    var tmphose =[];
    var hosescore =[];
    var spherecode = [];
	  // Determine atoms surrounding each atom and build HOSE code list
    for (var i=0; i<currentsphere.length;i++) {
	  	var  tmp_str = "";
	  	var  hs = 0;
	  	var atomid = currentsphere[i]; 
	  	console.log(atomid,":",atomlistvect[atomid]);
	  	if (bondtype[atomid] === 3) { tmp_str=tmp_str+"%"; hs += 300; }
	  	if (bondtype[atomid] === 1.5) { tmp_str=tmp_str+"*"; hs += 100; }
	  	if (bondtype[atomid] === 2) { tmp_str=tmp_str+"="; hs += 200; }
	  	if (atomlistvect[atomid] === "C") { tmp_str=tmp_str+"C"; hs += 24; }
	  	if (atomlistvect[atomid] === "O") { tmp_str=tmp_str+"O"; hs += 22; }
	  	if (atomlistvect[atomid] === "N") { tmp_str=tmp_str+"N"; hs += 20; }
	  	if (atomlistvect[atomid] === "P") { tmp_str=tmp_str+"P"; hs += 18; }
	  	if (atomlistvect[atomid] === "S") { tmp_str=tmp_str+"S"; hs += 16; }
	  	if (atomlistvect[atomid] === "Si") { tmp_str=tmp_str+"Q"; hs+= 14; }
	  	if (atomlistvect[atomid] === "B") { tmp_str=tmp_str+"B"; hs += 12; }
	  	if (atomlistvect[atomid] === "F") { tmp_str=tmp_str+"F"; hs += 10; }
	  	if (atomlistvect[atomid] === "Cl") { tmp_str=tmp_str+"X"; hs += 8; }
	  	if (atomlistvect[atomid] === "Br") { tmp_str=tmp_str+"Y"; hs += 6; }
	  	if (atomlistvect[atomid] === "I") { tmp_str=tmp_str+"I"; hs += 4; }

	  	tmphose[atomid] = tmp_str;
	  	hosescore[atomid] = hs;
	  	spherecode.push(tmphose,hosescore);
	  }
	  console.log(spherecode);
	}
	return spherecode;
}


function AdjList (graph) {
	 var AL =[];
	 for  (var j = 0; j < graph.length; j += 1) {
	 	var al=[];
        for (var i = 0; i < graph.length; i += 1) {
          if (graph[j][i]>0 && i!=j) { 
          	al.push(i);
 
          }
        }
        AL.push(al);
      }
    
    return AL;
}




function readmat(vect) {
    
	var matsize =  Math.sqrt(vect.size());
	var matrix= [];
	for (i =0;i<matsize;i++){   
		var step = i*matsize;
		var matrixrow = [];
		for (j =0;j<matsize;j++){  
		    matrixrow[j] = vect.get(j+step); 
		}
		matrix[i]=matrixrow;
	}	

	return matrix;
}



// compute thze bfs to store the spheres & distance
function bfs (graph, startNode, len) {
      var parents = [];
      var sphereid = [];
      var spheres = [];
      var bondtype = [];
      var queue = [];
      var visited = [];
      var current;
      queue.push(startNode);
      parents[startNode] = null;
      bondtype[startNode] = null;
      sphereid[startNode]=0;
      visited[startNode] = true;
      while (queue.length) {
        current = queue.shift();
        for (var i = 0; i < graph.length; i += 1) {
          if (i !== current && graph[current][i]>0 && !visited[i]) {  // add >0 to include the bond weight instead of boolean matrix!
            parents[i] = current;
            sphereid[i] = sphereid[current]+1;
            bondtype[i] = graph[current][i];
            visited[i] = true;
            queue.push(i);
          }
        }
      }
    // extract the resulting spheres atomindexes based on the distance to the core! 
    for (var j=1;j<len;j++) {
      	var spherescore = [];

      	for (var i =0;i<graph.length;i+=1){
	      	if(sphereid[i]===j) {
	      		spherescore.push(i); 
	      	}
	    }
		spheres.push(spherescore);

  	}
    return {spheres,parents,bondtype};
}

var numSpheres=5;

var C = mol.getNumAtoms();
console.log("number of Atoms:",C);

var atomneib = mol.getAtomNeighbors(5);	
var arrayatomid =readvector(atomneib);
var bondnum = mol.getBondNeighbors(5);	
var arraybondnum =readvector(bondnum);


console.log("atomids:",JSON.stringify(arrayatomid));
console.log("bondnum:",JSON.stringify(arraybondnum));




var g = mol.getAdjacencyMatrix(true);
var graph =readmat(g);
//var adjlist = AdjList(graph);
console.log(graph);
//console.log(adjlist);

var atomlist = mol.getSymbols(); // getAtomicNum()
var atomlistvect = readvector(atomlist);
console.log("atomelist:",JSON.stringify(atomlistvect));



var rings = mol.findSSSR();
var ringsval= [];

for (i =0;i<rings.size();i++){  
	ringsval[i] = rings.get(i); 
}
console.log("ring list:",ringsval);


console.log("--------------------------");


var  shortestpath = bfs(graph, 5, numSpheres);
console.log("Spheres list:",shortestpath);


var HSC = HosecodeSphere(shortestpath.spheres, atomlistvect, shortestpath.bondtype, ringsval);
console.log("Hose 1:",HSC);

/*
for (j=0;j<C;j++)
{
console.log(j);
/// only scans the carbons!

/*
var atomneib = mol.getAtomNeighbors(j);	
var arrayatomid =readvector(atomneib);
var bondnum = mol.getBondNeighbors(j);	
var arraybondnum =readvector(bondnum);


console.log("atomids:",JSON.stringify(arrayatomid));
console.log("bondnum:",JSON.stringify(arraybondnum));




var  shortestpath = bfs(graph, j, 8);

console.log("Spheres list:",shortestpath);
console.log("..........................");

}
*/
/*







  // Determine atoms surrounding each atom and build HOSE code list
  for (tmp_pt = up.first(); tmp_pt != NULL; tmp_pt = up.next()) {
    if (tmp_pt->element == "C") {  // only look at carbons
	    for (tmp_bond = bonds.first(); tmp_bond != 0; tmp_bond = bonds.next()) {
				if (tmp_bond->Find(tmp_pt) == true) {
				  tmp_str = "";
				  hs = 0;
				  tmp_pt2 = tmp_bond->otherPoint(tmp_pt);
				  if (tmp_bond->Order() == 3) { tmp_str.append("%"); hs += 300; }
				  if ((tmp_pt->aromatic == true) && (tmp_pt2->aromatic == true)) {
				    tmp_str.append("*"); hs += 100;
				  } else {
				    if (tmp_bond->Order() == 2) {
				      	if (tmp_bond->Dash() == 0) { tmp_str.append("="); hs += 200; }
				      	if (tmp_bond->Dash() == 1) { tmp_str.append("*"); hs += 100; }
				    }
				  }
				  if (tmp_pt2->element == "C") { tmp_str.append("C"); hs += 24; }
				  if (tmp_pt2->element == "O") { tmp_str.append("O"); hs += 22; }
				  if (tmp_pt2->element == "N") { tmp_str.append("N"); hs += 20; }
				  if (tmp_pt2->element == "P") { tmp_str.append("P"); hs += 18; }
				  if (tmp_pt2->element == "S") { tmp_str.append("S"); hs += 16; }
				  if (tmp_pt2->element == "Si") { tmp_str.append("Q"); hs+= 14; }
				  if (tmp_pt2->element == "B") { tmp_str.append("B"); hs += 12; }
				  if (tmp_pt2->element == "F") { tmp_str.append("F"); hs += 10; }
				  if (tmp_pt2->element == "Cl") { tmp_str.append("X"); hs += 8; }
				  if (tmp_pt2->element == "Br") { tmp_str.append("Y"); hs += 6; }
				  if (tmp_pt2->element == "I") { tmp_str.append("I"); hs += 4; }
				  if (tmp_pt2->ketos == 2) { tmp_str.append("$$"); }
				  if (tmp_pt2->ketos == 1) { tmp_str.append("$"); }
				  tmp_pt2->tmphose = tmp_str;
				  tmp_pt2->hosescore = hs;
				  sphere1.append(tmp_pt2);
				}
      		}
      tmp_str = "";
      do {
		hsmax = -1;
		for (tmp_pt3 = sphere1.first(); tmp_pt3 != 0; tmp_pt3 = sphere1.next()) {
	  		if (tmp_pt3->hosescore > hsmax) { 
	    			hsmax = tmp_pt3->hosescore; 
	    			ts = tmp_pt3->tmphose;
	    			tmp_pt2 = tmp_pt3;
	  		}
		}
		sphere1.remove(tmp_pt2);
		tmp_str.append(ts);
      	 } while (sphere1.count() > 0);
      tmp_str.append("(//)");
      hosecodes.append(tmp_str);
      tmp_pt->hosecode = tmp_str;
      cout << tmp_str << endl;
    }
  }





}
*/

mol.delete();


