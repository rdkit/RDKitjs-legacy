
var RDKit = require('rdkit');

var smi = 'ClC1(=C(C=CC(=C1)C(=O)C2(C(N)CCCC2))C)';
var smi1 = 'C1(=C(C=CC(=C1)C(=O)C2(C(N)CCCC2))C)';

console.log(smi);

var mol = RDKit.Molecule.fromSmiles(smi);

var obj2 = mol.getMorganUFPkeys(2);

var u =mol.getMorganFPlist(2);


console.log(readvector(obj2));


console.log(readvector(u));

var mol1 = RDKit.Molecule.fromSmiles(smi1);

console.log(mol.Similarity(mol1,'dice',0,0));

// gave the canonical order of atoms (new algo!)
var UniquerankAtoms = mol.CanonicalRankAtoms(true,true,true);

console.log(readvector(UniquerankAtoms));


function Node(value) {
    this.value = value;
    this.children = [];
    this.parent = null;

    this.setParentNode = function(node) {
        this.parent = node;
    }

    this.getParentNode = function() {
        return this.parent;
    }

    this.addChild = function(node) {
        node.setParentNode(this);
        this.children[this.children.length] = node;
    }

    this.getChildren = function() {
        return this.children;
    }

    this.removeChildren = function() {
        this.children = [];
    }
}


function HOSEsphere(value) {
    this.value = value;
    this.children = [];
    this.parent = null;

    this.setParentNode = function(node) {
        this.parent = node;
    }

    this.getParentNode = function() {
        return this.parent;
    }

    this.addChild = function(node) {
        node.setParentNode(this);
        this.children[this.children.length] = node;
    }

    this.getChildren = function() {
        return this.children;
    }

    this.removeChildren = function() {
        this.children = [];
    }
}





function readvector(vect) {
	var vectolength =  vect.size();
	var vectorval= [];
	for (i =0;i<vectolength;i++){   
	    vectorval[i] = vect.get(i); 
	}
	return vectorval;
}


function diff(a,b) {
	var nar = [];
	var ar = a.concat(b);

	for (var i in ar) {
		var f =ar[i];
		var t = 0;
		for (j=0; j<ar.length; j++) {
			if (ar[j]=== f) {
				t++;
			}
		}
		if (t===1) {
			nar.push(f);
		}

	}
	return nar;
};



//var atomicSymbols  = {"C" :90, "O" :89, "N" : 88, "S" : 87, "P" : 86, "Si" : 85, "B" : 84, "F" : 83, "Cl" : 82, "Br" : 81, ";" : 80, "I" : 79, "&" : 11, ",": 10};

function HC(Symbol,bondtype,isdejavu){
	  	var  tmp_str = "";
	  	var  hs = 0;

	  	if (isdejavu)  {Symbol = "&";}

		if (bondtype=== 3) { tmp_str=tmp_str+"%"; hs += 300; }
	  	if (bondtype=== 1.5) { tmp_str=tmp_str+"*"; hs += 100; }
	  	if (bondtype=== 2) { tmp_str=tmp_str+"="; hs += 200; }
	  	if (Symbol === "C") { tmp_str=tmp_str+"C"; hs += 90; }
	  	if (Symbol === "O") { tmp_str=tmp_str+"O"; hs += 89; }
	  	if (Symbol === "N") { tmp_str=tmp_str+"N"; hs += 88; }
	  	if (Symbol === "S") { tmp_str=tmp_str+"S"; hs += 87; }
	  	if (Symbol === "P") { tmp_str=tmp_str+"P"; hs += 86; }
	  	if (Symbol === "Si") { tmp_str=tmp_str+"Q"; hs+= 85; }
	  	if (Symbol === "B") { tmp_str=tmp_str+"B"; hs += 84; }
	  	if (Symbol === "F") { tmp_str=tmp_str+"F"; hs += 83; }
	  	if (Symbol === "Cl") { tmp_str=tmp_str+"X"; hs += 82;}
	  	if (Symbol === "Br") { tmp_str=tmp_str+"Y"; hs += 81;}
	  	if (Symbol === "I") { tmp_str=tmp_str+"I"; hs += 79; }
	  	if (Symbol === "&") { tmp_str=tmp_str+"&"; hs += 11; }

	  	if (!Symbol | Symbol ===null) {tmp_str = ","; hs = 10;}

	return {hs,tmp_str};

}


function HosecodeSphere(sphere, atomlistvect, bondtype, ringinfo, Neiarray, parents) { 
 var sphereScore=[];
 var child =childs(sphere, Neiarray, parents);
 var previousGroup=[];
 var spherecode = [];
 var dejavu = [];
 var S=[];

  dejavu[5]=true;

 for (var j=0;j<sphere.length;j++) {
    var currentsphere = sphere[j];
    console.log(currentsphere);
    var tmphose =[];
    var thehosescore =[];
    var Childgroup = [];
    var bondtypes = [];
    var symbol = [];
    var ranking = [];
    //console.log("prev:",previousGroup);

	  // Determine atoms surrounding each atom and build HOSE code list
    for (var i=0; i<currentsphere.length;i++) {
	  	var  hs,tmp_str;
	  	var atomid = currentsphere[i]; 

	  	if (!dejavu[atomid]) {
	  		isdejavu=false;
	  	}
	  	else {
	  		isdejavu = true;

	  	}

	  	dejavu[atomid]=true;
	    
	    symbol[i]= atomlistvect[atomid];
	  	bondtypes[i] = bondtype[atomid];

	  	var hose_code =HC(atomlistvect[atomid],bondtype[atomid],isdejavu);

	  	tmphose[i] = hose_code.tmp_str;
	  	thehosescore[i] = hose_code.hs;


	  	symbol[i]= atomlistvect[atomid];
	  	ranking[i]= i+1;
    	Childgroup.push(child[atomid]);

	  }

	  //spherecode.push(currentsphere,tmphose,thehosescore,previousGroup);
	  S.push(new Node({"bondtype":bondtypes,"atomid" :currentsphere,"symbols":symbol,"ranking":ranking,"childs": Childgroup, "hosescore": thehosescore, "hosecode":tmphose }));
	  
	  previousGroup =Childgroup;
	}
	return S;
}


function loopotherSpheres(Spheres,graph){
	var dejavu=[];
	var currentnodes = Spheres[0].value.atomid;
	for (var ni=0; ni< currentnodes.length; ni++) {
		dejavu[currentnodes[ni]]= true;

	}
	var initalscore  = Spheres[0].value.hosescore;
	var initalhose =  Spheres[0].value.hosecode;
		console.log("Sphere: 1 ------------");
		console.log(currentnodes,",score:",initalscore,"hosecode:",initalhose);

	for (var i=1; i<Spheres.length;i++) {
		console.log("Sphere:",i+1,"------------")
	 	var nextSphere = Spheres[i-1].value.childs;
	 	var currentchilds = Spheres[i].value.childs;
	 	var currentnodes = Spheres[i].value.atomid;
	 	var parentnodes = Spheres[i-1].value.atomid;


	 	var currentsymbol = Spheres[i].value.symbols;
	 	var currentbondtype = Spheres[i].value.bondtype;
	 	var parenthosescore = Spheres[i-1].value.hosescore;

	 	for (var j=0; j<nextSphere.length;j++)
	 		 {
	 		 	var childnode = nextSphere[j];

	 		 	for (var sj=0; sj<childnode.length;sj++){
	 		 		var myatomid =childnode[sj];
		 		    var previous = parentnodes[j];
	 		 		console.log("current atomid:",myatomid,"<=",previous);
		 		 	if (dejavu[myatomid]) {
		 		 		var isdejavu = true;
		 		 		var mybondtype = graph[previous][myatomid]; // look like this is not always the good value for dejavu cases!
		 		 		var mysymbol = '&';
		 		 		var hcode =HC(mysymbol,mybondtype,isdejavu);
		 		 		console.log("currentatomidindex:",myatomid,",bond:",mybondtype,";dejavu:",hcode);
		 		 	}
		 		 	else {
		 		 		dejavu[myatomid]=true;
		 		 		var isdejavu = false;
		 		 		var currentatomidindex = currentnodes.indexOf(myatomid);
		 		 		var mybondtype = currentbondtype[currentatomidindex];
		 		 		var mysymbol = currentsymbol[currentatomidindex];
		 		 		var hsval = HC(mysymbol,mybondtype,isdejavu);
		 		 		console.log("currentatomidindex:",currentatomidindex,",bond:",mybondtype,",symbol:",mysymbol,",HOSE:",hsval);
		 		 	}
	 		 	}
	 		 }
	  }
	  return 1;
}






function childs(sphere, Neiarray, parents) {
	child = [];
	var atomid = null;
	for (var j=0;j<sphere.length;j++) {
	    var currentsphere = sphere[j];
		  // Determine atoms surrounding each atom and build HOSE code list
	    for (var i=0; i<currentsphere.length;i++) {
		  	var atomid = currentsphere[i]; 
	    	child[atomid]=diff(Neiarray[atomid],parents[atomid]);	    
	    }
	}
	return child;
}

function ScoreandSort(spherecode) {
	console.log(spherecode);
	var linker = [];
	var cscore = [];
	var atoms = [];
	var dupindex = [];
	for(var i=0; i<spherecode.length/4; i++) {
			
		    for (var u=0; u<linker.length;u++){

		    	console.log("linker loop",linker[u]);
		    }
		    console.log("sphereid:",i+1);
		    // iteration on the next sphere to find duplicate and inject them to modify the score at the next level
		    cscore=spherecode[i*4+2];
		    atoms=spherecode[i*4];

			var dup = cscore.getDuplicates();
			var keys =Object.keys(dup);
			for (var j=0; j<keys.length;j++) {
				dupindex = dup[keys[j]];

				for (var k=0;k<dupindex.length;k++) {
					 linker[k] = atoms[dupindex[k]];
				}
			}
	}
	console.log("cummulative score:",cscore);
	return cscore;

}




Array.prototype.getDuplicates = function () {
    var duplicates = {};
    for (var i = 0; i < this.length; i++) {
        if(duplicates.hasOwnProperty(this[i])) {
            duplicates[this[i]].push(i);
        } else if (this.lastIndexOf(this[i]) !== i) {
            duplicates[this[i]] = [i];
        }
    }

    return duplicates;
};


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
var Neiarray = [];
var Bonarray = [];
// return the Atom & bond Neighbors
for (j=0;j<C;j++){
 	Neiarray.push(readvector(mol.getAtomNeighbors(j)));
 	Bonarray.push(readvector(mol.getBondNeighbors(j)));
}


var g = mol.getAdjacencyMatrix(true);
var graph =readmat(g);
//var adjlist = AdjList(graph);
// 
var atomlist = mol.getSymbols();
var atomlistvect = readvector(atomlist);

var rings = mol.findSSSR();
var ringsval= [];

for (i =0;i<rings.size();i++){  
	ringsval[i] = rings.get(i); 
}


console.log("--------------------------");


var  shortestpath = bfs(graph, 5, numSpheres);

var HSC = HosecodeSphere(shortestpath.spheres, atomlistvect, shortestpath.bondtype, ringsval, Neiarray,shortestpath.parents);

//// return the HSC structure of childs ,...
/*for (i=0; i<HSC.length;i++) {
console.log("S",i);
console.log("atomids:",HSC[i].value.atomid);
console.log("childids:",HSC[i].value.childs);
console.log("bondtypes:",HSC[i].value.bondtype);
console.log("symbol:",HSC[i].value.symbols);
console.log("rank:",HSC[i].value.ranking);
console.log("hose score:",HSC[i].value.hosescore);
console.log("hose code:",HSC[i].value.hosecode);
}*/


loopotherSpheres(HSC,graph);



mol.delete();

