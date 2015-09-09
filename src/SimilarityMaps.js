var num = require('numeric');
var RDKit = require('rdkit');


function Uni_Nor(X,Y,sigma,mux,muy)
{
	var Xmu=num.sub(X,num.rep([51,51],mux));
	var Ymu=num.sub(Y,num.rep([51,51],muy));
    var deno = num.rep([51,51],sigma*sigma);
 	var le = num.div(num.mul(Xmu,Xmu),deno);
 	var re = num.div(num.mul(Ymu,Ymu),deno); 
	var z = num.add(le,re);
	var denom = 6.29*sigma*sigma;

	return num.div(num.exp(num.neg(num.div(z,num.rep([51,51],2)))),num.rep([51,51],denom));
}

function calcAtomGaussians(sigma,step,weights,px,py)
{
    // default step 0.02 => 51 points
	var x = num.linspace(0,1,51);
	// meshgrid square grid X=Y'; 
	var X = num.rep([51,],x);
	var Y = num.transpose(X);
	//var Z = num.mul(Uni_Nor(X,Y,sigma,px[0],py[0]),num.rep([51,51],weights[0]));
	//console.log(num.sum(Z));
	//console.log(num.dim(Z));
	
	var Z=num.rep([51,51],0);
	console.log("Z is up");
	console.log(px.size());
	
	for (j=0;j<px.size();j=j+1){   
		Z= num.addeq(num.mul(Uni_Nor(X,Y,sigma,px[j],py[j]),num.rep([51,51],weights[0])));
	}
	return Z
}



var m=RDKit.Molecule.fromSmiles('Cn1cc(NC=O)cc1C(=O)Nc1cc(C(=O)Nc2cc(C(=O)NCCC(N)=[NH2+])n(C)c2)n(C)c1');
var w = m.getTPSAAtomContribs();
var ap = m.getAtomsPos2D();

var px= [];
var py= [];
var weigths = [];
for  (i=0;i<ap.size()-1;i=i+2)
{	
	weigths[i/2]=w.get(i/2);
	px[i]=ap.get(i)+1;
	py[i]=ap.get(i+1)+1;
}

console.log(px);
console.log(py);
/*
var xmin = 1000;
var xmax =-1000;
var ymin = 1000;
var ymax =-1000;
for  (i=0;i<x.size();i=i+1)
{	
 
 if (x[i]> xmax)
 	xmax = x[i];
 if (x[i]< xmin)
 	xmin = x[i];

 if (y[i]> ymax)
 	ymax = y[i];
 if (y[i]< ymin)
 	ymin = y[i];
}
*/


var t = ap.size();

// need to find why there is an issue to retreive the ap.get(t-1) value ?
var simmap=calcAtomGaussians(0.45,0.02,weigths,px,py);











