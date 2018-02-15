
var fs = require('fs');
var rdk = require('rdkit');
var performance = require('performance-now')



var filename = "abraham.txt";
var buf =   fs.readFilesync(filename, "utf8");
    
var t1= performance.now(); 

var lines = buf.split("\n");
for (var j = 0; j < lines.length; ++j) {
        var smi = lines[j];
		var mol = RDKit.Molecule.smilesToMol(smi);	
        var c = RDKit.calc_all_desc(mol);
        var g = JSON.stringify(c);
        console.log(j);
        mol.delete();
    }
var t2= performance.now();
console.log(t2-t1); 
