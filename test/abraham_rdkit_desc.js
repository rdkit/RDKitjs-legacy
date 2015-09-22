var fs = require('fs');
var RDKit = require('rdkit');
var now = require("performance-now");
var exec = require('child_process').exec;


var loadRDKit = RDKit.hello;
console.log('Hello');
var filename = "/Users/mbp/Documents/abraham.txt";
var buf =   fs.readFileSync(filename, "utf8");
var start = now();    
var lines = buf.split("\n");
var res = [];
for (var j = 0; j < lines.length; ++j) {
        var smi = lines[j];
        // need to add try ... catch check of kekulize process!!!
        try {
		var mol = RDKit.Molecule.fromSmiles(smi);	
        var c = RDKit.calc_all_desc(mol);
        var g = JSON.stringify(c);
        //console.log(g);
        res.push(g);
        //console.log(res.length);
        //console.log(ConvertToCSV(res));
      	mol.delete();
      }
      catch (err) {
      	// caution asynchr code!!!
      	console.log(j+1,smi);
        var cmd = 'obabel -:"'+ smi + '" -ocan';
      	exec(cmd, function(error, stdout, stderr) {
      			console.log("res:",stdout);
      			var mol = RDKit.Molecule.fromSmiles(stdout);	
        		var c = RDKit.calc_all_desc(mol);
        		var g = JSON.stringify(c);

      			mol.delete();
      			//console.log(g);
		});
      	
      }
  
    }
var end = now();
console.log((end-start).toFixed(3)) 
try {
fs.writeFileSync("/Users/mbp/Documents/rdkitjsdesc.json",res,'utf8');
}
catch (err){
	console.log(err);
}