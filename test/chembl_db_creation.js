
var fs = require('fs');
var rdk = require('rdkit');
//###################################################################
//###################################################################
//###################################################################
var util = require('util');
var os = require('os') ;
function monitor(name){
    console.log("#####################");
    console.log("Name:" + name);
    console.log(" memory",util.inspect(process.memoryUsage()));
    console.log(" os memory"," totalmem: "+os.totalmem(), " freemem: "+os.freemem());
    console.log("#####################");
}

monitor("start");


var Binary = require('mongodb').Binary;
var MongoClient = require('mongodb').MongoClient;
var url = 'mongodb://127.0.0.1/:27017/rdkit';
var moleculesCollection;
var databaseMongo;

MongoClient.connect(url, function(err, db) {
                    if(err){console.log(err);}
                    databaseMongo=db;
                    processFile();
                    moleculesCollection = databaseMongo.collection('molecules');
                    });


function processFile(){
    
var timestamp = (new Date()).getTime();

var parser = require('sdf-parser')


//var result = parser('/Users/guilllaumegodin/Downloads/chembl_18.sdf');

var LineInputStream = require('line-input-stream')

var stream = LineInputStream(fs.createReadStream('/Users/guilllaumegodin/Downloads/chembl_18.sdf', {flags:'r'}));
stream.setEncoding('ascii');
stream.setDelimiter('$$$$\n');


stream.on("error", function(err) {
          console.log(err);
          });

stream.on('line', processLine);

stream.on("end", function() {
          // No more data, all line events emitted before this event
          });

stream.on("close", function() {
          // Same as ReadableStream's close event
          });

var smiles = [];
    var i=0;
function processLine(line) {
    ++i;
    //console.log(line);
    try {
        var mol= rdk.Molecule.MolBlockToMol(line);
        
    
    var smile = mol.smilewrite();
    var binarymol = mol.MolToBinary();
    var mpf2 = mol.getMorganFP_GetOnBits(2,2048);
    
        var mpf = [];
        var mpflen = mpf2.size();
        for(var j=0;j < mpflen;j++){
            mpf.push(mpf2.get(j));
        }
        

        
        
        
        
    //console.log(mpf2.size());
    
    moleculesCollection.insert({
                               'smiles': smile,
                               'chembl_id': i,
                               'rdmol': Binary(binarymol),
                               'mfp': {'bits': mpf, 'count': mpflen}
                               });
    
    
    mol.delete();
    
        
    }catch (err) {
        console.log('err');
        return;
    }

    
    smiles.push(smile)
    if (smiles.length>100000) {
        //console.log(smiles);
        //console.log('done');
        var date1 = new Date(timestamp);
        var date2 = new Date();
        var timeDiff = Math.abs(date2.getTime() - date1.getTime());
        var diffSec = parseInt(timeDiff / (1000 ));
        var diffMin = parseInt(diffSec/60);
        diffSec = diffSec - (diffMin*60);
        var time = diffMin + "m" + diffSec + "s";
        console.log("success in: "+time);
       
        monitor("end");
        
        databaseMongo.close();

        process.exit();
    }
}
}
