var mongourl = 'mongodb://localhost:27017/rdkit';
//var collectionName = "molecules1";
var collectionName = "molecules2";
var start = 1;//start = 1
var end = 1000000;//include it
var stepToShow = 1000;
var fileToRead = '/Users/mbp/Document/db/mol1.sdf';//File to process
//var fileToRead = '/Users/mbp/Downloads/chembl20/chembl_20.sdf';//File to process

var errorFile = 'errors.txt';
var encoding = 'ascii';
var dolineInputStream = true;
 
//###################################################################
//###################################################################
//###################################################################
var fs = require('fs');
var rdk = require('rdkit');
var util = require('util');
var os = require('os');
var memwatch = require('memwatch-next');
var parse = require('sdf-parser');
 
var hd;
memwatch.on('leak', function(info) {
console.error(info);
if (!hd) {
   hd = new memwatch.HeapDiff();
} else {
   var diff = hd.end();
   console.error(util.inspect(diff, true, null));
   hd = null;
}
});
 
try{
    fs.unlinkSync(errorFile);
}catch(err){
 
}
 
function formatBytes(bytes) {
   var sign = bytes > 0 ? 1 : bytes == 0 ? 0 : -1;
   bytes = Math.abs(bytes);
   var decimals = 2;
   if(bytes == 0) return '0 Byte';
   var k = 1000;
   var dm = decimals + 1 || 3;
   var sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB'];
   var i = Math.floor(Math.log(bytes) / Math.log(k));
   return (bytes / Math.pow(k, i)).toPrecision(dm)*sign + ' ' + sizes[i];
}
 
function displayMemoryUsage(){
    displayTotalMemory();
    displayFreeMemory();
}
function displayTotalMemory(){
    console.log("\ttotalmem: \t" + formatBytes(os.totalmem()) + "(" + formatBytes(os.totalmem()-totalmem) + ")");
    totalmem = os.totalmem();
}
function displayFreeMemory(message){
    message = message ? message : "";
    console.log(message + "\tfreemem: \t" + formatBytes(os.freemem()) + " (" + formatBytes(os.freemem()-freemem) + ")");
    freemem = os.freemem();
}
function getPassedTime(){
    var date1 = new Date(currentTimestamp);
    var date2 = new Date();
    var timeDiff = Math.abs(date2.getTime() - date1.getTime());
    var diffSec = parseInt(timeDiff / (1000 ));
    var diffMin = parseInt(diffSec/60);
    diffSec = diffSec - (diffMin*60);
    var time = diffMin + "m" + diffSec + "s";
    currentTimestamp = date2.getTime();
    return "+"+time;
}
 
function monitor(name){
    displayFreeMemory(name + "\t" + getPassedTime());
}
 
function resumeTime(message){
    console.log(message);
    console.log("\ttotalmem: \t" + formatBytes(os.totalmem()) + "(" + formatBytes(os.totalmem()-beginTotalmem) + ")");
    console.log("\tfreemem: \t" + formatBytes(os.freemem()) + " (" + formatBytes(os.freemem()-beginFreemem) + ")");
}
 
var currentTimestamp = (new Date()).getTime();
var totalmem = beginTotalmem = os.totalmem();
var freemem = beginFreemem = os.freemem();
displayMemoryUsage();
 
 
var Binary = require('mongodb').Binary;
var MongoClient = require('mongodb').MongoClient;
 
var moleculesCollection;
var databaseMongo;
 
MongoClient.connect(mongourl,
    function(err, db) {
        if(err){fs.appendFileSync(errorFile, "\r\n"+JSON.stringify(err)+"\r\n");}
        databaseMongo=db;
        moleculesCollection = databaseMongo.collection(collectionName);
        //moleculesCollection.remove();//remove all

        processFile();
    }
);
 
function processFile(){
    
    var timestamp = (new Date()).getTime();
   
    if(dolineInputStream === true){
        var LineInputStream = require('line-input-stream');
        var stream = LineInputStream(fs.createReadStream(fileToRead, {flags:'r'}));
        stream.setEncoding(encoding);
        stream.setDelimiter('$$$$\n');
        stream.on("error", function(err) {
            fs.appendFileSync(errorFile, "\r\n"+JSON.stringify(err)+"\r\n");
        });
        stream.on('line', processLine);
        stream.on("end", function() {processEnd();});
        stream.on("close", function() {processEnd();});
    }else{
        var LineByLineReader = require('line-by-line');
        var lr = new LineByLineReader(fileToRead, {encoding: encoding});
        lr.on('error', function (err) {
            fs.appendFileSync(errorFile, "\r\n"+JSON.stringify(err)+"\r\n");
        });
        lr.on('end', function () {processEnd();});
        var data="";
        lr.on('line', function (line) {
            data = data + line +"\r\n";
            if(line.indexOf("$$$$") !== -1){
                processLine(data);
                data = "";
           }
        });
    }
 
 
    var i=1;
    var stepper=stepToShow;
    var nbsuccess = 0;
    var nberror = 0;
    function processLine(line) {
       
        if(i >= stepper){
            monitor("#step "+i);
            stepper+=stepToShow;
        }
 
        if(i >= start){//define point for start
          // try {

/*                var Lines = line.split('\n');
                goodlines = Lines[0]+"\n";
                cscart = false;
                for (ij=1;ij<Lines.length;ij++)
                {   
                    if (Lines[ij].indexOf("CsCart") !==-1) cscart = true;
                    if (ij!=1 & cscart) goodlines+=Lines[ij]+"\n";
                }
*/
                var result = parse(line);
                line = result.molecules[0].molfile.value;
                console.log(line);

                var mol= rdk.Molecule.MolBlockToMol(line);
                var smile = mol.toSmiles();
                var binarymol = mol.MolToBinary();

                // this is the folded example (AllChem.GetHashedMorganFingerprint(mol, radius=self.radius, nBits = self.length))
                // the unfolded (faster method is based on : AllChem.GetMorganFingerprint(mol, radius=self.radius) 
                // the bits are sorted using sorted(GetNonZeroElement().keys())
                var mpf2 = mol.getMorganUFPkeys(2);   // unfolded keys storage using getMorganUFPkeys   // folded =  mol.getMorganFP_GetOnBits(2,2048); 
                var mpf = [];
                var mpflen = mpf2.size();
                for(var j=0;j < mpflen;j++){
                    mpf.push(mpf2.get(j));
                }             

                moleculesCollection.insert(
                    {
                       'smiles': smile,
                       'chembl_id': i,
                       'rdmol': Binary(binarymol),
                       'mfp': {'bits': mpf, 'count': mpflen}
                    }
                );
                mol.delete();
                nbsuccess++;         
         /* }catch (err) {
                fs.appendFileSync(errorFile, "i="+i+"\r\n"+JSON.stringify(err)+"\r\n");
                nberror++;
                return;
            }
*/
        }
 
        if (i >= end) {
            processEnd();
        }
        i++;
    }
 
    function processEnd(){
        var date1 = new Date(timestamp);
        var date2 = new Date();
        var timeDiff = Math.abs(date2.getTime() - date1.getTime());
        var diffSec = parseInt(timeDiff / (1000 ));
        var diffMin = parseInt(diffSec/60);
        diffSec = diffSec - (diffMin*60);
        var time = diffMin + "m" + diffSec + "s";
        resumeTime("###i=" + i + " success: "+(nbsuccess-nberror)+" failed: "+nberror+" in "+time);
        databaseMongo.close();
        process.exit();
    }
}