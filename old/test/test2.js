var fs = require('fs');
var rdk = require('rdkit');

var parser = require('sdf-parser')


//var result = parser('/Users/guilllaumegodin/Downloads/chembl_18.sdf');
var timestamp = (new Date()).getTime();
var LineInputStream = require('line-input-stream')

var stream = LineInputStream(fs.createReadStream('chembl_18.sdf', {flags:'r'}));
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

function processLine(line) {
    
    //console.log(line);
    try {
    var mol= rdk.Molecule.MolBlockToMol(line);
        }
    
    catch (err) {console.log('err');return;}
    var smile = mol.smilewrite();
    //console.log(chemblid);
    //console.log:);
    
    mol.delete();
    
    
   smiles.push("a");

    if (smiles.length>10000) {
      var date1 = new Date(timestamp);
      var date2 = new Date();
      var timeDiff = Math.abs(date2.getTime() - date1.getTime());
      var diffSec = parseInt(timeDiff / (1000 ));
      var diffMin = parseInt(diffSec/60);
      diffSec = diffSec - (diffMin*60);
      var time = diffMin + "m" + diffSec + "s"; 
      console.log(" in "+time);
        //console.log(smiles);
        console.log('done');
        process.exit();
    }
}



/*
//var LineByLineReader = require('line-by-line'),
lr = new LineByLineReader('/Users/guilllaumegodin/Downloads/chembl_18.sdf', );

lr.on('error', function (err) {
 console.log(err);
});
var data ="";
var i=0;
var data="";
var tmp = [];

lr.on('line', function (line) {
lr.pause();
data = data + line +"\r\n";
 //  stop molecule identificator $$$$
if(line.indexOf("$$$$") !== -1){
  i++
  for(var j=0;j < tmp.length;j++){
   var dataToInsert="";
   if(tmp.length === j+1){
    dataToInsert = tmp[j] + data;
            //console.log("last tmp",dataToInsert);
      }
        else{
         dataToInsert = tmp[j] + "$$$$\r\n"; 
           //console.log("tmp",dataToInsert);
      }
    rdkitprocess(dataToInsert);

    }
    tmp = [];
    data = "";
}
      
      
else if(line.indexOf("CHEMBL") !== -1){
       tmp.push(data);
       data = "";
}



 if (i===100)  {
       lr.close();
       throw new Error('to many lines for me'); 
      }


    lr.resume();
   });

lr.on('end', function () {
    // All lines are read, file is closed now.

});



function rdkitprocess(myval){
    var rdk = require('rdkit');

            console.log(i);
            //fs.appendFileSync('message.txt', "i="+i+myval);
   try {
                var mol= rdk.Molecule.MolBlockToMol(myval);
                //var chemblid = mol.getProp('_MolFileInfo');
                var smile = mol.smilewrite();
                //console.log(chemblid);
                //console.log:);
                mol.delete();
            }
   catch (err) {
                console.error(err);
                fs.appendFileSync('error.txt', "i="+i+"\r\n"+JSON.stringify(myval)+"\r\n");
            
            }
        }*/