//Var to defined
var nbTodo = 2000000;
var waitOnError = 3000;//in ms
var fileToRead = 'chembl_18.sdf';//File to process
var nbretry = 1;



//###################################################################
//###################################################################
//###################################################################
var util = require('util');
var os = require('os') ;
var rdk = require('rdkit');
function monitor(name){
  console.log("#####################");
  console.log("Name:" + name);
  console.log(" memory",util.inspect(process.memoryUsage()));
  console.log(" os memory"," totalmem: "+os.totalmem(), " freemem: "+os.freemem());
  console.log("#####################");
}

monitor("start");

//###################################################################
//###################################################################
//###################################################################
//###################################################################


var fs = require('fs');
var LineByLineReader = require('line-by-line');
var lr = new LineByLineReader(fileToRead);
var MongoClient = require('mongodb').MongoClient;
var url = 'mongodb://gva6506.gva.europe.firmenich.com/:27017/rdkit';
var moleculesCollection;
var databaseMongo;

MongoClient.connect(url, function(err, db) {
  if(err){console.log(err);}
  databaseMongo=db;
  processFile();
  moleculesCollection = databaseMongo.collection('molecules');
});

function processFile(){

  //Static Var
  var i=0;
  var data="";
  var tmp = [];
  var nbsuccess = 0;
  var nberror = 0;
  var timestamp = (new Date()).getTime(); 

  try{
    fs.unlinkSync("message.txt");
  }catch(err){

  }
  try{
    fs.unlinkSync("error.txt");
  }catch(err){

  }

  lr.on('error', function (err) {
   console.log("lr error", err);
  });
  lr.on('end', function () {
    var date1 = new Date(timestamp);
    var date2 = new Date();
    var timeDiff = Math.abs(date2.getTime() - date1.getTime());
    var diffSec = parseInt(timeDiff / (1000 ));
    var diffMin = parseInt(diffSec/60);
    diffSec = diffSec - (diffMin*60);
    var time = diffMin + "m" + diffSec + "s"; 
    console.log("success: "+nbsuccess+" failed: "+nberror+" in "+time);
    monitor("end");
    
    //@todo
    //moleculesCollection.ensure_index('mfp.bits');
    //moleculesCollection.ensure_index('mfp.count');  
    databaseMongo.close();
  });

  lr.on('line', function (line) {
    if (i>=nbTodo) {
      lr.close();
      return;
    }

    line = line +"\r\n";  
    data = data + line;

    //  stop molecule identificator $$$$
    if(line.indexOf("$$$$") !== -1){
      i++;
      for(var j=0;j < tmp.length;j++){
        var dataToInsert="";
        if(tmp.length === j+1){
          dataToInsert = tmp[j] + data;
        }else{
          dataToInsert = tmp[j] + "$$$$\r\n"; 
        }
        var success = rdkitprocess(dataToInsert);
        if(success === true){
          // console.log("#" + i + " success");
          nbsuccess++;
        }else{
          // console.log("#" + i + " failed");
          nberror++;
        }
      }
      tmp = [];
      data = "";
    }else if(line.indexOf("CHEMBL") !== -1){
      tmp.push(data);
      data = "";
    } 

  });


  function rdkitprocess(myval){
    fs.appendFileSync('message.txt', "i="+i+myval);
    var attempt = 0;
    var success = false;
    while(attempt < nbretry && success !== true){
      success = callRdkit(myval);
      attempt++;
      // if(success !== true ){
      //   console.log("#" + i + " attempt failed: " + attempt + " with error:" + success);
      //   console.log("wait: " + waitOnError/1000 + "s")
      //   sleep(waitOnError);
      // }        
    }
    return success;
  }

  function savetomongo(mol){
    //@todo a reprendre
    var binarymol = mol.pickleMol();
    var mpf = mol.getMorganFingerprints(2,2048);
    mfp = AllChem.GetMorganFingerprint(rdmol, 2);
    
    //rdmol = Chem.Mol(molecule['rdmol']);
    mfp = list(AllChem.GetMorganFingerprintAsBitVect(rdmol, 2, nBits=2048).GetOnBits());
    moleculesCollection.insert({
      'smiles': Chem.MolToSmiles(rdmol, isomericSmiles=True),
      'chembl_id': rdmol.GetProp('chembl_id'),
      'rdmol': Binary(rdmol.ToBinary()),
      'mfp': {'bits': mfp, 'count': len(mfp)}
    });
  }

  function callRdkit(myval){
    
    try {
      var mol = rdk.Molecule.MolBlockToMol(myval);
      //savetomongo(mol);
      mol.delete();
      return true;
    }
    catch (err) {
      fs.appendFileSync('error.txt', "i="+i+"\r\n"+JSON.stringify(myval)+"\r\n");
      return err;
    }
  }
  function sleep(miliseconds) {
     var currentTime = new Date().getTime();
     while (currentTime + miliseconds >= new Date().getTime()) {
     }
  }
}