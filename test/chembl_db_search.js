var fs = require('fs');
var async = require('async');
 
var _ = require('underscore');
var gf = require('./global_functions');
var Binary = require('mongodb').Binary;
var MongoClient = require('mongodb').MongoClient;
var moleculesCollection;
var mfp_countsCollection;
var databaseMongo;
//###################################################################
//###################################################################
//###################################################################
 
var mongourl = 'mongodb://localhost:27017/rdkit';
var collectionName = "molecules";
var timestamp = (new Date()).getTime();
var limitRefresh = 200000000;
var limitCompared = 10;
var thresholdArray = [0.9,0.7,0.8,0.9];
 
// add index into db from mongo cmd
// db.molecules.createIndex({'mfp.bits':1})
// db.molecules.createIndex({'mfp.count':1})
MongoClient.connect(mongourl,
  function(err, db) {
    if(err){console.log(err);}
    databaseMongo=db;
    moleculesCollection = databaseMongo.collection(collectionName);
    mfp_countsCollection = databaseMongo.collection("mfp_counts");
  
    benchmarkFPfrequency(limitCompared, function(err,res){
        if(err){console.log(err,res);}else{
          for (var i = res.length - 1; i >= 0; i--) {
            console.log(res[i]);
          };
        }
        endProcess();
    });
    // countDo();
  }
);
 
 
function countDo(){//3min26s
  mfp_countsCollection.remove();//remove all
  var counts = {};
 
  moleculesCollection.find({}, {limit:limitRefresh, fields:{"mfp":1}}).toArray(
      function(err, molecule) {
        for (var imol = molecule.length - 1; imol >= 0; imol--) {
         
          for (var ibit = molecule[imol].mfp.bits.length - 1; ibit >= 0; ibit--) {
            var bit =  molecule[imol].mfp.bits[ibit];
            counts[bit] = counts[bit] ? counts[bit] + 1 : 1;
          };
        };
        for(var count in counts){
          mfp_countsCollection.insert({'_id': count, 'count': counts[count]});
        }
        endProcess();
      }
    );
}
 
 
 
function benchmarkFPfrequency(limitBatch, callbackEnd){
 
  moleculesCollection.find({}, {limit:limitBatch, fields:{smiles:1}}).toArray(
    function(err, molecules) {
 
      async.map(thresholdArray, //for each threshold
        function(threshold, callbackThreshold){
 
          async.map(molecules,  //for each molecules in the database
            function(mol, callbackmol){
              FPfrequency(mol.smiles, threshold, callbackmol);
            },
            function(err, resultsMol){
                if(err){callbackThreshold(err);return;}
                var sum = 0;
                var min = 999999999;
                var max = 0;
                for( var i = 0; i < resultsMol.length; i++ ){
                    sum += resultsMol[i].timeDiff ;
                    if(min > resultsMol[i].timeDiff){
                      min = resultsMol[i].timeDiff;
                    }
                    if(max < resultsMol[i].timeDiff){
                      max = resultsMol[i].timeDiff;
                    }
                }
                var avg = sum/resultsMol.length;
                callbackThreshold(null, {max:max,min:min,limitBatch:limitBatch,threshold:threshold,sum:sum,avg:avg});
            }
          );
 
        },
        function(err, results){
          if(err){callbackEnd(err);return;}
          var res = [];
          for (var i = results.length - 1; i >= 0; i--) {
           
            res.push(
              'limitBatch: ' + results[i].limitBatch +
              '\tthreshold: ' + results[i].threshold +
              '\tmin: ' + gf.toTimeString(results[i].min) +
              '\tmax: ' + gf.toTimeString(results[i].max) +
              '\tsum: ' + gf.toTimeString(results[i].sum) +
              '\tavg: ' + gf.toTimeString(results[i].avg)
              );
          };
          callbackEnd(null,res);
        }
      );
    }
  );
 
 
}
 
 
function FPfrequency(smile, threshold, callback){
  try{
    smile = gf.replaceAll(smile, ' 0\n', '');
    var rdkit = require('rdkit');
    var qmol = rdkit.Molecule.fromSmiles(smile);
    var qfp2 = qmol.getMorganFP_GetOnBits(2,2048);
    qmol.delete();
    var qfp = [];
    var mpflen = qfp2.size();
    for(var j=0;j < mpflen;j++){
        qfp.push(qfp2.get(j).toString());
    }
    var qn = mpflen;
    var qmin = parseInt(qn * threshold)
    var qmax = parseInt(qn / threshold)
    var limitmfp = qn - qmin + 1;
    var qtimestamp = (new Date()).getTime(); 
 
    mfp_countsCollection.find({_id: {$in: qfp}}, {sort:{'count':1},limit:limitmfp}).toArray(
      function(err, docs) {
        if(err){callback(err);return;}
        var reqbits = _.pluck(docs,'_id');
        reqbits = reqbits.map(function(el){return parseInt(el);})
        qfp = qfp.map(function(el){return parseInt(el);});
        var aggregate = [
          {'$match': {'mfp.count': {'$gte': qmin, '$lte': qmax}, 'mfp.bits': {'$in': reqbits}}},
          {'$project': {
              'tanimoto': {'$let': {
                  'vars': {'common': {'$size': {'$setIntersection': ['$mfp.bits', qfp]}}},
                  'in': {'$divide': ['$$common', {'$subtract': [{'$add': [qn, '$mfp.count']}, '$$common']}]}
              }},
              'smiles': 1,
             'chembl_id': 1
          }},
          {'$match':  {'tanimoto': {'$gte': threshold}}}
        ];
        var response = moleculesCollection.aggregate(aggregate).toArray(
          function(err,docs){
            if(err){callback(err);return;}
            var date1 = new Date(qtimestamp);
            var date2 = new Date();
            var timeDiff = Math.abs(date2.getTime() - date1.getTime());
            callback(null, {results:docs, timeDiff:timeDiff});
          }
        );
      }
    );
 
  }catch(err){
    callback(err);
  }
}
 
function endProcess(){
 
  var date1 = new Date(timestamp);
  var date2 = new Date();
  var timeDiff = Math.abs(date2.getTime() - date1.getTime());
  var time = gf.toTimeString(timeDiff);
  gf.resumeTime("### end in "+time);
  databaseMongo.close();
  process.exit();
}