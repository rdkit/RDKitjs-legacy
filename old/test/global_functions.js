var util = require('util');
var os = require('os') ;
 
var currentTimestamp = (new Date()).getTime();
var totalmem = beginTotalmem = os.totalmem();
var freemem = beginFreemem = os.freemem();
 
exports.formatBytes = formatBytes;
 
exports.displayMemoryUsage = function displayMemoryUsage(){
    displayTotalMemory();
    displayFreeMemory();
}
exports.displayTotalMemory = function displayTotalMemory(){
    console.log("\ttotalmem: \t" + formatBytes(os.totalmem()) + "(" + formatBytes(os.totalmem()-totalmem) + ")");
    totalmem = os.totalmem();
}
exports.displayFreeMemory = function displayFreeMemory(message){
    message = message ? message : "";
    console.log(message + "\tfreemem: \t" + formatBytes(os.freemem()) + " (" + formatBytes(os.freemem()-freemem) + ")");
    freemem = os.freemem();
}
exports.getPassedTime = function getPassedTime(){
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
 
exports.monitor = function monitor(name){
    displayFreeMemory(name + "\t" + getPassedTime());
}
 
exports.replaceAll = function replaceAll (string, find, replace) {
                return string.split(find).join(replace);
}
 
exports.toTimeString = function toTimeString(timeDiff){
  var diffMil = parseInt(timeDiff);
  var diffSec = parseInt(timeDiff / (1000 ));
  var diffMin = parseInt(diffSec/60);
  diffSec = diffSec - (diffMin*60);
  diffMil = diffMil - ( ( (diffMin*60) + diffSec) * 1000 );
  var time = diffMin + "m" + diffSec + "s" + diffMil + "ms";
  return time;
}
 
exports.resumeTime = function resumeTime(message){
    console.log(message);
    console.log("\ttotalmem: \t" + formatBytes(os.totalmem()) + "(" + formatBytes(os.totalmem()-beginTotalmem) + ")");
    console.log("\tfreemem: \t" + formatBytes(os.freemem()) + " (" + formatBytes(os.freemem()-beginFreemem) + ")");
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
 
module.export = exports;