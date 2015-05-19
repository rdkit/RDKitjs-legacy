(function(factory){
  if (typeof define === 'function' && define.amd) {
    define([], function(){
      return factory();
    });
  } else if (typeof exports === 'object') {
    module.exports = factory();
  } else {
    window.RDKit = factory();
  }
})(function(){

