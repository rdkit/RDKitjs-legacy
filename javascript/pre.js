(function(factory){
  if (typeof define === 'function' && define.amd) {
    define(function(){
      return factory();
    });
  } else if (typeof module === "object" && module.exports) {
    module.exports = factory();
  } else {
    window.RDKit = factory();
  }
})(function(){

