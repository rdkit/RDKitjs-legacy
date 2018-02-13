(function(factory) {
  function getCurrentPathBrowser() {
    sc = document.getElementsByTagName('script');

    for (idx = 0; idx < sc.length; idx++) {
      s = sc.item(idx);

      if (s.src && s.src.match(/rdkit\.js$/)) {
        return s.src;
      }
    }
  }
  function getModuleForBrowser() {
    var resolve;
    var promise = new Promise((res, rej) => {
      resolve = res;
    });
    var path = getCurrentPathBrowser();
    var Module = {
      locateFile(url) {
        return new URL(url, path).href;
      },
      onRuntimeInitialized() {
        resolve();
      },
      load() {
        return promise;
      }
    };
    return Module;
  }
  function getModuleForNode() {
    var resolve;
    var promise = new Promise((res, rej) => {
      resolve = res;
    });
    var Module = {
      locateFile(url) {
        return __dirname + '/' + url;
      },
      onRuntimeInitialized() {
        resolve();
      },
      load() {
        return promise;
      }
    };
    return Module;
  }
  if (typeof define === 'function' && define.amd) {
    define(function() {
      return factory(getModuleForBrowser());
    });
  } else if (typeof module === 'object' && module.exports) {
    module.exports = factory(getModuleForNode());
  } else {
    window.RDKit = factory(getModuleForBrowser());
  }
})(function(Module, module, exports, define) {
