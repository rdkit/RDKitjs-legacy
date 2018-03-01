'use strict';

const Runner = require('../utils/Runner');

class Compiler extends Runner {
  compile() {
    console.log(this.emscriptenPath);
  }
}

module.exports = Compiler;
