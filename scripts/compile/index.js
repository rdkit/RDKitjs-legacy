'use strict';

require('make-promises-safe');

const Compiler = require('./Compiler');

(async () => {
  const compiler = new Compiler();
  await checkExit(compiler.init());
  await checkExit(compiler.findEmscripten());
  await checkExit(compiler.compile());
})();

async function checkExit(promise) {
  const result = await promise;
  if (result === false) {
    console.error('An error occured');
    process.exit(1);
  }
}
