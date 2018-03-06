'use strict';

// eslint-disable-next-line import/no-unassigned-import
require('make-promises-safe');

const Installer = require('./Installer');

(async () => {
  const installer = new Installer();
  await checkExit(installer.init());
  await checkExit(installer.installDeps());
  await checkExit(installer.findEmscripten());
  await checkExit(installer.findMSBuild());
  await checkExit(installer.compileRdkit());
})();

async function checkExit(promise) {
  const result = await promise;
  if (result === false) {
    console.error('An error occured');
    process.exit(1);
  }
}
