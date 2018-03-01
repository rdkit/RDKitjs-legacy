'use strict';

const { join, resolve, relative } = require('path');

require('make-promises-safe');
const request = require('request');
const fs = require('fs-extra');
const tar = require('tar');

const deps = require('./deps');

(async () => {
  const projectDir = join(__dirname, '../..');
  const depsDir = join(projectDir, 'deps');

  await fs.ensureDir(depsDir);
  console.log('Checking dependencies');
  for (const dep of deps) {
    const depFolder = join(depsDir, dep.name, dep.version);
    const fileCheck = join(depFolder, dep.fileCheck);
    const depExists = await fs.exists(fileCheck);
    if (!depExists) {
      console.log(
        'Downloading ' + dep.name + ' version ' + dep.version + '...'
      );
      await getDep(dep.url, depFolder);
      console.log(
        `${dep.name} extracted to ${relative(projectDir, depFolder)}`
      );
    }
  }
})();

async function getDep(url, folder) {
  await fs.ensureDir(folder);
  const fetchStream = request(url);
  const tarStream = fetchStream.pipe(
    tar.x({
      strip: 1,
      C: folder
    })
  );
  await new Promise((res, rej) => {
    fetchStream.on('error', rej);
    tarStream.on('error', rej);
    tarStream.on('end', res);
  });
}
