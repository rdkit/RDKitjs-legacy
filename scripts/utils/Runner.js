'use strict';

const os = require('os');
const { join, normalize } = require('path');

const fs = require('fs-extra');

const deps = require('./deps');

class Runner {
  constructor() {
    this.projectDir = join(__dirname, '../..');
    this.isWindows = os.platform === 'win32';
    this.deps = {};
    for (const dep of deps) {
      this.deps[dep.name] = dep;
      dep.path = join(this.projectDir, 'deps', dep.name, dep.version);
    }
  }

  init() {
    if (this.isWindows) {
      console.error('Building on Windows in not supported yet');
      return false;
    }
    return true;
  }

  async findEmscripten() {
    console.log('Searching for emscripten');
    const emscriptenConfigFile = join(os.homedir(), '.emscripten');
    let emscriptenConfig;
    try {
      emscriptenConfig = await fs.readFile(emscriptenConfigFile, 'utf8');
    } catch (e) {
      console.error(
        `Could not find emscripten config file in ${emscriptenConfigFile}`
      );
      console.error(
        'Follow the installation instructions at https://kripken.github.io/emscripten-site/docs/getting_started/downloads.html#installation-instructions'
      );
      return false;
    }
    const split = emscriptenConfig.split(/[\r\n]+/);
    const root = split.find((line) => line.includes('EMSCRIPTEN_ROOT'));
    const end = root.substring(17);
    const path = normalize(end.substring(0, end.length - 1));
    console.log(`Found at ${path}`);
    this.emscriptenPath = path;
    return true;
  }
}

module.exports = Runner;
