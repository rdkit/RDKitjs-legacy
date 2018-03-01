'use strict';

const child_process = require('child_process');
const os = require('os');
const { join, normalize, relative } = require('path');

const request = require('request');
const fs = require('fs-extra');
const tar = require('tar');

const deps = require('./deps');

const isWindows = process.platform === 'win32';

class Installer {
  constructor() {
    this.projectDir = join(__dirname, '../..');
    this.deps = {};
  }

  init() {
    if (isWindows) {
      console.error('Building on Windows in not supported yet');
      return false;
    }
  }

  async installDeps() {
    const depsDir = join(this.projectDir, 'deps');

    await fs.ensureDir(depsDir);
    console.log('Checking dependencies');
    for (const dep of deps) {
      const depFolder = join(depsDir, dep.name, dep.version);
      this.deps[dep.name] = depFolder;
      const fileCheck = join(depFolder, dep.fileCheck);
      const depExists = await fs.exists(fileCheck);
      if (!depExists) {
        console.log(
          'Downloading ' + dep.name + ' version ' + dep.version + '...'
        );
        await getDep(dep.url, depFolder);
        console.log(
          `${dep.name} extracted to ${relative(this.projectDir, depFolder)}`
        );
      }
    }
    console.log('All dependencies installed');
  }

  async findEmscripten() {
    console.log('Searching for Emscripten');
    const emscriptenConfigFile = join(os.homedir(), '.emscripten');
    const emscriptenConfig = await fs.readFile(emscriptenConfigFile, 'utf8');
    const split = emscriptenConfig.split(/[\r\n]+/);
    const root = split.find((line) => line.includes('EMSCRIPTEN_ROOT'));
    const end = root.substring(17);
    const path = end.substring(0, end.length - 1);
    this.emscriptenPath = normalize(path);
    console.log('Found at ' + path);
  }

  async findMSBuild() {
    if (isWindows) {
      try {
        child_process.execSync('MSBuild /version');
      } catch (e) {
        console.error('MSBuild.exe not found');
        console.error('You must install Visual Studio');
        console.error(
          'If Visual Studio is already installed, run the script in a developer command prompt'
        );
        return false;
      }
    }
  }

  async compileRdkit() {
    console.log('Compiling RDKit');
    try {
      child_process.execSync('cmake --version');
    } catch (e) {
      console.error('cmake not found');
      return false;
    }

    const rdkitBuildDir = join(this.deps.rdkit, 'build');
    await fs.ensureDir(rdkitBuildDir);

    const cmakeCommand = [
      'cmake ..',
      '-DCMAKE_TOOLCHAIN_FILE=' +
        join(this.emscriptenPath, 'cmake/Modules/Platform/Emscripten.cmake'),
      '-DBoost_INCLUDE_DIR=' + this.deps.boost,
      '-DEIGEN3_INCLUDE_DIR=' + this.deps.eigen,
      '-DRDK_BUILD_PYTHON_WRAPPERS=OFF',
      '-DRDK_BUILD_CPP_TESTS=OFF',
      '-DRDK_BUILD_SLN_SUPPORT=OFF',
      '-DTHREADS_PTHREAD_ARG=OFF'
    ];

    child_process.execSync(cmakeCommand.join(' '), {
      cwd: rdkitBuildDir,
      stdio: 'inherit'
    });

    if (isWindows) {
      child_process.execSync(
        'MSBuild.exe /m:4 /p:Configuration=Release INSTALL.vcxproj',
        {
          cwd: rdkitBuildDir,
          stdio: 'inherit'
        }
      );
    } else {
      child_process.execSync('make', { cwd: rdkitBuildDir, stdio: 'inherit' });
    }
  }
}

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

module.exports = Installer;
