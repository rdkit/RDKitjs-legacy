'use strict';

const childProcess = require('child_process');
const os = require('os');
const { join, relative } = require('path');

const request = require('request');
const fs = require('fs-extra');
const tar = require('tar');

const Runner = require('../utils/Runner');

const deps = require('./deps');

const isWindows = process.platform === 'win32';
const cpus = os.cpus().length;
const parallel = cpus > 1 ? cpus - 1 : 1;

class Installer extends Runner {
  constructor() {
    super();
    this.deps = {};
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
        console.log(`Downloading ${dep.name} version ${dep.version} ...`);
        await getDep(dep.url, depFolder);
        console.log(
          `${dep.name} extracted to ${relative(this.projectDir, depFolder)}`
        );
      }
    }
    console.log('All dependencies installed');
  }

  findMSBuild() {
    if (isWindows) {
      try {
        childProcess.execSync('MSBuild /version');
      } catch (e) {
        console.error('MSBuild.exe not found');
        console.error('You must install Visual Studio');
        console.error(
          'If Visual Studio is already installed, run the script in a developer command prompt'
        );
        return false;
      }
    }

    return true;
  }

  async compileRdkit() {
    console.log('Compiling RDKit');
    try {
      childProcess.execSync('cmake --version');
    } catch (e) {
      console.error('cmake not found');
      return false;
    }

    const rdkitBuildDir = join(this.deps.rdkit, 'build');
    await fs.ensureDir(rdkitBuildDir);

    const cmakeCommand = [
      'cmake ..',
      `-DCMAKE_TOOLCHAIN_FILE=${join(
        this.emscriptenPath,
        'cmake/Modules/Platform/Emscripten.cmake'
      )}`,
      `-DBoost_INCLUDE_DIR=${this.deps.boost}`,
      `-DEIGEN3_INCLUDE_DIR=${this.deps.eigen}`,
      '-DRDK_BUILD_PYTHON_WRAPPERS=OFF',
      '-DRDK_BUILD_CPP_TESTS=OFF',
      '-DRDK_BUILD_SLN_SUPPORT=OFF',
      '-DTHREADS_PTHREAD_ARG=OFF'
    ];

    childProcess.execSync(cmakeCommand.join(' '), {
      cwd: rdkitBuildDir,
      stdio: 'inherit'
    });

    if (isWindows) {
      childProcess.execSync(
        `MSBuild.exe /m:${parallel} /p:Configuration=Release INSTALL.vcxproj`,
        {
          cwd: rdkitBuildDir,
          stdio: 'inherit'
        }
      );
    } else {
      childProcess.execSync(`make -j${parallel}`, {
        cwd: rdkitBuildDir,
        stdio: 'inherit'
      });
    }

    return true;
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
