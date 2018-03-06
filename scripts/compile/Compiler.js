'use strict';

const childProcess = require('child_process');
const { join } = require('path');

const fs = require('fs-extra');

const Runner = require('../utils/Runner');

const debug = process.argv.includes('--debug');

class Compiler extends Runner {
  async compile() {
    const emscriptenRoot = join(this.emscriptenPath, '../..');
    const emscriptenScript = join(emscriptenRoot, 'emsdk_env.sh');
    const out = join(this.projectDir, 'out');
    await fs.ensureDir(out);

    const jsConfigFlags = [
      'WASM=1',
      'NODEJS_CATCH_EXIT=0',
      'MODULARIZE=1',
      "EXPORT_NAME='\"'rdk'\"'"
    ];

    if (debug) {
      jsConfigFlags.push('DISABLE_EXCEPTION_CATCHING=0');
    }

    const jsConfig = jsConfigFlags.map((s) => `-s ${s}`).join(' ');

    const rdkitPath = this.deps.rdkit.path;

    const includes = [join(rdkitPath, 'Code'), this.deps.boost.path]
      .map((s) => `-I${s}`)
      .join(' ');

    const libFiles = [
      'libRDKitAlignment.so',
      'libRDKitDistGeometry.so',
      'libRDKitDistGeomHelpers.so',
      'libRDKitEigenSolvers.so',
      'libRDKitFileParsers.so',
      'libRDKitForceField.so',
      'libRDKitForceFieldHelpers.so',
      'libRDKitGraphMol.so',
      'libRDKitMolAlign.so',
      'libRDKitMolTransforms.so',
      'libRDKitRDGeneral.so',
      'libRDKitRDGeometryLib.so',
      'libRDKitSmilesParse.so',
      'libRDKitSubstructMatch.so'
    ]
      .map((s) => join(rdkitPath, 'build/lib', s))
      .join(' ');

    const optimize = debug ? '-O1' : '-Os';

    const emcc = [
      'emcc --bind',
      optimize,
      jsConfig,
      '-o out/rdkit.js',
      includes,
      libFiles,
      'src/rdkit.cc'
    ];

    childProcess.execSync(
      `source '${emscriptenScript}'  && ${emcc.join(' ')}`,
      { shell: '/bin/bash', stdio: 'inherit' }
    );
  }
}

module.exports = Compiler;
