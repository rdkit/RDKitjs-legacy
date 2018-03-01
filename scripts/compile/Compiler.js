'use strict';

const childProcess = require('child_process');
const { join } = require('path');

const fs = require('fs-extra');

const Runner = require('../utils/Runner');

class Compiler extends Runner {
  async compile() {
    const emscriptenRoot = join(this.emscriptenPath, '../..');
    const emscriptenScript = join(emscriptenRoot, 'emsdk_env.sh');
    const out = join(this.projectDir, 'out');
    await fs.ensureDir(out);

    const jsConfig = [
      'WASM=1',
      'NODEJS_CATCH_EXIT=0',
      'MODULARIZE=1',
      "EXPORT_NAME='\"'RDK'\"'"
    ]
      .map((s) => `-s ${s}`)
      .join(' ');

    const rdkitPath = this.deps.rdkit.path;

    const includes = [join(rdkitPath, 'Code'), this.deps.boost.path]
      .map((s) => `-I${s}`)
      .join(' ');

    const libFiles = [
      'libRDKitSmilesParse.so',
      'libRDKitGraphMol.so',
      'libRDKitRDGeneral.so',
      'libRDKitFileParsers.so',
      'libRDKitDistGeometry.so',
      'libRDKitForceFieldHelpers.so',
      'libRDKitRDGeometryLib.so',
      'libRDKitDistGeomHelpers.so',
      'libRDKitForceField.so',
      'libRDKitEigenSolvers.so',
      'libRDKitAlignment.so',
      'libRDKitSubstructMatch.so'
    ]
      .map((s) => join(rdkitPath, 'build/lib', s))
      .join(' ');

    const emcc = [
      'emcc --bind -Os',
      jsConfig,
      '-o out/rdkit.js',
      includes,
      libFiles,
      'rdkit.cc'
    ];

    childProcess.execSync(`source '${emscriptenScript}'  && ${emcc.join(' ')}`);
  }
}

module.exports = Compiler;
