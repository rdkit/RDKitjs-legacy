'use strict';

const bindings = require('../dist/rdkit.js');

const Mol = require('./Mol');

module.exports = {
  Mol,
  load: () => bindings.load()
};
