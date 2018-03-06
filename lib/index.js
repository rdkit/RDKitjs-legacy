'use strict';

const bindings = require('../dist/rdkit.js');

const Mol = require('./Mol');

const api = {
  Mol
};

module.exports = bindings.load().then(() => api);
