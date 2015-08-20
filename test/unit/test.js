'use strict';

var RDKit = require('../rdkit');

describe('test', function () {
    it('should work', function () {
        RDKit.hello.should.equal('world');
    });
});