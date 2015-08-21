'use strict';

var RDKit = require('../rdkit');

describe('test', function () {
    it('should work', function () {
        RDKit.hello.should.equal('world');
    });
});


describe('calc_all_desc', function () {
    it('calc_all_desc', function () {
		var smi = 'CCCCCOC(CO)';
		var mol = RDKit.Molecule.fromSmiles(smi);	
		var c = RDKit.calc_all_desc(mol);
        c.should.equal('');
    });
});


describe('calc_all_desc', function () {
    it('calc_all_desc', function () {
		var smi = 'CCCCCOC(CO)';
		var mol = RDKit.Molecule.fromSmiles(smi);	
		var c = RDKit.calc_all_desc(mol);
        c.should.equal('');
    });
});

