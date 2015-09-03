'use strict';

var RDKit = require('../rdkit');

describe('RDKit loading', function () {
    it('JS functions Hello', function () {
        RDKit.hello.should.equal('world');

    });
});

    // static constructors
//    newmolecule();
describe('Molecule Creations', function () {
    it('Murko', function () {
        var smi = 'C1CCCC1OC(CO)';
        var murcko = RDKit.Molecule.MurckofromSmiles(smi);
        murcko.sdwrite().should.equal('\n     RDKit          \n\n  5  5  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  1  1  0\nM  END\n$$$$\n');
        murcko.delete();
    });


    it.skip('Mol2BlockToMol', function () {
        var molBlock = '';
        var mol = RDKit.Molecule.Mol2BlockToMol(molBlock);
        mol.sdwrite().should.equal();
        mol.delete();
    });

    it('MolBlockToMol', function () {
        var molBlock = '\n     RDKit          \n\n  9  9  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\n  5  1  1  0\nM  END\n$$$$\n';
        var mol = RDKit.Molecule.MolBlockToMol(molBlock);
        mol.smilewrite().should.equal('OCCOC1CCCC1 0\n');
        mol.delete();

    });

    it('fromSmarts', function () {
        var smarts ='[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]';
        var mol = RDKit.Molecule.fromSmarts(smarts);
        mol.sdwrite().should.equal('\n     RDKit          \n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  6  0\n  2  3  2  0\n  2  4  6  0\nV    1 [N&X3,N&X4&+]\nV    2 [C&X3]\nV    3 [O&X1]\nV    4 [O&X2&H1,O&X1&-]\nM  END\n$$$$\n');
        mol.delete();
    });

  it('fromSmiles', function () {
        var smi = 'C1CCCC1OC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);
        mol.sdwrite().should.equal('\n     RDKit          \n\n  9  9  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\n  5  1  1  0\nM  END\n$$$$\n');
        mol.delete();
   });

    it('molFromPickle', function () {
        var pickle='ï¾­Þ\u0000\u0000\u0000\u0000\u0007\u0000\u0000\u0000\u0002\u0000\u0000\u0000\u0000\u0000\u0000\u0000\t\u0000\u0000\u0000\t\u0000\u0000\u0000\u0001\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0003\u0001\b\u0000 \u0000\u0000\u0000\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\b\u0000`\u0000\u0000\u0000\u0001\u0001\u000b\u0000\u0001\u0000\u0001\u0002\u0000\u0002\u0003\u0000\u0003\u0004\u0000\u0004\u0005\u0000\u0005\u0006\u0000\u0006\u0007\u0000\u0007\b\u0000\u0004\u0000\u0000\u0014\u0001\u0005\u0000\u0001\u0002\u0003\u0004\u0017\u0000\u0000\u0000\u0000\u0016';
        var mol = RDKit.Molecule.molFromPickle(pickle);
        mol.sdwrite().should.equal('\n     RDKit          \n\n  9  9  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\n  5  1  1  0\nM  END\n$$$$\n');
        mol.delete();
    });


        // Pickle molecule representation
    it('MolToBinary', function () {
         var smi = 'C1CCCC1OC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);
        var p= mol.MolToBinary();
        p.should.equal('ï¾­Þ\u0000\u0000\u0000\u0000\u0007\u0000\u0000\u0000\u0002\u0000\u0000\u0000\u0000\u0000\u0000\u0000\t\u0000\u0000\u0000\t\u0000\u0000\u0000\u0001\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0003\u0001\b\u0000 \u0000\u0000\u0000\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\u0006\u0000`\u0000\u0000\u0000\u0002\u0002\b\u0000`\u0000\u0000\u0000\u0001\u0001\u000b\u0000\u0001\u0000\u0001\u0002\u0000\u0002\u0003\u0000\u0003\u0004\u0000\u0004\u0005\u0000\u0005\u0006\u0000\u0006\u0007\u0000\u0007\b\u0000\u0004\u0000\u0000\u0014\u0001\u0005\u0000\u0001\u0002\u0003\u0004\u0017\u0000\u0000\u0000\u0000\u0016');
        mol.delete();
    });

});



describe('Compute all descriptors', function () {
    it('calc_all_desc', function () {
		var smi = 'CCCCCOC(CO)';
		var mol = RDKit.Molecule.fromSmiles(smi);	
		var c = RDKit.calc_all_desc(mol);
        c=JSON.stringify(c);
        c.should.eql('{"mw":132,"exactMW":132.115029752,"formula":"C7H16O2","frsp3":1,"mqn":[7,0,0,0,0,0,0,0,0,2,0,0,8,0,0,0,0,0,6,4,2,1,1,0,0,2,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"peovsa":[9.843390348640755,0,0,0,0,19.765380445542643,6.4208216229260096,6.606881964512918,13.213763929025836,0,0,0,0,0],"smrvsa":[9.843390348640755,0,0,0,26.186202068468653,19.820645893538753,0,0,0,0],"tpsa":29.46,"slogpvsa":[0,24.927173288379457,4.736862953800049,0,26.186202068468653,0,0,0,0,0,0,0],"logp":1.1855,"mr":37.42979999999999,"labuteASA":56.83668269591208,"lipinskiHBD":1,"lipinskiHBA":2,"numHeterocycles":0,"numRings":0,"numHeteroatoms":2,"numHBA":2,"numHBD":1,"numAliphaticCarbocycles":0,"numAliphaticHeterocycles":0,"numSaturatedCarbocycles":0,"numSaturatedHeterocycles":0,"numAliphaticRings":0,"numAromaticRings":0,"numSaturatedRings":0,"numAmideBonds":0,"numRotatableBonds":6,"numAromaticHeterocycles":0,"chi0n":6.098102573083107,"chi1n":3.6006848163930116,"chi2n":2.043086014632321,"chi3n":1.1278531854030205,"chi4n":0.6207359402846874,"chi0v":6.098102573083107,"chi1v":3.6006848163930116,"chi2v":2.043086014632321,"chi3v":1.1278531854030205,"chi4v":0.6207359402846874,"kappa1":8.92,"kappa2":7.919999999999998,"kappa3":7.920000000000001,"hallKierAlpha":-0.08}');
        mol.delete();
    });

/*    it('computeGasteigerCharges', function () {

        var smi="CCCCC(C)C";
        var mol = RDKit.Molecule.fromSmiles(smi);
        mol.computeGasteigerCharges();
        mol.delete();

    });*/
});


describe('FingerPrints', function () { 
    it('getRDKFP', function () {
        var smi = 'CCCCCOC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);   
        mol.getRDKFP().should.equal('00000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000100000000000000000001000000000000000000000000001000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000100000000000000000000000000010000000000001000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000010001000001000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000100000000000000000000000000100000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000010000000000001000000000000000000000000000000000000000000000100000100000000010000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000010000000000000000000000000000000100000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000100000100000000000000000000000000000000001000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000001000000000000010000000000000010001000000010000000000000000000000000000000000000000000000000000000000000000000000000000000');
        mol.delete();
    });

    it('getMorganFP', function () {
        var smi = 'CCCCCOC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);   
        mol.getMorganFP(2,2048).should.equal('00000000000001000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000');
        mol.delete();
    });


    it('getMorganFP_GetOnBits', function () {
        var smi = 'CCCCCOC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);   
        var mpf2 = mol.getMorganFP_GetOnBits(2,2048);
        var mpf = [];
        var mpflen = mpf2.size();
        for(var j=0;j < mpflen;j++){
                mpf.push(mpf2.get(j));
            }
        mpf.should.eql([ 13, 80, 222, 294, 473, 591, 691, 695, 794, 807, 1057, 1069, 1249, 1325, 1369, 1444, 1911, 1992 ]);
        mol.delete();
    });
    

    it('getLayeredFP', function () {
        var smi = 'CCCCCOC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);   
        mol.getLayeredFP(2,2,2048).should.equal('00000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000');
        mol.delete();
    });


    it('getMACCSFP', function () {
        var smi = 'CCCCCOC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);   
        mol.getMACCSFP().should.equal('00000000000000000000000000000000000000000000000000000000000000000000000010000000001000100001000000000000100011000011101000000010110010000011000000010000010101011000100');
        mol.delete();
    });

   
    it('getPatternFP', function () {
        var smi = 'CCCCCOC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);   
        mol.getPatternFP().should.equal('00000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000100001000000000000000000000000000000000010000000000000000000000000000000000010001000000000000000000001000000000000000000000100000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100001000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000010000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000001000000000000000000000000100000000000000000000001000000000000000010001000010000000000000000000000000000001000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001011000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000');
        mol.delete();
    });
});

    
    
    
    // 3D Force Field minimization
describe('3D Force Field minimization', function () {
    it('MMFFoptimizeMolecule', function () {
        var smi = 'CCCCCOC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);  
        mol.addHs(); 
        mol.EmbedMolecule();
        mol.MMFFoptimizeMolecule();
        mol.removeHs()
        mol.sdwrite().should.equal('\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n    3.1421    1.8386   -0.4923 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6872    0.3939   -0.6248 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2560    0.2090   -0.1205 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8032   -1.2480   -0.2553 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6112   -1.4612    0.2765 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5285   -0.7132   -0.5206 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8815   -0.8099   -0.0574 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1425    0.2334    1.0272 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8418    1.5355    0.5278 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\nM  END\n$$$$\n');
        mol.delete();
    });

    it('MMFFoptimizeMolecule with paramaters', function () {
        var smi = 'CCCCCOC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);  
        mol.addHs(); 
        mol.EmbedMolecule();    
        mol.MMFFoptimizeMoleculearg(1000, 'MMFF94');
        mol.removeHs()
        mol.sdwrite().should.equal('\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n   -3.8778    0.3120   -0.5134 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0248   -0.5928    0.3621 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6103   -0.0585    0.5998 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7594   -0.0088   -0.6737 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6675    0.4592   -0.3946 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3211   -0.4950    0.4391 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6477   -0.1175    0.8000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.6604   -0.4181   -0.3038 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.9785   -0.0672    0.1025 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\nM  END\n$$$$\n');
        mol.delete();
    });


    it('MMFFOptimizeMoleculeConfs', function () { 
        var smi = 'CCCCCOC(CO)';
        this.timeout(50000);

        for (i=0;i<50;i++){
            var mol = RDKit.Molecule.fromSmiles(smi);  
            mol.addHs(); 
            mol.EmbedMultipleConfsarg(2,200,2015);   
            mol.MMFFOptimizeMoleculeConfs(8,200,'MMFF94');
            process.stdout.write('.');
            mol.delete();
        }
        mol = RDKit.Molecule.fromSmiles(smi);  
        mol.addHs(); 
        mol.removeHs()
        mol.EmbedMultipleConfsarg(3,1000,2015);   
        mol.MMFFOptimizeMoleculeConfs(8,1000,'MMFF94');
        mol.sdwriteConfs().should.equal('\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n    4.2049   -0.1343   -0.0700 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.7636    0.1674   -0.4518 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8530    0.1802    0.7765 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4134    0.5914    0.4625 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3280   -0.4062   -0.4258 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6681    0.0257   -0.6665 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.5743   -0.3641    0.3649 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.9713    0.1350    0.0157 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.8966   -0.2437    1.0272 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\nM  END\n$$$$\n\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n   -3.8544   -0.0601   -0.7140 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0532    1.0931   -0.1298 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7683    0.6501    0.5730 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7396    0.0371   -0.3759 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5628   -0.2933    0.3465 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4641   -0.8645   -0.6010 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.7225   -1.2189   -0.0276 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.6792   -0.0304    0.0174 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.9275   -0.4269    0.5719 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\nM  END\n$$$$\n\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n   -3.4103    0.9186   -0.5602 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9693    0.3926    0.7969 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9105   -0.7087    0.7105 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5777   -0.2243    0.1392 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4764   -1.3299    0.1717 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6805   -0.9212   -0.4768 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5109   -0.0983    0.3386 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8150    0.1867   -0.3998 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6449    1.0642    0.3526 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\nM  END\n$$$$\n');
        mol.delete();
    });


    it('UFFOptimizeMolecule', function () {
        var smi = 'CCCCCOC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);  
        mol.addHs(); 
        mol.EmbedMolecule(); 
        mol.UFFOptimizeMolecule();
        mol.removeHs()
        mol.sdwrite().should.equal('\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n   -3.2713    0.4667   -1.0876 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9791   -0.8786   -0.4235 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8749   -0.7887    0.6411 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4844   -0.5591    0.0307 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6043   -0.5961    1.1078 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8803   -0.4433    0.5144 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3128    0.8990    0.6042 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.6179    1.0667   -0.1689 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6397    0.3074    0.4191 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\nM  END\n$$$$\n');
        mol.delete();
    });


    it('EmbedMolecule', function () {
        var smi = 'CCCCCOC(CO)';
        var mol = RDKit.Molecule.fromSmiles(smi);  
        mol.addHs(); 
        mol.EmbedMoleculearg(1000,1);
        mol.removeHs()
        mol.sdwrite().should.equal('\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n   -3.4616    1.4289    0.4314 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9477    0.0284    0.1317 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4668    0.0949   -0.1176 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8657   -1.2686    0.0740 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6165   -1.2723   -0.1495 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2728   -0.1862    0.4022 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6359   -0.3660    0.1169 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.3801    0.9281    0.0006 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5567    0.7909    0.7537 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  1  0\n  8  9  1  0\nM  END\n$$$$\n');
        mol.delete();
    });

});    





/*
    
    

    
    
  
describe('', function () {
    it('should work', function () {
        RDKit.  getPath();

    });
});

      describe('', function () {
    it('should work', function () {

    smilewrite();

    });
});

  describe('', function () {
    it('should work', function () {
sdwriteConfs();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.
    compute2DCoords();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.    
    Drawing2D();
       });
}); 
    // similarity

describe('', function () {
    it('should work', function () {
        RDKit.
    TanimotoSimilarityfromSmile (smilesref);
    });
});
describe('', function () {
    it('should work', function () {
        RDKit.
    DiceSimilarityfromSmile (smilesref);
    });
});
describe('', function () {
    it('should work', function () {
        RDKit.
    TverskySimilarityfromSmile( smilesref,a, b);
        });
});
describe('', function () {
    it('should work', function () {
        RDKit.
    findSSSR(res);
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.
    AlignMolConformers();
    });
});
describe('', function () {
    it('should work', function () {
        RDKit.
    AlignMol(smilesref);
        });
});

    // molecule manipulation & cleaning, ...
describe('', function () {
    it('should work', function () {
        RDKit.
    addHs();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.
    removeHs();
    });
});


  describe('', function () {
    it('should work', function () {
sanitizeMol();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.
    cleanUp();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.
    Kekulize();
    });
});

    
    // descriptors
    describe('', function () {
    it('should work', function () {
describe('', function () {
    it('should work', function () {
        RDKit.  getMW();

    });
});

describe('', function () {
    it('should work', function () {
        RDKit.ExactMW();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Formula();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Chi0v();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Chi1v();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Chi2v();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Chi3v();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Chi4v();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Chi0n();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Chi1n();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Chi2n();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Chi3n();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Chi4n();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.HallKierAlpha();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Kappa1();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Kappa2();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.Kappa3();
    });
});

describe('', function () {
    it('should work', function () {
        RDKit.logp_mr();
    });
});

    
    
describe('', function () {
    it('should work', function () {
    LipinskiHBA();
});
});

describe('', function () {
    it('should work', function () {
    LipinskiHBD();
});
});

describe('', function () {
    it('should work', function () {
    NumRotatableBonds();
});
});

describe('', function () {
    it('should work', function () {
    NumHBD();
});
});

describe('', function () {
    it('should work', function () {
    NumHBA();
});
});

describe('', function () {
    it('should work', function () {
    NumHeteroatoms();
});
});

describe('', function () {
    it('should work', function () {
    NumAmideBonds();
});
});

    FractionCSP3();
describe('', function () {
    it('should work', function () {
    NumRings();
});
});

describe('', function () {
    it('should work', function () {
    NumAromaticRings();
});
});

describe('', function () {
    it('should work', function () {
    NumAliphaticRings();
});
});

describe('', function () {
    it('should work', function () {
    NumSaturatedRings();
});
});

describe('', function () {
    it('should work', function () {
    NumHeterocycles();
});
});

describe('', function () {
    it('should work', function () {
    NumAromaticHeterocycles();
});
});

describe('', function () {
    it('should work', function () {
    NumAromaticCarbocycles ();
});
});

describe('', function () {
    it('should work', function () {
    NumSaturatedHeterocycles();
});
});

describe('', function () {
    it('should work', function () {
    NumSaturatedCarbocycles();
});
});

describe('', function () {
    it('should work', function () {
    NumAliphaticHeterocycles();
});
});

describe('', function () {
    it('should work', function () {
    NumAliphaticCarbocycles();
});
});

describe('', function () {
    it('should work', function () {
    LabuteASA();
});
});

    
    TPSA();
  describe('', function () {
    it('should work', function () {
SlogP_VSA();
});
});
  describe('', function () {
    it('should work', function () {
SMR_VSA();
});
});

   describe('', function () {
    it('should work', function () {
    PEO_VSA();
});
});

describe('', function () {
    it('should work', function () {
        RDKit.
    MQNs();
});
});


describe('', function () {
    it('should work', function () {
        RDKit.  getSubstructMatches(smilesref);
    });
});


describe('', function () {
    it('should work', function () {

    HasSubstructMatchStr(smilesref);
    });
});    

  

describe('', function () {
    it('should work', function () {
        RDKit.  getProp(key);
    });
});



  describe('', function () {
    it('should work', function () {
setProp(key, value);
});
});
    describe('', function () {
    it('should work', function () {
describe('', function () {
    it('should work', function () {
        RDKit.  getNumAtoms();

    });
});

    describe('', function () {
    it('should work', function () {
describe('', function () {
    it('should work', function () {
        RDKit.  getNumConformers();
    });
});


describe('', function () {
    it('should work', function () {
        RDKit.  getConformer(id);
    });
});


 describe('', function () {
    it('should work', function () {
     hasProp(key);
    });
});   


describe('', function () {
    it('should work', function () {
        RDKit.  getproplist();
    });
});


    
    
    
    // atom & bond manipulations
  describe('', function () {
    it('should work', function () {
    addAtom (atomid);
    });
});
    // this is in development stage caution not working for the moment!!!!
    describe('', function () {
    it('should work', function () {
    addBond (beginAtomIdx, endAtomIdx,bondtypeid);
    });
});


  describe('', function () {
    it('should work', function () {
setBondDir (Bondid, bonddirid);
});
});
    

*/
