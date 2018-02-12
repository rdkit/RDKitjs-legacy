var RDKit = rdk(Module);

// custom JS methods

RDKit.hello = 'world';


function calc_all_desc(mol){
    var mw = mol.getMW();
    var fr_sp3 = mol.FractionCSP3();
    var ExactMW = mol.ExactMW();
    var Formula = mol.Formula();
    var Chi0v = mol.Chi0v();
    var Chi1v = mol.Chi1v();
    var Chi2v = mol.Chi2v();
    var Chi3v = mol.Chi3v();
    var Chi4v = mol.Chi4v();
    var Chi0n = mol.Chi0n();
    var Chi1n = mol.Chi1n();
    var Chi2n = mol.Chi2n();
    var Chi3n = mol.Chi3n();
    var Chi4n = mol.Chi4n();
    var HallKierAlpha = mol.HallKierAlpha();
    var Kappa1 = mol.Kappa1();
    var Kappa2 = mol.Kappa2();
    var Kappa3 = mol.Kappa3();
    
    var t = new RDKit.VectorDouble();
    t = mol.logp_mr(); 
    var logp = t.get(0);
    var mr = t.get(1);

    var LipinskiHBA = mol.LipinskiHBA();
    var LipinskiHBD = mol.LipinskiHBD();
    var NumRotatableBonds = mol.NumRotatableBonds();
    var NumHBD = mol.NumHBD();
    var NumHBA = mol.NumHBA();
    var NumHeteroatoms = mol.NumHeteroatoms();
    var NumAmideBonds = mol.NumAmideBonds();
    var NumRings = mol.NumRings();
    var NumAromaticRings = mol.NumAromaticRings();
    var NumAliphaticRings = mol.NumAliphaticRings();
    var NumSaturatedRings = mol.NumSaturatedRings();
    var NumHeterocycles = mol.NumHeterocycles();
    var NumAromaticHeterocycles = mol.NumAromaticHeterocycles();
    var NumAromaticCarbocycles = mol.NumAromaticCarbocycles();
    var NumSaturatedHeterocycles = mol.NumSaturatedHeterocycles();
    var NumSaturatedCarbocycles = mol.NumSaturatedCarbocycles();
    var NumAliphaticHeterocycles = mol.NumAliphaticHeterocycles();
    var NumAliphaticCarbocycles = mol.NumAliphaticCarbocycles();
    var LabuteASA = mol.LabuteASA();
    var TPSA = mol.TPSA();

    // get array values
   var slogp = new RDKit.VectorDouble();
    slogp = mol.SlogP_VSA();    
    SlogP_VSA = [];
    for (i =0;i<slogp.size();i++)
    {
        SlogP_VSA.push(slogp.get(i));
    }


    var smr = new RDKit.VectorDouble();
    smr = mol.SMR_VSA();    
    SMR_VSA = [];
    for (i =0;i<smr.size();i++)
    {
        SMR_VSA.push(smr.get(i));
    }


    var peo = new RDKit.VectorDouble();
    peo = mol.PEO_VSA();    
    PEO_VSA = [];
    for (i =0;i<peo.size();i++)
    {
        PEO_VSA.push(peo.get(i));
    }


    var mqn = new RDKit.VectorUint();
    mqn = mol.MQNs();   
    MQNs = [];
    for (i =0;i<mqn.size();i++)
    {
        MQNs.push(mqn.get(i));
    }
   
    return {mw:mw,exactMW:ExactMW,formula:Formula,frsp3:fr_sp3,mqn:MQNs,peovsa:PEO_VSA,smrvsa:SMR_VSA, 
        tpsa:TPSA,slogpvsa:SlogP_VSA,logp:logp,mr:mr,labuteASA:LabuteASA,lipinskiHBD:LipinskiHBD,lipinskiHBA:LipinskiHBA,
        numHeterocycles:NumHeterocycles,numRings:NumRings,numHeteroatoms:NumHeteroatoms,numHBA:NumHBA,
        numHBD:NumHBD,numAliphaticCarbocycles:NumAliphaticCarbocycles,numAliphaticHeterocycles:NumAliphaticHeterocycles,
        numSaturatedCarbocycles:NumSaturatedCarbocycles,numSaturatedHeterocycles:NumSaturatedHeterocycles,
        numAliphaticRings:NumAliphaticRings,numAromaticRings:NumAromaticRings,numSaturatedRings:NumSaturatedRings,
        numAmideBonds:NumAmideBonds, 
        numRotatableBonds:NumRotatableBonds,numAromaticHeterocycles:NumAromaticHeterocycles,
        chi0n:Chi0n,chi1n:Chi1n,chi2n:Chi2n,chi3n:Chi3n,chi4n:Chi4n,
        chi0v:Chi0v,chi1v:Chi1v,chi2v:Chi2v,chi3v:Chi3v,chi4v:Chi4v,kappa1:Kappa1,kappa2:Kappa2,kappa3:Kappa3,
        hallKierAlpha:HallKierAlpha
   }; 
};

/*
// Abraham Descriptors Model from publication : 
var smartsA =['[C][OX2H]','[c][OX2H]','[C][NX3;H2]','[c][NX3;H2;!$(NC=O)]','[C][NX3;H1;!R][C]','[C][NX3;H1;R][C]','[c][NX3;H1;!$(NC=O)][C]','[c][nX3;H1][c]','[CX3](=O)[OX1H0-,OX2H1]','[CX3](=[OX1])[NX3H2]','[CX3](=[OX1])[NX3;H1][C]','[CX3](=[OX1])[NX3;H1][c]','[$([SX4](=[OX1])(=[OX1])([!O])[NH,NH2,NH3+]),$([SX4+2]([OX1-])([OX1-])([!O])[NH,NH2,NH3+])]','[NX3;H1]C(=[OX1])[NX3;H1]','[NX3;H0]C(=[OX1])[NX3;H1]','[NX3;H1]C(=[OX1])O','[NX3;H1]C(=N)[NX3;H0]','[C]#[CH]','P[OH,O-]','[CH][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]','[CH]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]','[CX4]([CX3](=O)[OX1H0-,OX2H1])[CX4][CX3](=O)[OX1H0-,OX2H1]','[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX3](=O)[OX1H0-,OX2H1]','[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[OH]','[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][OH]','[nX3;H1]:n','[nX3;H1]:c:n','[OX2;H1]CC[O,N]','[OX2;H1]C[C,N]=[O,S]','[OX2;H1]c1ccccc1[O,NX3]','[OX2;H1]c1ccccc1C=[O,S]','[OX2;H1]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]','[NH,NH2,NH3+]CC[O,N]','[NH,NH2,NH3+]c1ccccc1[O,N]','[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]','[OX2H]c1ccccc1[Cl,Br,I]','[OX1]=[C,c]~[C,c]C[OH]','[OH]c1cccc2cccnc12','[OH]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1','[OH]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1','[NH,NH2,NH3+]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1','[NH,NH2,NH3+]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1','[CX3](=O)([OX1H0-,OX2H1])c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1','[CX3](=O)([OX1H0-,OX2H1])c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1','[OH]c1c([CX4])cccc1[CX4]','[NH,NH2,NH3+]c1c([CX4])cccc1[CX4]','[OH]c1c(C[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cccc1','[OH]c1cc([CX3](=O)[OX1H0-,OX2H1])ccc1','[OH]c1ccc([CX3](=O)[OX1H0-,OX2H1])cc1','[OH]c1cc([$([CH](=O)),$(C(=O)C)])ccc1','[OH]c1ccc([$([CH](=O)),$(C(=O)C)])cc1'];
var smartsB =['[CX4H3]','[CX4H2]','[CX4H1]','[CX4H0]','*=[CX3H2]','[$(*=[CX3H1]),$([cX3H1](a)a)]','[$(*=[CX3H0]),$([cX3H0](a)(a)A)]','c(a)(a)a','*#C','[C][NX3;H2]','[c][NX3;H2]','[C][NX3;H1][C]','[c][NX3;H1]','[c][nX3;H1][c]','[C][NX3;H0](C)[C]','[c][NX3;H0](C)[C]','[c][nX3;H0][c]','*=[Nv3;!R]','*=[Nv3;R]','[nX2H0,nX3H1+](a)a','N#C[A;!#1]','N#C[a;!#1]',                           '[$([A;!#1][NX3](=O)=O),$([A;!#1][NX3+](=O)[O-])]','[$([a;!#1][NX3](=O)=O),$([a;!#1][NX3+](=O)[O-])]','[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]','[OH]','[OX2;H0;!R]','[OX2;H0;R]','[oX2](a)a','*=O','[SX2](*)*','[sX2](a)a','*=[SX1]','*=[SX3]','[$([#16X4](=[OX1])(=[OX1])([!#8])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([!#8])[OX2H0])]','[S,s]',                                                                                                                                                                                     '[P,p]','FA','Fa','Cl','Br','I','[CX3;!R](=[OX1])[OX2H0]','[CX3;R](=[OX1])[OX2H0;R]','P(=[OX1])(O)(O)O','[CX3](=[OX1])([OX2H0])[OX2H0]','[CX3](=O)[OX1H0-,OX2H1]','nC=[OX1]','[N;!R]C=[OX1]','[N;R][C;R]=[OX1]','[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]','NC(=[OX1])N','[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]','[CX3](=[OX1])[NX3][CX3](=[OX1])','C1(=[OX1])C=CC(=[OX1])C=C1','[$([CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])]','[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]','*1~*2~*(~*3~*(~*~*~*~*3)~*1)~*~*~*1~*2~*~*~*1',                                                                                                                                                                             '[OX2;H1]CC[O,N]','[OX2;H1]C[C,N]=[O,S]','[OX2;H1]c1ccccc1[O,NX3]','[OX2;H1]c1ccccc1C=[O,S]','[OX2;H1]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]','[NH,NH2,NH3+]CC[O,N]','[NH,NH2,NH3+]c1ccccc1[O,N]','[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]','[OX2H]c1ccccc1[Cl,Br,I]','[CX4]([OH])[CX4][OH]','n:n','o:n','n:c:n','o:c:n','n:c:c:n','[F,Cl,Br,I,N,O,S]-c:c-[F,Cl,Br,I,N,O,S]','[F,Cl,Br,I,N,O,S]-c:c:c-[F,Cl,Br,I,N,O,S]','[F,Cl,Br,I,N,O,S]-c:c:c:c-[F,Cl,Br,I,N,O,S]','P(=[OX1])N','Nc:n','[$(cC[OH]);!$(c[CX3](=O)[OX1H0-,OX2H1])]','[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]','[OX2]-c:c-[OX2]'];
var functionA = ['Alipahtic -OH 1','phenol -OH 2','Ali -NH2 3','Aro -NH2 4','NH 5','NH 6','Aniline NH 7','NH pyrrole 8','Acide 9','Amide I 10','Amine II Ali 11','Amine II Aro 12','Thiamide 13','urea 1 14','urea 2 15','carbamate 16','guanidine 17','Alkyne 18','Phosphoric acid 19','RRCHX 20','RCHX2 21','diacid 22','acid special 23','alcohol special 24','special special 25','pyrazole type N 26','imidazole type N 27','H-bond 1 59','H-bond 2 60','H-bond 3 61','H-bond 4 62','H-bond 5 63','H-bond 6 64','H-bond 7 65','H-bond 8 66','H-bond 9 67','H-bond 10 37  37','8-OH quinoline 38','3-X phenol 39','4-X phenol 40','3-X aniline 41 ','4-X aniline 42','3 X benzoic acid 43','4 X benzoic acid 44','2,6 dialkyl phenol 45','2,6 dialkyl aniline 46','2 CX phenol 47','3 CO2H phenol 48','4 CO2H phenol 49','3 C=O phenol 50','4 C=O phenol 51'];
var functionB =['CH3 sp3 1','CH2 sp3 2','CH sp3 3 ','-C- 4','CH2= 5','=CH- 6','=C- 7','C fused aromatic 8','C sp tb 9','Amine 1 sp3 alip 10','Amine 1 sp3 aro 11','Amine 2 sp3 alip 12','Amine 2 sp3 aro 13','HN pyrrole 14','Amine 3 sp3 alip 15','Amine 3 sp3 aro 16','RN pyrrole 17','N sp2 alip 18','N sp2 cyclic 19','Pyridine 20','Nitrile ali 21','Nitrile aro 22','NITRO ali 23','NITRO aro 24','NITRATE 25','any -OH 26','O sp3 noncyclic 27','O sp3 cyclic 28','-O- aromatic 29','=O sp2 30','-S- sp3 31','-S- aromatic 32','S= sp2 33','--S= Sp2 34','sulfonate 35','soufre 36','phosphore 37','F alipathic 38','F aromatic 39','Cl 40','Br 41 ','I 42','Non cycli ester 43','Lactone 44','Phosphate 45','carbonate 46','carboxylique acid 47','Amide Aromatic 48','Amide non cyclic aliphatic 49','Lactam 50','sulfonamide 51','urea 52','carbamate 53','imide 54','quinone 55','-CX2- 56','--XCX-CX-- 57','Steroide 58','H-bond 1 59','H-bond 2 60','H-bond 3 61','H-bond 4 62','H-bond 5 63','H-bond 6 64','H-bond 7 65','H-bond 8 66','H-bond 9 67','1,2 diol 68','pyrazine interaction 69','isoxazole interaction 70','pyrimidine interaction 71','oxazole interaction 72','pyrazine interaction 73','ortho intercation 74','meta intercation 75','para interaction 76','phosphamide 77','2 amino pyridine 78','benzyl alcohol 79','N-oxide 80','1,2-dimethoxy 81'];
var Acoef = [0.345,0.543,0.177,0.247,0.087,0.321,0.194,0.371,0.243,0.275,0.281,-0.091,0.356,-0.165,-0.119,-0.105,0.17,0.082,0.493,0.019,0.05,-0.362,0.118,0.1,0.051,0.194,0.042,-0.089,-0.161,-0.251,-0.418,-0.45,-0.155,0,-0.093,-0.11,-0.601,-0.475,0.119,0.176,0.08,0.084,0.085,0.055,-0.162,-0.181,0.195,-0.203,0.096,0.185,0.203];
var Ecoef = [-0.104,0,0.089,0.187,-0.045,0.068,0.18,0.3,0.04,0.085,0.163,0.138,0.192,-0.03,0.22,0.346,0.083,0.117,0.121,0.046,0,0,0.2,0.21,0,0.061,0.014,0.013,-0.125,-0.041,0.33,0.116,0.364,0.413,0,0.465,0.295,-0.18,-0.23,0.023,0.196,0.533,-0.113,0,-0.1,0,-0.192,0.221,0,0.061,-0.111,-0.11,0,0,0,-0.017,0.012,0.285,0.029,0,-0.069,0,0,0,0,0,-0.1,-0.043,0.092,-0.113,0,0.052,0,0,0,0,-0.08,0.185,0,0,0];
var Scoef = [-0.075,    0,  0.036,  0.071,  -0.085, 0.05,   0.101,  0.121,  0.034,  0.175,  0.383,  0.265,  0.311,  0.221,  0.323,  0.295,  0.265,  0.125,  0.254,  0.223,  0.694,  0.39,   0,  -0.231, -0.476, 0.247,  0.185,  0.185,  0,  0.37,   0.189,  0,  0.618,  1.065,  -0.505, 0.643,  0.703,  -0.042, 0,  0.082,  0.161,  0.198,  -0.225, 0.36,   -0.24,  -0.19,  -0.412, -0.076, 0.175,  -0.1,   -0.569, -0.553, -0.588, -0.51,  -0.411, -0.05,  0,  1.029,  -0.067, -0.095, -0.237, -0.344, -0.276, -0.102, 0,  -0.14,  -0.12,  0.052,  0.024,  0.047,  -0.04,  0.087,  -0.051, -0.043, -0.038, 0,  -0.452, 0.098,  0,  0.434,  0.38];
var BHcoef = [0.007,0,0.011,0.037,0.019,0.011,0,0.019,0.028,0.481,0.275,0.541,0.415,0.316,0.653,0.321,0.392,0.2,0.596,0.321,0.242,0.103,-0.476,-0.525,-0.204,0.307,0.211,0.331,0.047,0.334,0.168,0.043,0.071,0.448,-0.188,0,1.183,-0.036,0,0,-0.011,0,-0.206,-0.214,-0.394,-0.267,-0.308,-0.095,-0.287,-0.231,-0.446,-0.076,-0.252,-0.148,-0.051,-0.014,0.013,0.267,0,-0.068,-0.079,-0.387,-0.126,0,-0.059,-0.045,-0.13,0,-0.132,-0.157,-0.098,-0.17,-0.089,0.031,-0.035,-0.023,-0.668,-0.042,0.131,-0.408,-0.216];
var BOcoef= [0,0,0.02,0.047,0.024,0.012,0,0.018,0.032,0.486,0.326,0.543,0.426,0.267,0.655,0.338,0.338,0.202,0.589,0.3,0.245,0.093,-0.595,-0.533,-0.202,0.311,0.226,0.33,0.06,0.339,0.175,0.083,0.069,0.319,-0.19,0,1.189,-0.033,0,0,0,0,-0.223,-0.169,-0.408,-0.298,-0.312,-0.038,-0.292,-0.242,-0.443,-0.054,-0.251,-0.149,-0.05,-0.016,0.01,0.218,0,-0.09,-0.122,-0.403,-0.12,0,-0.027,-0.069,-0.13,-0.018,-0.094,-0.141,-0.113,-0.184,-0.073,0.025,-0.033,-0.025,-0.668,-0.057,0.129,-0.405,-0.218];
var Lcoef = [0.321,0.499,0.449,0.443,0.244,0.469,0.624,0.744,0.332,0.781,0.949,0.568,0.912,1.25,0.4,0.869,0.794,-0.235,-0.24,0.574,0.757,0.732,0.278,0.347,0,0.672,0.36,0.359,0.057,0.495,1.258,0.848,0.954,2.196,0,0.554,2.051,-0.143,-0.147,0.669,1.097,1.59,-0.39,0.406,-0.483,0,-0.369,0,0.603,0.583,0,0,0,0,0,-0.111,0.054,0.488,-0.072,-0.337,0,-0.303,-0.364,0.062,0,0.169,-0.4,0.1,-0.179,0,0.042,0.209,-0.058,-0.081,-0.026,0,0,0.149,-0.145,0,0];
    

// Abraham ODT model clusters from publication : 
var smartsODT = ['[SX2H1]','[CX3](=O)[OX2H1]', '[$([CX3H][#6]),$([CX3H2])]=[OX1]','[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]',  '[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]']
var functionODT= ['Mercaptan','CARBOXYLIC ACID','ALDEHYDE','ESTER','UNSATURATED'];

// Polarizability, McGowan volume & molecular refraction models 
var functionPol =  ['CH3 sp3 1','CH2 sp3 2','CH sp3 3 ','F desc 5','Cl desc 6','Br desc 7','I desc 8','NITRO 9','N desc 10','O desc 11','S sulfone 12','soufre any 13','phosphore any 14'];
var smartsPol = ['[CX4H3]','[CX4H2]','[CX4H1]','F','Cl','Br','I','[$([NX3](=O)=O),$([NX3+](=O)[O-])]','[N,n]','[O,o]','[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]','[S,s]','[P,p]'];
var Polcoef =  [10.152, 8.765, 5.702, 3.833, 16.557, 24.123, 38.506, 10.488, 6.335, 4.307, 15.726, 22.366, 11.173];



// binary fragments
function fragpresent (smarts, functions, mol) {
  var desc = [];
  for (i=0; i<smarts.length;i++) {
     var v = mol.GetSubstructMatchesNumber(smarts[i]);
     if (v>0) desc[i]=1;
     else desc[i]=0;
    }
  return desc
};


// count fragments
function fragnum(smarts, functions, mol) {
  var desc = [];
  for (i=0; i<smarts.length;i++) {
     var v = mol.GetSubstructMatchesNumber(smarts[i]);
     //console.log(smarts[i],functions[i],v);

     desc[i]=v;
 }
  return desc;
};


function calcV_Pol_MR(atomnumbers,mol) {
  var nC=0, nB=0, nO=0,nN=0,nS=0,nI=0,nCl=0,nBr=0,nP=0,nF=0,nTe=0,nSn=0,nSb=0,nSe=0,nAs=0,nGe=0,nSi=0;
  for (j=0; j<atomnumbers.length;j++)
  {
    if (atomnumbers[j]=== 6)       nC++
    else if (atomnumbers[j]=== 5)  nB++
    else if (atomnumbers[j]=== 8)  nO++
    else if (atomnumbers[j]=== 7)  nN++
    else if (atomnumbers[j]=== 16) nS++
    else if (atomnumbers[j]=== 53) nI++
    else if (atomnumbers[j]=== 17) nCl++
    else if (atomnumbers[j]=== 35) nBr++
    else if (atomnumbers[j]=== 15) nP++
    else if (atomnumbers[j]=== 9)  nF++
    else if (atomnumbers[j]=== 52) nTe++
    else if (atomnumbers[j]=== 50) nSn++
    else if (atomnumbers[j]=== 51) nSb++
    else if (atomnumbers[j]=== 34) nSe++
    else if (atomnumbers[j]=== 33) nAs++
    else if (atomnumbers[j]=== 32) nGe++
    else if (atomnumbers[j]=== 14) nSi++
  }

  // count number of Hydrogens
  mol.removeHs();
  var nH1 = mol.getNumAtoms();
  mol.addHs();
  var nH2 =  mol.getNumAtoms();
  var nH = nH2-nH1;   
  //console.log("number of H:",nH);
  var nRing = mol.NumRings();
  //console.log("number of Rings:",nRing);





  // Abraham MH, McGowan JC., Chromatographia. 1987;23:243-6.
  var Sumx = 0.0871*nH +0.1635*nC + 0.1439*nN + nO*0.1243 + nP*0.2487 + nS*0.2291 + nF*0.1048 + nCl*0.2095 + nBr*0.2621 + nI*0.3453;
  Sumx = Sumx+ 0.2683*nSi + 0.1832*nB + 0.3102*nGe +0.2942*nAs + 0.2781*nSe + 0.3935*nSn + 0.3774*nSb + 0.3614*nTe ;
  // count number of bond
  var nBond = nH +nC + nN + nO+ nP+ nS + nF + nCl+ nBr+ nI+nSi + nB + nGe +nAs + nSe + nSn + nSb + nTe - 1 + nRing;
  // McGowan Volume 
  var V= Sumx-0.0656*nBond;

  //console.log("McGowan Volume:",V);

  // Polarizability 
  var Pol = -1.529;
  var Polfinger =fragnum(smartsPol, functionPol, mol);
  Polfinger[12]=Polfinger[12]-Polfinger[11];  // remove the sulfones in nS count!
  Polfinger[9]=Polfinger[9]-Polfinger[8]; // remove the Nitro in nN count!


  for (j=0;j<Polfinger.length;j++) 
  { 
    if (Polfinger[j]>0) { Pol += Polcoef[j]*Polfinger[j];
      //console.log(j,Pol,Polcoef[j],Polfinger[j]);
    }
  }

  // Polarizability
  //console.log("Pol:",Pol);

  // Molecular refraction
  var MR=4/3*Math.PI*Pol;
  //console.log("MR:",MR);

  return {nC,nB,nO,nN,nS,nI,nCl,nBr,nP,nF,nTe,nSn,nSb,nSe,nAs,nGe,nSi,nH,nRing,V,Pol,MR};
};
  

function Abrahams(smartsA, functionA, smartsB, functionB, smartsODT,functionODT,mol){

var Afinger = fragnum(smartsA, functionA, mol);
var Bfinger = fragnum(smartsB, functionB, mol);
var ODTcluster = fragpresent(smartsODT, functionODT, mol);
var A = 0.003;
var Bh= 0.071;
var Bo= 0.064;
var L = 0.13;
var S = 0.277;
var E = 0.24;

for (j=0;j<Afinger.length;j++) 
{ 
  if (Afinger[j]>0) { A += Acoef[j]*Afinger[j];
    //console.log(j,Aval,Acoef[j],Afinger[j]);
  }
}

// need to reduce FRAGMENT 36 proportion based on others sulphur contributions in : (31 OR 32 OR 33 OR 34 OR 35)
Bfinger[35]= Bfinger[35]-(Bfinger[30]+Bfinger[31]+Bfinger[32]+Bfinger[33]+Bfinger[34]);


for (k=0;k<Bfinger.length;k++) 
{
  if (Bfinger[k]>0) {
    Bh += BHcoef[k]*Bfinger[k];
    Bo += BOcoef[k]*Bfinger[k];
    L += Lcoef[k]*Bfinger[k];
    S += Scoef[k]*Bfinger[k];
    E += Ecoef[k]*Bfinger[k];
  }
}

//  logODT model linear from Abraham publication
var M = ODTcluster[0];
var AC = ODTcluster[1];
var AL = ODTcluster[2];
var UE = ODTcluster[3]*ODTcluster[4];
var logODT = 1.826- 0.882*E-0.408*S -0.999*A-2.196*Bh-0.578*L-4.065*M-1.805*AL-1.424*AC-1.290*UE; 


return {A,Bh,Bo,L,S,E,logODT};
};

RDKit.Abrahams = Abrahams;
RDKIt.calcV_Pol_MR = calcV_Pol_MR;

*/


RDKit.calc_all_desc = calc_all_desc;



// END custom JS methods

return RDKit;
});
