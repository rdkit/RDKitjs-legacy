

var fs = require('fs');
var RDKit = require('rdkit');


var smartsA ={'[C][OX2H]','[c][OX2H]','[C][NX3;H2]','[c][NX3;H2;!$(NC=O)]','[C][NX3;H1;!R][C]','[C][NX3;H1;R][C]','[c][NX3;H1;!$(NC=O)][C]','[c][nX3;H1][c]','[CX3](=O)[OX1H0-,OX2H1]','[CX3](=[OX1])[NX3H2]','[CX3](=[OX1])[NX3;H1][C]','[CX3](=[OX1])[NX3;H1][c]','[$([SX4](=[OX1])(=[OX1])([!O])[NH,NH2,NH3+]),$([SX4+2]([OX1-])([OX1-])([!O])[NH,NH2,NH3+])]','[NX3;H1]C(=[OX1])[NX3;H1]','[NX3;H0]C(=[OX1])[NX3;H1]','[NX3;H1]C(=[OX1])O','[NX3;H1]C(=N)[NX3;H0]'),'[C]#[CH]','P[OH,O-]','[CH][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]','[CH]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]','[CX4]([CX3](=O)[OX1H0-,OX2H1])[CX4][CX3](=O)[OX1H0-,OX2H1]','[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX3](=O)[OX1H0-,OX2H1]','[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[OH]','[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][OH]','[nX3;H1]:n','[nX3;H1]:c:n','[OX2;H1]CC[O,N]','[OX2;H1]C[C,N]=[O,S]','[OX2;H1]c1ccccc1[O,NX3]','[OX2;H1]c1ccccc1C=[O,S]','[OX2;H1]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]','[NH,NH2,NH3+]CC[O,N]','[NH,NH2,NH3+]c1ccccc1[O,N]','[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]','[OX2H]c1ccccc1[Cl,Br,I]','[OX1]=[C,c]~[C,c]C[OH]','[OH]c1cccc2cccnc12','[OH]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1','[OH]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1','[NH,NH2,NH3+]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1','[NH,NH2,NH3+]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1','[CX3](=O)([OX1H0-,OX2H1])c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1','[CX3](=O)([OX1H0-,OX2H1])c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1','[OH]c1c([CX4])cccc1[CX4]','[NH,NH2,NH3+]c1c([CX4])cccc1[CX4]','[OH]c1c(C[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cccc1','[OH]c1cc([CX3](=O)[OX1H0-,OX2H1])ccc1','[OH]c1ccc([CX3](=O)[OX1H0-,OX2H1])cc1','[OH]c1cc([$([CH](=O)),$(C(=O)C)])ccc1','[OH]c1ccc([$([CH](=O)),$(C(=O)C)])cc1'};

var functionA = {'Alipahtic -OH 1','phenol -OH 2','Ali -NH2 3','Aro -NH2 4','NH 5','NH 6','Aniline NH 7','NH pyrrole 8','Acide 9','Amide I 10','Amine II Ali 11','Amine II Aro 12','Thiamide 13','urea 1 14','urea 2 15','carbamate 16','guanidine 17','Alkyne 18','Phosphoric acid 19','RRCHX 20','RCHX2 21','diacid 22','acid special 23','alcohol special 24','special special 25','pyrazole type N 26','imidazole type N 27','H-bond 1 59','H-bond 2 60','H-bond 3 61','H-bond 4 62','H-bond 5 63','H-bond 6 64','H-bond 7 65','H-bond 8 66','H-bond 9 67','H-bond 10 37  37','8-OH quinoline 38','3-X phenol 39','4-X phenol 40','3-X aniline 41 ','4-X aniline 42','3 X benzoic acid 43','4 X benzoic acid 44','2,6 dialkyl phenol 45','2,6 dialkyl aniline 46','2 CX phenol 47','3 CO2H phenol 48','4 CO2H phenol 49','3 C=O phenol 50','4 C=O phenol 51'};

var Acoef = [0.345,0.543,0.177,0.247,0.087,0.321,0.194,0.371,0.243,0.275,0.281,-0.091,0.356,-0.165,-0.119,-0.105,0.17,0.082,0.493,0.019,0.05,-0.362,0.118,0.1,0.051,0.194,0.042,-0.089,-0.161,-0.251,-0.418,-0.45,-0.155,0,-0.093,-0.11,-0.601,-0.475,0.119,0.176,0.08,0.084,0.085,0.055,-0.162,-0.181,0.195,-0.203,0.096,0.185,0.203];

var desc = [];
function A(smartsA, functionA, mol) {
  for (i=0; i<smarts.length;i++) { i < NB_fonction:
     var func = RDKit.Molecule.MolFromSmarts(smartsA[i]);
     var v = mol.GetSubstructMatchesNumber(func);
     desc[i]=v;
  return desc
}



smi = 'O=[N+]([O-])c1c(cc(O)cc1)C';
var mol = RDKit.Molecule.fromSmiles(smi);	
var Afinger =A(smartsA, functionA, mol);
var Aval = 0.003;
for (j=0;j<Afinger.length;j++) 
{
	if (Afinger[j]>0) Aval += Acoef[j]*Afinger[j];

}
console.log(Aval);
