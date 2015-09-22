from rdkit import Chem
import cPickle,sys
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
import time
start_time = time.time()

nm=[x[0] for x in Descriptors._descList]

calc = MoleculeDescriptors.MolecularDescriptorCalculator(nm)
suppl = Chem.SmilesMolSupplier('/Users/mbp/Documents/abraham.txt',titleLine=False)

w = open('rdkitdesc.txt',"wb")
errors = open('rdkiterror.txt',"wb")

line = "Molid"
for names in nm:
	line = line+ ","+ names 
line = line+"\n";
w.write(line)
nDone=0
for mol in suppl: 
    nDone += 1
    line = str(nDone)    
    if mol is None: 
    	errors.write(line + ",Error\n" )
    	continue
    descrs = calc.CalcDescriptors(mol)
    for dec in descrs:
        line = line+ "," +str(dec)
    line = line+"\n";
    w.write(line)
w.close()
errors.close()
print("--- %s seconds ---" % (time.time() - start_time))
