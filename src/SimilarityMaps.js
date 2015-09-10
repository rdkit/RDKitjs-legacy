var num = require('numeric');
var RDKit = require('rdkit');


function bivariate_normal(X,Y,sigma,mux,muy)
{
	var lensize = num.dim(X)[0];
	var rsize = [lensize,lensize];
	var Xmu = num.sub(X,mux);
	var Ymu = num.sub(Y,muy);
 	var le = num.div(num.pow(Xmu,2),sigma*sigma);
 	var re = num.div(num.pow(Ymu,2),sigma*sigma); 
	var z = num.add(le,re);
	var denom = 2*Math.PI*sigma*sigma;
	var res = num.div(num.exp(num.neg(num.div(z,2))),denom);
	return res;
}

function calcAtomGaussians(sigma,step,weigths,px,py)
{   
	var x = num.linspace(0,1,51);     // default step 0.02 => 51 points
	var v = [];
	var slen = x.length;
	var X = num.rep([slen,],x);
	var Y = num.transpose(X); // meshgrid square grid X=Y'; 
	var Z = num.rep([slen,slen],0);
	for (j=0;j<px.length;j=j+1){   
		if (Math.abs(weigths[j])>0){
			v = bivariate_normal(X,Y,sigma/10,px[j],py[j]); // rescale also the sigma value... caution
			console.log("vsum:",num.sum(v));
			Z = num.add(Z,num.mul(v,weigths[j]));
			console.log("Zsum:",num.sum(Z));
		}
	}
	return {X,Y,Z};
}

var m=RDKit.Molecule.fromSmiles('COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21');
var w = m.getTPSAAtomContribs();
var ap = m.getAtomsPos2D();

var px= [];
var py= [];
var weigths = [];
for  (i=0;i<ap.size()-1;i=i+2)
{	
	weigths[i/2]=w.get(i/2);
	px[i/2]=(ap.get(i));  // rescaling = center origin to left corner origin & dpi factor 
	py[i/2]=(ap.get(i+1)); 
}


xmin = Math.min.apply(Math, px);
xmax = Math.max.apply(Math, px);

ymin = Math.min.apply(Math, py);
ymax = Math.max.apply(Math, py);
// minx = arr.reduce(function(a, b, i, px) {return Math.min(a,b)});
console.log("xmin:",xmin,",xmax:",xmax);
console.log(px);
console.log("ymin:",ymin,",ymax:",ymax);
console.log(py);


// need to find why there is an issue to retreive the ap.get(t-1) value ?
var simmap=calcAtomGaussians(0.45,0.02,weigths,px,py);

console.log(JSON.stringify(simmap));







/*def bivariate_normal(X, Y, sigmax=1.0, sigmay=1.0, mux=0.0, muy=0.0, sigmaxy=0.0):
    Xmu = X-mux
    Ymu = Y-muy
    rho = sigmaxy/(sigmax*sigmay)
    z = Xmu**2/sigmax**2 + Ymu**2/sigmay**2 - 2*rho*Xmu*Ymu/(sigmax*sigmay)
    denom = 2*np.pi*sigmax*sigmay*np.sqrt(1-rho**2)
    return np.exp(-z/(2*(1-rho**2))) / denom


def MolToMPL(mol,size=(300,300),kekulize=True, wedgeBonds=True, imageType=None, fitImage=False, options=None, **kwargs):
  if not mol:
    raise ValueError('Null molecule provided')
  from rdkit.Chem.Draw.mplCanvas import Canvas
  canvas = Canvas(size)
  if options is None:
    options = DrawingOptions()
    options.bgColor=None
  if fitImage:
      drawingOptions.dotsPerAngstrom = int(min(size) / 10)
  options.wedgeDashedBonds=wedgeBonds
  drawer = MolDrawing(canvas=canvas, drawingOptions=options)
  omol=mol
  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)
    
  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)
  
  drawer.AddMol(mol,**kwargs)
  omol._atomPs=drawer.atomPs[mol]
  for k,v in iteritems(omol._atomPs):
    omol._atomPs[k]=canvas.rescalePt(v)
  canvas._figure.set_size_inches(float(size[0])/100,float(size[1])/100)
  return canvas._figure

  dpi = 300
  def rescalePt(self,p1): 
      return [float(p1[0])/self._dpi,float(self.size[1]-p1[1])/self._dpi] 



     class Canvas(CanvasBase): 
    def __init__(self, size, name='', imageType='png'): 
       self._name = name 
      self.size=size 
 23      dpi = max(size[0],size[1]) 
 24      figsize=(int(float(size[0])/dpi),int(float(size[1])/dpi)) 
 25      self._figure = figure(figsize=figsize) 
 26      self._axes = self._figure.add_axes([0,0,2.5,2.5]) 
 27      self._axes.set_xticklabels('') 
 28      self._axes.set_yticklabels('') 
 29      self._dpi = dpi 


def calcAtomGaussians(mol,a=0.03,step=0.02,weights=None):
  import numpy
  from matplotlib import mlab
  x = numpy.arange(0,1,step)
  y = numpy.arange(0,1,step)
  X,Y = numpy.meshgrid(x,y)
  if weights is None:
    weights=[1.]*mol.GetNumAtoms()
  Z = mlab.bivariate_normal(X,Y,a,a,mol._atomPs[0][0], mol._atomPs[0][1])*weights[0] # this is not bivariate case ... only univariate no mixtures #matplotlib.mlab.bivariate_normal(X, Y, sigmax=1.0, sigmay=1.0, mux=0.0, muy=0.0, sigmaxy=0.0)
  for i in range(1,mol.GetNumAtoms()):
    Zp = mlab.bivariate_normal(X,Y,a,a,mol._atomPs[i][0], mol._atomPs[i][1])
    Z += Zp*weights[i]
  return X,Y,Z

def GetSimilarityMapFromWeights(mol, weights, colorMap=cm.PiYG, scale=-1, size=(250, 250), sigma=None,  #@UndefinedVariable  #pylint: disable=E1101
                                coordScale=1.5, step=0.01, colors='k', contourLines=10, alpha=0.5, **kwargs):
  if mol.GetNumAtoms() < 2: raise ValueError("too few atoms")
  fig = Draw.MolToMPL(mol, coordScale=coordScale, size=size, **kwargs)
  if sigma is None:
    if mol.GetNumBonds() > 0:
      bond = mol.GetBondWithIdx(0)
      idx1 = bond.GetBeginAtomIdx()
      idx2 = bond.GetEndAtomIdx()
      sigma = 0.3 * math.sqrt(sum([(mol._atomPs[idx1][i]-mol._atomPs[idx2][i])**2 for i in range(2)]))
    else:
      sigma = 0.3 * math.sqrt(sum([(mol._atomPs[0][i]-mol._atomPs[1][i])**2 for i in range(2)]))
    sigma = round(sigma, 2)
  x, y, z = Draw.calcAtomGaussians(mol, sigma, weights=weights, step=step)
  # scaling
  if scale <= 0.0: maxScale = max(math.fabs(numpy.min(z)), math.fabs(numpy.max(z)))
  else: maxScale = scale
  # coloring
  fig.axes[0].imshow(z, cmap=colorMap, interpolation='bilinear', origin='lower', extent=(0,1,0,1), vmin=-maxScale, vmax=maxScale)
  # contour lines
  # only draw them when at least one weight is not zero
  if len([w for w in weights if w != 0.0]):
      fig.axes[0].contour(x, y, z, contourLines, colors=colors, alpha=alpha, **kwargs)
  return fig
*/






