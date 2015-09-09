def bivariate_normal(X, Y, sigmax=1.0, sigmay=1.0, mux=0.0, muy=0.0, sigmaxy=0.0):
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

def calcAtomGaussians(mol,a=0.03,step=0.02,weights=None):
  import numpy
  from matplotlib import mlab
  x = numpy.arange(0,1,step)
  y = numpy.arange(0,1,step)
  X,Y = numpy.meshgrid(x,y)
  if weights is None:
    weights=[1.]*mol.GetNumAtoms()
  Z = mlab.bivariate_normal(X,Y,a,a,mol._atomPs[0][0], mol._atomPs[0][1])*weights[0] # this is not bivariate case ... only univariate no mixtures
  for i in range(1,mol.GetNumAtoms()):
    Zp = mlab.bivariate_normal(X,Y,a,a,mol._atomPs[i][0], mol._atomPs[i][1])
    Z += Zp*weights[i]
  return X,Y,Z
  

#matplotlib.mlab.bivariate_normal(X, Y, sigmax=1.0, sigmay=1.0, mux=0.0, muy=0.0, sigmaxy=0.0)


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
