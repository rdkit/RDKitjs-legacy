RDKitjs
=======
Port RDKit to js using emscripten

First of all we need to install emscripten and RDKit. To do this please read the following articles.

There are two possibles installation guides (for Ubuntu 14.04):  
http://baoilleach.blogspot.ch/2015/02/cheminformaticsjs-rdkit.html  

http://gmrand.blogspot.ch/2015/03/howto-install-rdkit-and-emscripten-on.html  

For OSX
http://gmrand.blogspot.ch/2015/05/howto-install-rdkit-and-emscripten-on.html

P.S.: due to a bug report we suggest to apply this patch to rdkit

Patching of RDKit files for emscripten binding
===============================================

* adding class MMFFMolProperties in rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/Builder.h
```bash  
  namespace RDKit {
    class ROMol;
    namespace MMFF {
+     class MMFFMolProperties;
```

* adding AtomTyper.h in rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/MMFF.h
```bash  
  #include <ForceField/ForceField.h>
+ #include "AtomTyper.h"
  #include "Builder.h"
```

  
Install RDKitjs
==================

git clone https://github.com/thegodone/RDKitjs.git

cd RDKitjs

* mkdir build
* cd build

### Commandline for OSX

```bash
cmake ..  \
-DCMAKE_TOOLCHAIN_FILE=/usr/local/Cellar/emscripten/1.32.0/libexec/cmake/Modules/Platform/Emscripten.cmake  \
-DRDKIT_INCLUDE_DIR=/Users/marco/Toolchain/rdkit-Release_2015_03_1/Code/  \
-DBoost_INCLUDE_DIR=/usr/local/Cellar/boost155/1.55.0_1/include/  \
-DRDKIT_LIB_DIR=/Users/marco/Toolchain/rdkit-Release_2015_03_1/lib/ \
-DEMSCRIPTEN_BIN=/usr/local/Cellar/emscripten/1.32.0/bin/
```

### Commandline for Linux

```bash
cmake .. \
-DCMAKE_TOOLCHAIN_FILE=${HOME}/Toolchain/emsdk_portable/emscripten/master/cmake/Modules/Platform/Emscripten.cmake \
-DRDKIT_INCLUDE_DIR=${HOME}/Toolchain/rdkit-Release_2015_03_1/Code \
-DBoost_INCLUDE_DIR=${HOME}/Toolchain/boost.1.57.0/include \
-DRDKIT_LIB_DIR=${HOME}/Toolchain/rdkit-Release_2015_03_1/lib \
-DEMSCRIPTEN_BIN=${HOME}/Toolchain/emsdk_portable/emscripten/master
```

* finally $ make

Then you should obtain one new file => src/rdkit.js 

Want to try it without compilation ? It's already possible
=================
Thers is a first example in the visualizer project there:  

http://www.lactame.com/visualizer/head/index.html?viewURL=http://visualizer.epfl.ch/x/urhBCrzbZ0r75WWbEmRp/view.json&dataURL=http://visualizer.epfl.ch/x/urhBCrzbZ0r75WWbEmRp/data.json  

You can draw a molecule in the botton module which will generate the 3D model using MMFF force field. 


another example of javascript can be found there:  
https://iwatobipen.wordpress.com/2015/05/21/rdkit-in-javascript/  
thanks to iwatobipen!  
source code: https://github.com/iwatobipen/rdkit_javascript  

You can also use nodejs:  
the current module is available for npm / nodejs:  
https://www.npmjs.com/package/rdkit  

Current stability status  
===============
This project is not stable but lot of basic RDKit functions are already ported look at the test/exemple.txt file for a example of function availables  

to have a complete list of available function look at the EMSCRIPTEN_BINDINGS section in rdmol.h  

there are two type of methods:   
* create a molecule (class_function)  
* apply functions on a created molecule (function)  

Future requests / help / missing RDKit functions
================
If you want to contribute or need RDKit functions not already mapped please add a comment in the issues of this project.  

Guillaume Godin  

