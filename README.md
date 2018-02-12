RDKitjs
=======

Port RDKit to js using emscripten with ASM & WASM 


Emscripten installation & activation
===============================================

Install emscripten (1.37.33 or higher).
https://kripken.github.io/emscripten-site/docs/getting_started/downloads.html


RDKit emscripten binding
===============================================

install RDKit 2015.3.1, 2016.3.5, 2016.9.4 (old version of Libraries). 
or 
install RDKit 2017.3.3, 2017.9.3 (New version of Libraries). 

you need to add eigen3 include files for 3Ddescriptors (since RDkit 2017.3)

```bash
wget https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz
wget https://github.com/rdkit/rdkit/archive/Release_2017_03_3.tar.gz
wget https://github.com/rdkit/rdkit/archive/Release_2016_09_4.tar.gz
wget https://github.com/rdkit/rdkit/archive/Release_2016_03_5.tar.gz
wget http://downloads.sourceforge.net/project/rdkit/rdkit/Q1_2015/RDKit_2015_03_1.tgz
```

Choose your version based on required RDkit functionnalities and compile RDKit first for emscripten:
```bash
tar -xf RDKit_20YY_0M_V.tgz
cd rdkit-Release_20YY_0M_V
mkdir build
cd build


cmake  .. \
-DCMAKE_TOOLCHAIN_FILE=/path/to/emsdk/emscripten/1.37.33/cmake/Modules/Platform/Emscripten.cmake \
-DRDK_BUILD_PYTHON_WRAPPERS=OFF \
-DRDK_BUILD_CPP_TESTS=OFF \
-DRDK_BUILD_SLN_SUPPORT=OFF \
-DBoost_INCLUDE_DIR=/path/to/boost/boost@1.60/1.60.0/include/ \
-DEIGEN3_INCLUDE_DIR=/path/to/eigen/3.3.4/include/eigen3/ \
-DTHREADS_PTHREAD_ARG=OFF
```

based on your number of CPU:
```bash
make -jX
```

Install RDKitjs
==================

clone this repository
```bash

git clone https://github.com/thegodone/RDKitjs.git


cd /path/to/RDKitjs/
rm -rf build
mkdir build
cd build
```

### for OSX
```bash

cmake .. \
-DCMAKE_TOOLCHAIN_FILE=/path/to/emscripten/1.37.33/cmake/Modules/Platform/Emscripten.cmake \
-DRDKIT_INCLUDE_DIR=/path/to/rdkitversion/rdkit-Release_2017_09_3/Code/ \
-DBoost_INCLUDE_DIR=/path/to/boost@1.60/1.60.0/include/  \
-DRDKIT_LIB_DIR=/path/to/rdkitversion/rdkit-Release_2017_09_3/build/lib/ \
-DRDKitNEWLIB=ON \
-DEMSCRIPTEN_BIN=/path/to/emscripten/1.37.33/
```

### for Linux
```bash
cmake .. \
-DCMAKE_TOOLCHAIN_FILE=/path/to/emscripten/1.37.33/cmake/Modules/Platform/Emscripten.cmake \
-DRDKIT_INCLUDE_DIR=/path/to/rdkitversion/rdkit-Release_2015_03_1/Code \
-DBoost_INCLUDE_DIR=/path/to/boost/include \
-DRDKIT_LIB_DIR=/path/to/rdkitversion/rdkit-Release_2015_03_1/lib \
-DEMSCRIPTEN_BIN=/path/to/emscripten/1.37.33/
```



```bash
rm src/rdkit.js

make

cp src/rdkit.js ../build/src/rdkit.js

# create a web loading file
npm run build-test

# copy web loading file into the localhost folder
cp ../test/rdkit.js /path/to/localhostcode/rdkitjs/rdkit.js

```




Please take a look at those references:

=> Ubuntu 14.04:  
http://baoilleach.blogspot.ch/2015/02/cheminformaticsjs-rdkit.html  

http://gmrand.blogspot.ch/2015/03/howto-install-rdkit-and-emscripten-on.html  

=> OSX:
http://gmrand.blogspot.ch/2015/05/howto-install-rdkit-and-emscripten-on.html



P.S.: CAUTION due to a bug report we suggest to apply this patch to rdkit

ONLY for RDKIT version 2015.03.1 : Patching of RDKit files for emscripten binding 
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

