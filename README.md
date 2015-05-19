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

%%%% pacth of RDKit files for emscripten
*** ./Code/GraphMol/ForceFieldHelpers/MMFF/Builder.h       Wed May 13 09:32:12 2015
  namespace RDKit {
    class ROMol;
    namespace MMFF {
+     class MMFFMolProperties;
  
*** /build/common/rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/MMFF.h        Wed May 13 10:09:29 2015
  #include <ForceField/ForceField.h>
+ #include "AtomTyper.h"
  #include "Builder.h"

  
Install RDKitjs
==================

git clone https://github.com/thegodone/RDKitjs.git

cd RDKitjs

* mkdir build
* cd build

### Commandline for OSX

```bash
cmake .. \
-DCMAKE_TOOLCHAIN_FILE=/usr/local/Cellar/emscripten/1.32.0/libexec/cmake/Modules/Platform/Emscripten.cmake  \
-DRDKIT_INCLUDE_DIR=/Users/marco/Toolchain/rdkit-Release_2015_03_1/Code/ \
-DBoost_INCLUDE_DIR=/usr/local/Cellar/boost155/1.55.0_1/include/ \
-DRDKIT_LIB_DIR=/Users/marco/Toolchain/rdkit-Release_2015_03_1/lib/
```

### Commandline for Linux

```bash
cmake .. \
-DCMAKE_TOOLCHAIN_FILE=${HOME}/Toolchain/emsdk_portable/emscripten/master/cmake/Modules/Platform/Emscripten.cmake \
-DRDKIT_INCLUDE_DIR=${HOME}/Toolchain/rdkit-Release_2014_09_2/Code \
-DBoost_INCLUDE_DIR=${HOME}/Toolchain/boost.1.57.0/include/ \
-DRDKIT_LIB_DIR=${HOME}/Toolchain/rdkit-Release_2015_03_1/lib/
-DEMSCRIPTEN_BIN=${HOME}/Toolchain/emsdk_portable/emscripten/master

make
```


Then you should obtain one new file => src/rdkit.js 

the current module is available for npm/nodje:
https://www.npmjs.com/package/rdkit


