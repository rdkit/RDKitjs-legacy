# RDKitjs
port RDKit to js using emscritpen

there are two possibles installation guides (for Ubuntu 14.04):  
http://baoilleach.blogspot.ch/2015/02/cheminformaticsjs-rdkit.html  

http://gmrand.blogspot.ch/2015/03/howto-install-rdkit-and-emscripten-on.html  

### Install emscripten
http://kripken.github.io/emscripten-site/docs/getting_started/downloads.html


I. Setup your development environmemt (for mac os x):  
Required for compilation of emscripten:  
emscripten (manually or emsdk)  
clang (of emscripten same version than emscripten)  
node   
python 2.7.X (for mac os x 2.7.9 is available).  
boost_1_57_0 
java jre 
git  
cmake (at least 3.0)  
rdkit Release 2014.09.03 or 2015.03.01 (both working! for new depiction use 2015).

crutial part: you need to compile RDKit & boost with clang/emscripten:

open a terminal (with sudo privileges):  $ or # 
for boost_1_57_0: 
Download the boost and umpact than go to the boost_1_57_0 folder
$ cd boost_1_57_0


$ ./bootstrap.sh --prefix=/path/to/build/boost/boost.1.57  
$ ./b2 --with=all -j5 install  

for rdkit:  
example   
$ wget http://downloads.sourceforge.net/project/rdkit/rdkit/Q3_2014/RDKit_2014_09_2.tar.gz  
$ tar -xf RDKit_2014_09_2.tar.gz  
$ cd rdkit-Release_2014_09_2/  
$ mkdir build  
$ cd build  

Edit the CMakeLists.txt of RDKit and comment the line 41-42:  
   #include(TestBigEndian)  
   #TEST_BIG_ENDIAN(RDK_BIG_ENDIAN)  


you need to compile RDKit with special flags:
cmake  -DCMAKE_TOOLCHAIN_FILE=/path/to/emscripten/emscripten/cmake/Modules/Platform/Emscripten.cmake -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_CPP_TESTS=OFF -DRDK_BUILD_SLN_SUPPORT=OFF -DBoost_INCLUDE_DIR=/path/to/build/boost/boost.1.57/include/  -DTHREADS_PTHREAD_ARG=OFF ..  

then execute the folowing commands (the first one takes 5-10 min)  
$ make -j5 (of course if you have only 2 cores replace -j5 by -j2)  
$ make install  


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



After the patch you can compile the code:  
to recompile your own rdkit.cpp:  
$ cd /path/to/rdkit  
$ sh build.sh (don't forget to rename your paths to the rdkit & boost/include & path/to/build/boost/boost.1.57/include/ 
/path/to/build/rdkit-Release_2014_09_2/
/path/to/emscripten/em++ in the build.sh -- remark build_fast.sh use -O1 faster)

you should obtain one new file => rdkit.js 















