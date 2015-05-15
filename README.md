# RDKitjs
port RDKit to js using emscritpen

there are two possibles installation guides (for Ubuntu 14.04):
http://baoilleach.blogspot.ch/2015/02/cheminformaticsjs-rdkit.html

http://gmrand.blogspot.ch/2015/03/howto-install-rdkit-and-emscripten-on.html

### Install emscripten
http://kripken.github.io/emscripten-site/docs/getting_started/downloads.html


I. Setup your development environmemt (for mac os x):
Required for compilation of emscripten
emscripten (master or skd)
clang (of emscripten same version than emscripten)
node 
python 2.7.X (for mac os x 2.7.9 is available).
boost_1_57_0
git
cmake (a least 3.0)
rdkit Release 2014.09.03 or 2015.03.01 (both working!).

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
$ make -j5 (but if you have only 2 cores use -j2)
$ make install

II. link rdmol folder with the emscripten/rdkit/boost folders
the best solution is to create 3 symbolic links into rdmol folder  
$ ln -fs /path/to/build/boost/boost.1.57/include/ path/to/rdmol/include

ln -fs /emscripten/rdkit-Release_2015_03_1/code/ /emscripten/rdmol/build/code

$ ln -fs /path/to/build/rdkit-Release_2014_09_2/code/ path/to/rdmol/code
ln -fs /emscripten/rdkit-Release_2015_03_1/code/ /emscripten/rdmol/build/code
$ ln -fs /path/to/build/rdkit-Release_2014_09_2/lib/ path/to/rdmol/lib
ln -fs /emscripten/rdkit-Release_2015_03_1/lib/ /emscripten/rdmol/build/lib

to recompile your own rdmol.cpp:
$ cd /path/to/rdmol
$ mkdir build
$ cd build
$ path/to/emscripten/em++  --bind -o rdmol.js ../rdmol.cpp -Icode -Iinclude lib/libGraphMol.so lib/libDescriptors.so lib/libRDGeneral.so lib/libRDGeometryLib.so lib/libSmilesParse.so lib/libDataStructs.so lib/libFingerprints.so lib/libSubgraphs.so  -O2 

you should obtain 2 files => rdmol.js & rdmol.js.mem















