-------- Update 15/11/2018
Dear All,

Greg and I, consider this project is important for the community.
We will transfer it to RDKit repo in coming months.

---------

# RDKitjs

Port of RDKit to JavaScript using emscripten and WebAssembly.
> Note, this has been added to the RDKit repo. Please see the blog post [here](http://rdkit.blogspot.com/2019/11/introducing-new-rdkit-javascript.html), which has code [here](https://github.com/rdkit/rdkit/tree/Release_2019_09_1/Code/MinimalLib/demo). 
# Installation

```bash
npm install rdkit
```

# Usage

Using NodeJS

```js
// applications.js
const RDKit = require('RDKit');
await RDKit.load();

var my_mol = "COc1ccc(CCN(C(=O)CCCBr)C2c3cc(NC(C)=O)c([N+](=O)[O-])cc3OC(C)(C)C2O)cc1";

function atoms(smi) {
  var mol = RDKit.Mol.fromSmiles(smi);
  var num_atoms = mol.getNumAtoms();
  console.log( num_atoms );
}

atoms(my_mol)

// returns 37
```

# Development

## Prerequisites

For now, this library can only be compiled from a Unix system (no Windows).
You need to have on your system:

* Node.js 8 or later
* cmake
* make

## Emscripten installation & activation

Install emscripten using those instructions: https://kripken.github.io/emscripten-site/docs/getting_started/downloads.html#linux-and-mac-os-x  
Do not forget to activate it with `./emsdk activate latest`.

## Installation of dependencies

```bash
npm install
npm run install-deps
```

This will install all required dependencies (including RDKit itself) and compile
the RDKit library files.

## Compilation of WebAssembly module

```bash
npm run build
```

## Execute tests

```bash
# compiles the project and executes tests
npm run test

# executes tests without compiling (faster if compilation was already done)
npm run test-only
```

# Want to try it without compilation ? It's already possible

There is the first example in the visualizer project there:

https://www.cheminfo.org/Chemistry/Cheminformatics/RDKit_demo/index.html

You can draw a molecule in the bottom module which will generate the 3D model using MMFF force field.

another example of javascript can be found there:  
https://iwatobipen.wordpress.com/2015/05/21/rdkit-in-javascript/  
thanks to iwatobipen!  
source code: https://github.com/iwatobipen/rdkit_javascript

Try RDKitjs interactively with [RunKit](https://npm.runkit.com/rdkit)

You can also use nodejs:  
the current module is available for npm / nodejs:  
https://www.npmjs.com/package/rdkit

# Current stability status

This project is not stable and still under heavy development.

# Feature requests / help / missing RDKit functions

If you want to contribute or need RDKit functions not already mapped please add
a comment in the issues of this project.
