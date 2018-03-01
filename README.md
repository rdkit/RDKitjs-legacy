# RDKitjs

Port of RDKit to JavaScript using emscripten and WebAssembly.

# Installation

```bash
npm install rdkit
```

# Usage

TODO

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

Thers is a first example in the visualizer project there:

https://www.cheminfo.org/Chemistry/Cheminformatics/RDKit_demo/index.html

You can draw a molecule in the botton module which will generate the 3D model using MMFF force field.

another example of javascript can be found there:  
https://iwatobipen.wordpress.com/2015/05/21/rdkit-in-javascript/  
thanks to iwatobipen!  
source code: https://github.com/iwatobipen/rdkit_javascript

You can also use nodejs:  
the current module is available for npm / nodejs:  
https://www.npmjs.com/package/rdkit

# Current stability status

This project is not stable and still under heavy development.

# Feature requests / help / missing RDKit functions

If you want to contribute or need RDKit functions not already mapped please add
a comment in the issues of this project.
