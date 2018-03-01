'use strict';

const deps = [
  {
    name: 'rdkit',
    version: '2017_09_3',
    fileCheck: 'license.txt',
    url: 'https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz'
  },
  {
    name: 'boost',
    version: '1.66.0',
    fileCheck: 'boost/version.hpp',
    url:
      'https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz'
  },
  {
    name: 'eigen',
    version: '3.3.4',
    fileCheck: 'README.md',
    url: 'http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz'
  }
];

module.exports = deps;
