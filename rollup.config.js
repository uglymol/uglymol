
//import buble from 'rollup-plugin-buble';
const version = require('./package.json').version

const banner = `/*!
 * UglyMol v${version}. Macromolecular Viewer for Crystallographers.
 * Copyright 2014 Nat Echols
 * Copyright 2016 Diamond Light Source Ltd
 * Copyright 2016 Marcin Wojdyr
 * Released under the MIT License.
 */`;

export default {
  entry: 'src/all.js',
  plugins: [/*buble() */],
  format: 'umd',
  dest: 'uglymol.js',
  moduleName: 'UM',
  external: ['three'],
  globals: { three: 'THREE' },
  intro: `exports.VERSION = '${version}';\n`,
  banner,
};
