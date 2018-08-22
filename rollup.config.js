
import buble from 'rollup-plugin-buble';
const version = require('./package.json').version

const banner = `/*!
 * UglyMol v${version}. Macromolecular Viewer for Crystallographers.
 * Copyright 2014 Nat Echols
 * Copyright 2016 Diamond Light Source Ltd
 * Copyright 2016 Marcin Wojdyr
 * Released under the MIT License.
 */`;

let build = {
  input: 'src/all.js',
  plugins: [],
  output: {
    file: 'uglymol.js',
    format: 'umd',
    name: 'UM',
    intro: `var VERSION = exports.VERSION = '${version}';\n`,
    banner,
    sourcemap: true,
    indent: false,
  },
};

if (process.env.TARGET !== 'dev') {
  // disable arrow b/c https://gitlab.com/Rich-Harris/buble/issues/158
  const transforms = { arrow: false, dangerousForOf: true };
  build.plugins.push(buble({transforms}));
}

//export default build;
export default [build];
