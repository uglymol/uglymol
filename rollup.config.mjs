
import buble from '@rollup/plugin-buble';
import sucrase from '@rollup/plugin-sucrase';
import { createRequire } from 'node:module';
const require = createRequire(import.meta.url);
const version = require('./package.json').version

const banner = `/*!
 * UglyMol v${version}. Macromolecular Viewer for Crystallographers.
 * Copyright 2014 Nat Echols
 * Copyright 2016 Diamond Light Source Ltd
 * Copyright 2016 Marcin Wojdyr
 * Released under the MIT License.
 */`;

let build = {
  input: 'src/all.ts',
  //treeshake: false,
  plugins: [
    sucrase({
      include: ['**/src/*.ts'],
      transforms: ['typescript']
    })
  ],
  output: {
    file: 'uglymol.js',
    format: 'umd',
    name: 'UM',
    intro: `var VERSION = exports.VERSION = '${version}';\n`,
    banner,
    sourcemap: false,
    indent: false,
  },
};

if (process.env.TARGET !== 'dev') {
  // disable arrow b/c https://github.com/bublejs/buble/issues/208
  const transforms = { arrow: false, dangerousForOf: true };
  build.plugins.push(buble({transforms}));
}

//export default build;
export default [build];
