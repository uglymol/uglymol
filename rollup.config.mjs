
import sucrase from '@rollup/plugin-sucrase';
import terser from '@rollup/plugin-terser';
import { createRequire } from 'node:module';
const require = createRequire(import.meta.url);
const version = require('./package.json').version;

const banner = `/*!
 * GemmiMol v${version}. Macromolecular Viewer for Crystallographers.
 * Copyright 2014 Nat Echols
 * Copyright 2016 Diamond Light Source Ltd
 * Copyright 2016 Marcin Wojdyr
 * Released under the MIT License.
 */`;

const basePlugins = [
  sucrase({
    include: ['**/src/*.ts'],
    transforms: ['typescript'],
  }),
];

const output = {
  format: 'umd',
  name: 'GM',
  intro: `var VERSION = exports.VERSION = '${version}';\n`,
  banner,
  sourcemap: false,
  indent: false,
};

const build = {
  input: 'src/all.ts',
  plugins: basePlugins,
  output: {...output, file: 'gemmimol.js'},
};

const minified = {
  input: 'src/all.ts',
  plugins: [
    ...basePlugins,
    terser({
      compress: true,
      mangle: true,
      format: {
        comments: /^!/,
      },
    }),
  ],
  output: {...output, file: 'gemmimol.min.js'},
};

export default [build, minified];
