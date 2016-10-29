
//import buble from 'rollup-plugin-buble';
const version = require('./package.json').version

const banner = `/*!
 * UglyMol v${version}. Macromolecular Viewer for Crystallographers.
 * Copyright 2014 Nat Echols
 * Copyright 2016 Diamond Light Source Ltd
 * Copyright 2016 Marcin Wojdyr
 * Released under the MIT License.
 */`;

// glsl() copied on from three.js/rollup.config.js
function glsl () {
  return {
    transform ( code, id ) {
      if ( !/\.glsl$/.test( id ) ) return;

      var transformedCode = 'export default ' + JSON.stringify(
        code
          .replace( /[ \t]*\/\/.*\n/g, '' )
          .replace( /[ \t]*\/\*[\s\S]*?\*\//g, '' )
          .replace( /\n{2,}/g, '\n' )
      ) + ';';
      return {
        code: transformedCode,
        map: { mappings: '' }
      }
    }
  };
}

function three_import() {
  return {
    transform(code, id) {
      //if ( !/\.glsl$/.test( id ) ) return;
      return {
        code: code.replace("THREE from 'three'",
                           "THREE from '../tools/three-imports.js'"),
        map: { mappings: '' }
      };
    }
  };
}

let build = {
  entry: 'src/all.js',
  plugins: [/*buble() */],
  format: 'umd',
  dest: 'uglymol.js',
  moduleName: 'UM',
  external: ['three'],
  globals: { three: 'THREE' },
  intro: `exports.VERSION = '${version}';\n`,
  banner,
  sourceMap: true,
};

// build with included three.js subset: BUNDLE_DEPS=1 rollup -c
if (process.env.BUNDLE_DEPS) {
  build.plugins.push(glsl(), three_import());
  build.external = [];
  build.globals = {};
  build.dest = 'uglymol-nodeps.js';
  console.log('\nYou may run next:\n' +
              'uglifyjs uglymol-nodeps.js -cm > uglymol-nodeps.min.js\n');
}

export default build;
