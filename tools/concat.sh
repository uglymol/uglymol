#!/bin/sh -eu

[ -e src/elmap.js ] || { echo "Run me from top-level uglymol dir"; exit 1; }

cat > uglymol.js << EOF
/* UglyMol - macromolecular viewer for crystallographers, a fork of xtal.js.
 * https://uglymol.github.io

$(cat LICENSE)
*/

var UGLYMOL_VERSION = '$(node -p "require('./package.json').version")'
var THREE = THREE || require('three'); // eslint-disable-line
EOF

cat src/unitcell.js \
    src/model.js \
    src/isosurface.js \
    src/elmap.js \
    src/lines.js \
    src/viewer.js \
  | grep -v '|| require(' \
  | grep -v 'if (typeof module !== ' \
  >> uglymol.js

cat >> uglymol.js << EOF
if (typeof module !== 'undefined') module.exports = {
  VERSION: UGLYMOL_VERSION,
  UnitCell: UnitCell,
  Model: Model,
  ElMap: ElMap,
  isosurface: isosurface,
  Viewer: Viewer
};
EOF

wc uglymol.js
