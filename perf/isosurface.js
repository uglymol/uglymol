
var util = util || require('./util'); // eslint-disable-line
var UM = UM || require('../uglymol'); // eslint-disable-line

(function () {  // namespace is needed for perf.html
'use strict';

var dmap_buf = util.open_as_array_buffer('1mru.omap');
var map = new UM.ElMap();
map.from_dsn6(dmap_buf.slice(0));
map.extract_block(15, [25, 26, 35]);

util.bench('isosurface', function () {
  map.isomesh_in_block(1.5);
});
})();
