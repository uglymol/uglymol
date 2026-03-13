'use strict';

var util = util || require('./util');
var GM = GM || require('../gemmimol');

const dsn6_buf = util.open_as_array_buffer('1mru.omap');
const map2_buf = util.open_as_array_buffer('1mru.map');
const map0_buf = util.open_as_array_buffer('1mru_m0.map');

let map;

function print_map_stats() {
  console.log('    mean/rms: ' + [map.stats.mean, map.stats.rms]);
}

var setup = util.load_gemmi().then(function (gemmi) {
  util.bench('ElMap#from_dsn6', function () {
    map = new GM.ElMap();
    map.from_dsn6(dsn6_buf.slice(0));
  }, {onComplete: print_map_stats});

  util.bench('ElMap#from_ccp4 mode0', function () {
    map = new GM.ElMap();
    map.from_ccp4(map0_buf.slice(0), true, gemmi);
  });

  util.bench('ElMap#from_ccp4 mode2', function () {
    map = new GM.ElMap();
    map.from_ccp4(map2_buf.slice(0), true, gemmi);
  }, {onComplete: print_map_stats});

  util.bench('ElMap#extract_block', function () {
    map.extract_block(15, [25, 26, 35]);
  });
});

if (typeof module !== 'undefined') {
  module.exports = setup;
}
