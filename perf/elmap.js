'use strict';

var util = util || require('./util'); // eslint-disable-line
var UM = UM || require('../uglymol'); // eslint-disable-line

var dsn6_buf = util.open_as_array_buffer('1mru.omap');
var map2_buf = util.open_as_array_buffer('1mru.map');
var map0_buf = util.open_as_array_buffer('1mru_m0.map');

var map;

function print_map_stats() {
  console.log('    mean/rms: ' + [map.mean, map.rms]);
}

util.bench('ElMap#from_dsn6', function () {
  map = new UM.ElMap();
  map.from_dsn6(dsn6_buf.slice(0));
},
{onComplete: print_map_stats});

util.bench('ElMap#from_ccp4 mode0', function () {
  map = new UM.ElMap();
  map.from_ccp4(map0_buf.slice(0));
});

util.bench('ElMap#from_ccp4 mode2', function () {
  map = new UM.ElMap();
  map.from_ccp4(map2_buf.slice(0));
},
{onComplete: print_map_stats});

util.bench('ElMap#extract_block', function () {
  map.extract_block(15, [25, 26, 35]);
});

