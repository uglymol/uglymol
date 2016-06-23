'use strict';

var util = util || require('./util'); // eslint-disable-line
var ElMap = ElMap || require('../src/elmap'); // eslint-disable-line

var dmap_buf = util.open_as_array_buffer('1mru.omap');
var cmap_buf = util.open_as_array_buffer('1mru_2mFo-DFc.ccp4');

var map;

function print_map_stats() {
  console.log('mean/rms: ' + [map.mean, map.rms]);
}

util.bench('ElMap#from_dsn6', function () {
  map = new ElMap();
  map.from_dsn6(dmap_buf.slice(0));
},
{onComplete: print_map_stats});

util.bench('ElMap#from_ccp4', function () {
  map = new ElMap();
  map.from_ccp4(cmap_buf.slice(0));
},
{onComplete: print_map_stats});

util.bench('ElMap#extract_block', function () {
  map.extract_block(15, [25, 26, 35]);
});

