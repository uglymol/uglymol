'use strict';

var ElMap = require('../src/elmap');
var isosurface = require('../src/isosurface');
var util = require('./util');

var dsn6_file = '1mru.omap';
var dmap_buf = util.open_as_array_buffer(dsn6_file);

var suite = util.new_benchmark_suite();
var map = new ElMap();
map.from_dsn6(dmap_buf.slice(0));
map.extract_block(15, [25, 26, 35]);
var geometry;

suite.add('isosurface', function () {
  var bl = map.block;
  geometry = isosurface(bl.points, bl.values, bl.size, 1.0);
});

suite.run();
