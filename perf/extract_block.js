'use strict';

var util = require('./util');
var ElMap = require('../src/elmap');

var cmap_buf = util.open_as_array_buffer('1mru_2mFo-DFc.ccp4');

var suite = util.new_benchmark_suite();
var map = new ElMap();
map.from_ccp4(cmap_buf.slice(0));

suite.add('ElMap#extract_block', function () {
  map.extract_block(15, [25, 26, 35]);
});

suite.run();
