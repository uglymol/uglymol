
var ElMap = require('../src/elmap');
var Viewer = require('../src/viewer');
var util = require('../perf/util');

describe('Viewer', function () {
  'use strict';
  var viewer = new Viewer('viewer');
  var emap = new ElMap();
  var cmap_buf = util.open_as_array_buffer('1mru.map');
  emap.from_ccp4(cmap_buf);
  viewer.add_map(emap, false);
  viewer.toggle_cell_box();
  viewer.controls.update();
  viewer.update_camera();
  viewer.shift_clip();
  viewer.change_isolevel_by(0, 0.1);
});
