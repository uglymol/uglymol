
var ElMap = require('../src/elmap');
var Model = require('../src/model');
var Viewer = require('../src/viewer');
var util = require('../perf/util');

describe('Viewer', function () {
  'use strict';
  var viewer = new Viewer('viewer');
  var emap = new ElMap();
  var cmap_buf = util.open_as_array_buffer('1mru.map');
  emap.from_ccp4(cmap_buf);
  var pdb_string = util.open_as_utf8('1mru.pdb');
  var model = new Model();
  model.from_pdb(pdb_string);
  it('misc calls', function () {
    viewer.add_map(emap, false);
    viewer.toggle_map_visibility(viewer.map_bags[0], false);
    viewer.toggle_map_visibility(viewer.map_bags[0], true);
    viewer.toggle_cell_box();
    viewer.controls.update();
    viewer.update_camera();
    viewer.shift_clip();
    viewer.change_isolevel_by(0, 0.1);
    viewer.center_next_residue();
    viewer.set_model(model);
    viewer.center_next_residue();
    viewer.set_selection(model.atoms[0]);
  });

  it('keydown', function () {
    for (var i = 0; i < 5; i++) { // try all color schemes
      viewer.keydown({keyCode: 67, shiftKey: false}); // c
    }
    viewer.keydown({keyCode: 67, shiftKey: true}); // C
    //          d   f   m   n   [    ]
    var keys = [68, 70, 77, 78, 219, 221];
    for (var j = 0; j < keys.length; j++) {
      viewer.keydown({keyCode: keys[j]});
    }
  });
});
