
var UM = require('../uglymol');
var util = require('../perf/util');

describe('Viewer', function () {
  'use strict';
  var viewer = new UM.Viewer('viewer');
  var emap = new UM.ElMap();
  var cmap_buf = util.open_as_array_buffer('1mru.map');
  emap.from_ccp4(cmap_buf);
  var pdb_string = util.open_as_utf8('1mru.pdb');
  var model = new UM.Model();
  model.from_pdb(pdb_string);
  it('misc calls (1mru)', function () {
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
    viewer.recenter();
    viewer.recenter([11, 22, 33]);
    viewer.select_atom(model.atoms[1]);
  });

  pdb_string = util.open_as_utf8('1yk4.pdb');
  model.from_pdb(pdb_string);
  it('misc calls (1yk4)', function () {
    viewer.set_model(model);
    viewer.config.hydrogens = true;
    viewer.recenter();
  });

  it('keydown', function () {
    function press(codes) {
      for (var i = 0; i < codes.length; i++) {
        var code = codes[i];
        var shift = false;
        if (typeof code === 'string') {
          shift = (code !== code.toLowerCase());
          code = code.toUpperCase().charCodeAt(0);
        }
        viewer.keydown({keyCode: code, shiftKey: shift});
      }
    }
    press(['c', 'b', 'C', 'b', 'b', 'b', 'c', 'c', 'c', 'c']); // colors
    press(['d', 'f', 'm', 'n', 219/*[*/, 221/*]*/]);
    press(['y', 'w', 'u', 220/*backslash*/, 'y', 'w', 'u', 220, 'r']);
    press([99/*numpad 3*/, 110/*decimal point (Mac)*/]);
    press([107/*add*/, 109/*subtract*/]);
    press([32/*space*/, 999/*dummy*/]);
    press(['t', 't', 'q', 'q']);
  });
});
