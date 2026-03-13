
var GM = require('../gemmimol');
var util = require('../perf/util');

describe('Viewer', () => {
  'use strict';
  var viewer = new GM.Viewer('viewer');
  var emap = new GM.ElMap();
  var cmap_buf = util.open_as_array_buffer('1mru.map');
  var gemmi;
  var model;
  var model2;
  beforeAll(function () {
    return util.load_gemmi().then(function (loaded) {
      gemmi = loaded;
      emap.from_ccp4(cmap_buf, true, gemmi);
      return Promise.all([
        util.load_models_from_gemmi('1mru.pdb'),
        util.load_models_from_gemmi('1yk4.pdb'),
      ]);
    }).then(function (models) {
      model = models[0][0];
      model2 = models[1][0];
    });
  });
  it('misc calls (1mru)', () => {
    viewer.add_map(emap, false);
    viewer.toggle_map_visibility(viewer.map_bags[0], false);
    viewer.toggle_map_visibility(viewer.map_bags[0], true);
    viewer.toggle_cell_box();
    viewer.controls.update();
    viewer.update_camera();
    viewer.shift_clip();
    viewer.change_isolevel_by(0, 0.1);
    viewer.center_next_residue();
    viewer.add_model(model);
    viewer.center_next_residue();
    viewer.recenter();
    viewer.recenter([11, 22, 33]);
    viewer.select_atom({bag: viewer.model_bags[0], atom: model.atoms[1]});
  });

  it('misc calls (1yk4)', () => {
    viewer.add_model(model2);
    viewer.config.hydrogens = true;
    viewer.recenter();
  });

  it('keydown', () => {
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
