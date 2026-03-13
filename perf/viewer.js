
var util = util || require('./util'); // eslint-disable-line
var GM = GM || require('../gemmimol'); // eslint-disable-line

(function () {  // namespace is needed for perf.html
'use strict';

const cmap_buf = util.open_as_array_buffer('1mru.map');

let viewer = new GM.Viewer({});
let emap = new GM.ElMap();

var setup = Promise.all([util.load_gemmi(), util.load_models_from_gemmi('1mru.pdb')])
  .then(function (result) {
    var gemmi = result[0];
    var models = result[1];
    emap.from_ccp4(cmap_buf, true, gemmi);
    let model = models[0];
    viewer.add_model(model);

    util.bench('add trace', function () {
      viewer.config.render_style = 'trace';
      viewer.set_model_objects(viewer.model_bags[0]);
      viewer.clear_model_objects(viewer.model_bags[0]);
    });

    util.bench('add bonds', function () {
      viewer.config.render_style = 'lines';
      viewer.set_model_objects(viewer.model_bags[0]);
      viewer.clear_model_objects(viewer.model_bags[0]);
    });

    util.bench('add ribbon', function () {
      viewer.config.render_style = 'ribbon';
      viewer.set_model_objects(viewer.model_bags[0]);
      viewer.clear_model_objects(viewer.model_bags[0]);
    });

    util.bench('add_map+clear', function () {
      viewer.add_map(emap, false);
      viewer.clear_el_objects(viewer.map_bags.pop());
    });
  });

// makes sense only when pdb has hydrogens
//util.bench('get_visible_atoms', function () {
//  viewer.model_bags[0].get_visible_atoms();
//});

if (typeof module !== 'undefined') {
  module.exports = setup;
}
})();
