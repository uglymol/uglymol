
var util = util || require('./util'); // eslint-disable-line
var UM = UM || require('../uglymol'); // eslint-disable-line

(function () {  // namespace is needed for perf.html
'use strict';

const pdb_string = util.open_as_utf8('1mru.pdb');
const cmap_buf = util.open_as_array_buffer('1mru.map');

let viewer = new UM.Viewer({});
let model = UM.modelsFromPDB(pdb_string)[0];
viewer.add_model(model);
let emap = new UM.ElMap();
emap.from_ccp4(cmap_buf);

util.bench('add trace', function () {
  viewer.config.render_style = 'trace';
  viewer.set_atomic_objects(viewer.model_bags[0]);
  viewer.clear_atomic_objects(viewer.model_bags[0]);
});

util.bench('add bonds', function () {
  viewer.config.render_style = 'lines';
  viewer.set_atomic_objects(viewer.model_bags[0]);
  viewer.clear_atomic_objects(viewer.model_bags[0]);
});

util.bench('add ribbon', function () {
  viewer.config.render_style = 'ribbon';
  viewer.set_atomic_objects(viewer.model_bags[0]);
  viewer.clear_atomic_objects(viewer.model_bags[0]);
});

util.bench('add_map+clear', function () {
  viewer.add_map(emap, false);
  viewer.clear_el_objects(viewer.map_bags.pop());
});

// makes sense only when pdb has hydrogens
//util.bench('get_visible_atoms', function () {
//  viewer.model_bags[0].get_visible_atoms();
//});
})();
