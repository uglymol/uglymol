'use strict';

var util = util || require('./util'); // eslint-disable-line
var Model = Model || require('../src/model'); // eslint-disable-line
var ElMap = ElMap || require('../src/elmap'); // eslint-disable-line
var Viewer = Viewer || require('../src/viewer'); // eslint-disable-line

(function () {
  var pdb_string = util.open_as_utf8('1mru.pdb');
  var cmap_buf = util.open_as_array_buffer('1mru.map');

  var viewer = new Viewer();
  var model = new Model();
  model.from_pdb(pdb_string);
  viewer.set_model(model);
  var emap = new ElMap();
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
})();
