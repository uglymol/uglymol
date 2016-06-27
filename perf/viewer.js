'use strict';

var util = util || require('./util'); // eslint-disable-line
var Model = Model || require('../src/model.js'); // eslint-disable-line
var Viewer = Viewer || require('../src/viewer.js'); // eslint-disable-line

(function () {
  var pdb_string = util.open_as_utf8('1mru.pdb');

  var viewer = new Viewer();
  var model = new Model();
  model.from_pdb(pdb_string);
  viewer.set_model(model);

  util.bench('make_trace', function () {
    viewer.config.render_style = 'trace';
    viewer.set_atomic_objects(viewer.model_bags[0]);
    viewer.clear_atomic_objects(viewer.model_bags[0]);
  });

  util.bench('make_bonds', function () {
    viewer.config.render_style = 'lines';
    viewer.set_atomic_objects(viewer.model_bags[0]);
    viewer.clear_atomic_objects(viewer.model_bags[0]);
  });
})();
