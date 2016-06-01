'use strict';

var util = require('./util');
var Model = require('../src/model.js');
var Viewer = require('../src/viewer.js');

var suite = util.new_benchmark_suite();
var pdb_string = util.open_as_utf8('1mru.pdb');

var viewer = new Viewer();
var model = new Model();
model.from_pdb(pdb_string);
viewer.add_model_bag(model, '1mru');

suite.add('make_trace', function () {
  viewer.model_config.render_style = 'trace';
  viewer.set_atomic_objects(viewer.models[0]);
});

suite.add('make_bonds', function () {
  viewer.model_config.render_style = 'lines';
  viewer.set_atomic_objects(viewer.models[0]);
});

suite.run();
