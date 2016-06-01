'use strict';

var util = require('./util');
var Model = require('../src/model.js');

var suite = util.new_benchmark_suite();
var pdb_string = util.open_as_utf8('1mru.pdb');

var model;

suite.add('Model#from_pdb', function () {
  model = new Model();
  model.from_pdb(pdb_string);
});

suite.add('only calculate_connectivity', function () {
  model.calculate_connectivity();
});

suite.run();
