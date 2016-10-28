'use strict';

var util = util || require('./util'); // eslint-disable-line
var UM = UM || require('../uglymol'); // eslint-disable-line

var pdb_string = util.open_as_utf8('1mru.pdb');
var model;

util.bench('Model#from_pdb', function () {
  model = new UM.Model();
  model.from_pdb(pdb_string);
});

util.bench('only calculate_connectivity', function () {
  model.calculate_connectivity();
});

