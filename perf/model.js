'use strict';

var util = util || require('./util'); // eslint-disable-line
var UM = UM || require('../uglymol'); // eslint-disable-line

const pdb_string = util.open_as_utf8('1mru.pdb');
let model;

util.bench('Model#from_pdb', function () {
  model = UM.modelsFromPDB(pdb_string)[0];
});

util.bench('only calculate_connectivity', function () {
  model.calculate_connectivity();
});

