'use strict';

var util = util || require('./util'); // eslint-disable-line
var GM = GM || require('../gemmimol'); // eslint-disable-line

let model;

var setup = Promise.all([util.load_gemmi(), util.load_models_from_gemmi('1mru.pdb')])
  .then(function (result) {
    var gemmi = result[0];
    var models = result[1];
    var buffer = util.open_as_array_buffer('1mru.pdb');
    model = models[0];

    util.bench('Gemmi read_structure', function () {
      var st = gemmi.read_structure(buffer, '1mru.pdb');
      st.delete();
    });

    util.bench('only calculate_connectivity', function () {
      model.calculate_connectivity();
    });
  });

if (typeof module !== 'undefined') {
  module.exports = setup;
}
