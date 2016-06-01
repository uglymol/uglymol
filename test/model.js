
var assert = require('chai').assert;
var fs = require('fs');
var Model = require('../src/model');

// O(n^2) loop, for testing purposes only
function get_connectivity_simple(atoms) {
  'use strict';
  var connectivity = [];
  var i;
  for (i = 0; i < atoms.length; i++) {
    connectivity.push([]);
  }
  for (i = 0; i < atoms.length; i++) {
    for (var j = i + 1; j < atoms.length; j++) {
      if (atoms[i].is_bonded_to(atoms[j])) {
        connectivity[i].push(j);
        connectivity[j].push(i);
      }
    }
  }
  return connectivity;
}

describe('Model', function () {
  'use strict';
  var model = new Model();
  //var fn = __dirname + '/../data/1YJP.pdb';
  var fn = __dirname + '/../data/2gkg.pdb';
  var pdb_string = fs.readFileSync(fn, {encoding: 'utf8'});
  model.from_pdb(pdb_string);
  it('atoms', function () {
    for (var i = 0; i < model.atoms.length; i++) {
      var atom = model.atoms[i];
      assert.equal(atom.i_seq, i);
    }
  });
  it('bonds', function () {
    var atoms = model.atoms;
    var simple_conn = get_connectivity_simple(atoms);
    // console.log(simple_conn);
    var model_conn = [];
    for (var i = 0; i < atoms.length; i++) {
      model_conn.push(atoms[i].bonds.sort(function (a, b) { return a - b; }));
    }
    assert.deepEqual(model_conn, simple_conn);
  });
});

