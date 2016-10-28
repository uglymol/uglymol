
var assert = require('chai').assert;
var util = require('../perf/util');
var Model = require('../uglymol').Model;

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
  var pdb_string = util.open_as_utf8('1YJP.pdb');
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
  it('next_residue', function () {
    var a1 = model.next_residue();  // first residue
    assert.equal(a1.resseq, 1);
    assert.equal(a1.name, 'CA');
    var atom_label = a1.long_label();
    assert.equal(atom_label.indexOf('CA /1'), 0);
    var next_res_atom = model.next_residue(a1);
    assert.equal(next_res_atom.resseq, 2);
    assert.equal(next_res_atom.name, 'CA');
    assert.equal(model.next_residue(next_res_atom, true), a1);
    var last_res_atom = model.next_residue(a1, true);
    assert.equal(model.next_residue(last_res_atom), a1);
  });
  it('get_nearest_atom', function () {
    var a1 = model.next_residue();  // first residue
    var atms = [a1, model.next_residue(a1), model.next_residue(a1, true)];
    for (var i = 0; i < atms.length; i++) {
      var a = atms[i];
      var nearest = model.get_nearest_atom(a.xyz[0], a.xyz[1]+0.4, a.xyz[2]);
      assert.equal(a, nearest);
    }
  });
});

