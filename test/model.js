
var util = require('../perf/util');
var modelsFromPDB = require('../uglymol').modelsFromPDB;

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

describe('Model', () => {
  'use strict';
  var pdb_string = util.open_as_utf8('1YJP.pdb');
  var model = modelsFromPDB(pdb_string)[0];
  it('atoms', () => {
    for (var i = 0; i < model.atoms.length; i++) {
      var atom = model.atoms[i];
      expect(atom.i_seq).toEqual(i);
    }
  });
  it('bonds', () => {
    var atoms = model.atoms;
    var simple_conn = get_connectivity_simple(atoms);
    // console.log(simple_conn);
    var model_conn = [];
    for (var i = 0; i < atoms.length; i++) {
      model_conn.push(atoms[i].bonds.sort(function (a, b) { return a - b; }));
    }
    expect(model_conn).toEqual(simple_conn);
  });
  it('next_residue', () => {
    var a1 = model.next_residue();  // first residue
    expect(a1.resseq).toEqual(1);
    expect(a1.name).toEqual('CA');
    var atom_label = a1.long_label();
    expect(atom_label.indexOf('CA /1')).toEqual(0);
    var next_res_atom = model.next_residue(a1);
    expect(next_res_atom.resseq).toEqual(2);
    expect(next_res_atom.name).toEqual('CA');
    expect(model.next_residue(next_res_atom, true)).toEqual(a1);
    var last_res_atom = model.next_residue(a1, true);
    expect(model.next_residue(last_res_atom)).toEqual(a1);
  });
  it('get_nearest_atom', () => {
    var a1 = model.next_residue();  // first residue
    var atms = [a1, model.next_residue(a1), model.next_residue(a1, true)];
    for (var i = 0; i < atms.length; i++) {
      var a = atms[i];
      var nearest = model.get_nearest_atom(a.xyz[0], a.xyz[1]+0.4, a.xyz[2]);
      expect(a).toEqual(nearest);
    }
  });
});

