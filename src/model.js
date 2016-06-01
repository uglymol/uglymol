
var UnitCell = UnitCell || require('./unitcell'); // eslint-disable-line

var Model = (function () {
'use strict';

var AMINO_ACIDS = [
  'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
  'LYS', 'MET', 'MSE', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'UNK'
];
var NUCLEIC_ACIDS = [
  'DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U', 'rA', 'rC', 'rG', 'rU',
  'Ar', 'Cr', 'Gr', 'Ur'
];

var NOT_LIGANDS = ['HOH'].concat(AMINO_ACIDS, NUCLEIC_ACIDS);

function Model() {
  this.atoms = [];
  this.unit_cell = null;
  this.space_group = null;
  this.has_hydrogens = false;
  this.lower_bound = null;
  this.upper_bound = null;
}

Model.prototype.from_pdb = function (pdb_string) {
  var lines = pdb_string.split('\n');
  var chain_index = 0;  // will be ++'ed for the first atom
  var last_chain = null;
  var atom_i_seq = 0;
  //var last_atom = null;
  for (var i = 0; i < lines.length; i++) {
    var line = lines[i];
    var rec_type = line.substring(0, 6);
    if (rec_type === 'ATOM  ' || rec_type === 'HETATM') {
      var new_atom = new Atom();
      new_atom.from_pdb_line(line);
      new_atom.i_seq = atom_i_seq++;
      if (!this.has_hydrogens && new_atom.element === 'H') {
        this.has_hydrogens = true;
      }
      if (new_atom.chain !== last_chain) {
        chain_index++;
      }
      new_atom.chain_index = chain_index;
      last_chain = new_atom.chain;
      //last_atom = new_atom;
      this.atoms.push(new_atom);
    } else if (rec_type === 'ANISOU') {
      //last_atom.set_uij_from_anisou(line);
    } else if (rec_type === 'CRYST1') {
      var a = parseFloat(line.substring(6, 15));
      var b = parseFloat(line.substring(15, 24));
      var c = parseFloat(line.substring(24, 33));
      var alpha = parseFloat(line.substring(33, 40));
      var beta = parseFloat(line.substring(40, 47));
      var gamma = parseFloat(line.substring(47, 54));
      //var sg_symbol = line.substring(55, 66);
      this.unit_cell = new UnitCell(a, b, c, alpha, beta, gamma);
    } else if (rec_type.substring(0, 3) === 'TER') {
      last_chain = null;
    }
  }
  if (this.atoms.length === 0) {
    throw Error('No atom records found.');
  }
  this.calculate_bounds();
  this.calculate_connectivity();
};

Model.prototype.calculate_bounds = function () {
  var lower = this.lower_bound = [Infinity, Infinity, Infinity];
  var upper = this.upper_bound = [-Infinity, -Infinity, -Infinity];
  for (var i = 0; i < this.atoms.length; i++) {
    var atom = this.atoms[i];
    for (var j = 0; j < 3; j++) {
      var v = atom.xyz[j];
      if (v < lower[j]) { lower[j] = v; }
      if (v > upper[j]) { upper[j] = v; }
    }
  }
  // with a margin
  for (var k = 0; k < 3; ++k) {
    lower[k] -= 0.001;
    upper[k] += 0.001;
  }
};

Model.prototype.next_residue = function (atom, backward) {
  var len = this.atoms.length;
  var start = (atom ? atom.i_seq : 0) + len;  // +len to avoid idx<0 below
  for (var i = (atom ? 1 : 0); i < len; i++) {
    var idx = (start + (backward ? -i : i)) % len;
    var a = this.atoms[idx];
    if (!a.is_main_conformer()) { continue; }
    if ((a.name === 'CA' && a.element === 'C') || a.name === 'P') {
      return a;
    }
  }
};

Model.prototype.extract_trace = function () {
  var segments = [];
  var current_segment = [];
  var last_atom = null;
  var i;
  for (i = 0; i < this.atoms.length; i++) {
    var atom = this.atoms[i];
    if (atom.altloc !== '' && atom.altloc !== 'A') { continue; }
    if ((atom.name === 'CA' && atom.element === 'C') || atom.name === 'P') {
      var start_new = true;
      if (last_atom !== null && last_atom.chain_index === atom.chain_index) {
        var dxyz2 = atom.distance_sq(last_atom);
        if ((atom.name === 'CA' && dxyz2 <= 5.5*5.5) ||
            (atom.name === 'P' && dxyz2 < 7.5*7.5)) {
          current_segment.push(atom);
          start_new = false;
        }
      }
      if (start_new) {
        if (current_segment.length > 2) {
          segments.push(current_segment);
        }
        current_segment = [atom];
      }
      last_atom = atom;
    }
  }
  if (current_segment.length > 2) {
    segments.push(current_segment);
  }
  //console.log(segments.length + " segments extracted");
  return segments;
};

Model.prototype.get_center = function () {
  var xsum = 0, ysum = 0, zsum = 0;  // eslint-disable-line
  var n_atoms = this.atoms.length;
  for (var i = 0; i < n_atoms; i++) {
    var xyz = this.atoms[i].xyz;
    xsum += xyz[0];
    ysum += xyz[1];
    zsum += xyz[2];
  }
  return [xsum / n_atoms, ysum / n_atoms, zsum / n_atoms];
};


// Single atom and associated labels
function Atom() {
  this.hetero = false;
  this.name = '';
  this.altloc = '';
  this.resname = '';
  this.chain = '';
  this.chain_index = null;
  this.resseq = '';
  this.icode = null;
  this.xyz = [0, 0, 0];
  this.occ = 1.0;
  this.b = 0;
  this.element = '';
  this.charge = 0;
  this.i_seq = null;
  this.is_ligand = null;
  this.bonds = [];
}

// http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
Atom.prototype.from_pdb_line = function (pdb_line) {
  if (pdb_line.length < 66) {
    throw Error('ATOM or HETATM record is too short: ' + pdb_line);
  }
  var rec_type = pdb_line.substring(0, 6);
  if (rec_type === 'HETATM') {
    this.hetero = true;
  } else if (rec_type !== 'ATOM  ') {
    throw Error('Wrong record type: ' + rec_type);
  }
  this.name = pdb_line.substring(12, 16).trim();
  this.altloc = pdb_line.substring(16, 17).trim();
  this.resname = pdb_line.substring(17, 20).trim();
  this.chain = pdb_line.substring(20, 22).trim();
  this.resseq = pdb_line.substring(22, 26);
  this.icode = pdb_line.substring(26, 27).trim();
  var x = parseFloat(pdb_line.substring(30, 38));
  var y = parseFloat(pdb_line.substring(38, 46));
  var z = parseFloat(pdb_line.substring(46, 54));
  this.xyz = [x, y, z];
  this.occ = parseFloat(pdb_line.substring(54, 60));
  this.b = parseFloat(pdb_line.substring(60, 66));
  if (pdb_line.length >= 78) {
    this.element = pdb_line.substring(76, 78).trim().toUpperCase();
  }
  if (pdb_line.length >= 80) {
    this.charge = pdb_line.substring(78, 80).trim();
  }
  this.is_ligand = (NOT_LIGANDS.indexOf(this.resname) === -1);
};

Atom.prototype.resid = function () {
  return this.resseq + this.icode;
};

Atom.prototype.b_as_u = function () {
  // B = 8 * pi^2 * u^2
  return Math.sqrt(this.b / (8 * 3.14159 * 3.14159));
};

Atom.prototype.distance_sq = function (other) {
  var dx = this.xyz[0] - other.xyz[0];
  var dy = this.xyz[1] - other.xyz[1];
  var dz = this.xyz[2] - other.xyz[2];
  return dx*dx + dy*dy + dz*dz;
};

Atom.prototype.distance = function (other) {
  return Math.sqrt(this.distance_sq(other));
};

Atom.prototype.midpoint = function (other) {
  return [(this.xyz[0] + other.xyz[0]) / 2,
          (this.xyz[1] + other.xyz[1]) / 2,
          (this.xyz[2] + other.xyz[2]) / 2];
};

Atom.prototype.is_hydrogen = function () {
  return this.element === 'H' || this.element === 'D';
};

Atom.prototype.is_s_or_p = function () {
  return this.element === 'S' || this.element === 'P';
};

Atom.prototype.is_ion = function () {
  return this.element === this.resname;
};

Atom.prototype.is_water = function () {
  return this.resname === 'HOH';
};

Atom.prototype.is_same_residue = function (other, ignore_altloc) {
  return other.resseq === this.resseq && other.icode === this.icode &&
         other.chain === this.chain && other.resname === this.resname &&
         (ignore_altloc || other.altloc === this.altloc);
};

Atom.prototype.is_same_conformer = function (other) {
  return this.altloc === '' || other.altloc === '' ||
         this.altloc === other.altloc;
};

Atom.prototype.is_main_conformer = function () {
  return this.altloc === '' || this.altloc === 'A';
};

Atom.prototype.is_bonded_to = function (other) {
  /** @const */ var MAX_DIST_SQ = 1.99 * 1.99;
  /** @const */ var MAX_DIST_H_SQ = 1.3 * 1.3;
  /** @const */ var MAX_DIST_SP_SQ = 2.2 * 2.2;

  if (!this.is_same_conformer(other)) { return false; }
  if (this.element === 'H' && other.element === 'H') { return false; }
  var dxyz2 = this.distance_sq(other);
  if (dxyz2 > MAX_DIST_SP_SQ) { return false; }
  if (this.element === 'H' || other.element === 'H') {
    return dxyz2 <= MAX_DIST_H_SQ;
  }
  return dxyz2 <= MAX_DIST_SQ || this.is_s_or_p() || other.is_s_or_p();
};

Atom.prototype.long_label = function () {
  var a = this;
  return a.name + ' /' + a.resseq + ' ' + a.resname + '/' + a.chain +
         ' - occ: ' + a.occ.toFixed(2) + ' bf: ' + a.b.toFixed(2) +
         ' ele: ' + a.element + ' pos: (' + a.xyz[0].toFixed(2) + ',' +
         a.xyz[0].toFixed(2) + ',' + a.xyz[2].toFixed(2) + ')';
};


// Partition atoms into boxes for quick neighbor searching.
function Cubicles(atoms, box_length, lower_bound, upper_bound) {
  var i;
  this.boxes = [];
  this.xdim = Math.ceil((upper_bound[0] - lower_bound[0]) / box_length);
  this.ydim = Math.ceil((upper_bound[1] - lower_bound[1]) / box_length);
  this.zdim = Math.ceil((upper_bound[2] - lower_bound[2]) / box_length);
  //console.log("Cubicles: " + this.xdim + "x" + this.ydim + "x" + this.zdim);
  var nxyz = this.xdim * this.ydim * this.zdim;
  for (i = 0; i < nxyz; i++) {
    this.boxes.push([]);
  }

  this.find_box_id = function (x, y, z) {
    var xstep = Math.floor((x - lower_bound[0]) / box_length);
    var ystep = Math.floor((y - lower_bound[1]) / box_length);
    var zstep = Math.floor((z - lower_bound[2]) / box_length);
    var box_id = (zstep * this.ydim + ystep) * this.xdim + xstep;
    return (box_id >= 0 && box_id < this.boxes.length) ? box_id : null;
  };

  for (i = 0; i < atoms.length; i++) {
    var xyz = atoms[i].xyz;
    var box_id = this.find_box_id(xyz[0], xyz[1], xyz[2]);
    if (box_id === null) {
      throw Error('wrong cubicle');
    }
    this.boxes[box_id].push(i);
  }
}

Cubicles.prototype.get_nearby_atoms = function (box_id) {
  var indices = [];
  var xydim = this.xdim * this.ydim;
  var uv = Math.max(box_id % xydim, 0);
  var u = Math.max(uv % this.xdim, 0);
  var v = Math.floor(uv / this.xdim);
  var w = Math.floor(box_id / xydim);
  console.assert((w * xydim) + (v * this.xdim) + u === box_id);
  for (var iu = u-1; iu <= u+1; iu++) {
    if (iu < 0 || iu >= this.xdim) { continue; }
    for (var iv = v-1; iv <= v+1; iv++) {
      if (iv < 0 || iv >= this.ydim) { continue; }
      for (var iw = w-1; iw <= w+1; iw++) {
        if (iw < 0 || iw >= this.zdim) { continue; }
        var other_box_id = (iw * xydim) + (iv * this.xdim) + iu;
        if (other_box_id >= this.boxes.length || other_box_id < 0) {
          throw Error('Box out of bounds: ID ' + other_box_id);
        }
        var box = this.boxes[other_box_id];
        for (var i = 0; i < box.length; i++) {
          indices.push(box[i]);
        }
      }
    }
  }
  return indices;
};

Model.prototype.calculate_connectivity = function () {
  var atoms = this.atoms;
  var cubes = new Cubicles(atoms, 3.0, this.lower_bound, this.upper_bound);
  var cnt = 0;
  for (var i = 0; i < cubes.boxes.length; i++) {
    var box = cubes.boxes[i];
    if (box.length === 0) { continue; }
    var nearby_atoms = cubes.get_nearby_atoms(i);
    for (var a = 0; a < box.length; a++) {
      var atom_id = box[a];
      for (var k = 0; k < nearby_atoms.length; k++) {
        var j = nearby_atoms[k];
        if (j > atom_id && atoms[atom_id].is_bonded_to(atoms[j])) {
          atoms[atom_id].bonds.push(j);
          atoms[j].bonds.push(atom_id);
          cnt++;
        }
      }
    }
  }
  //console.log(atoms.length + ' atoms, ' + cnt + ' bonds.');
  this.cubes = cubes;
};

Model.prototype.get_nearest_atom = function (x, y, z) {
  var box_id = this.cubes.find_box_id(x, y, z);
  var indices = this.cubes.get_nearby_atoms(box_id);
  var nearest = null;
  var min_d2 = Infinity;
  for (var i = 0; i < indices.length; i++) {
    var atom = this.atoms[indices[i]];
    var dx = atom.xyz[0] - x;
    var dy = atom.xyz[1] - y;
    var dz = atom.xyz[2] - z;
    var d2 = dx*dx + dy*dy + dz*dz;
    if (d2 < min_d2) {
      nearest = atom;
      min_d2 = d2;
    }
  }
  return nearest;
};

return Model;
})();

if (typeof module !== 'undefined') { module.exports = Model; }
