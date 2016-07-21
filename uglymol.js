/* UglyMol - macromolecular viewer for crystallographers, a fork of xtal.js.
 * https://uglymol.github.io

The MIT License (MIT)

Copyright (c) 2014 Nat Echols
Copyright (c) 2016 Diamond Light Source Ltd
Copyright (c) 2016 Marcin Wojdyr

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

var UGLYMOL_VERSION = '0.1.0'
var THREE = THREE || require('three'); // eslint-disable-line

var UnitCell = (function () {
'use strict';

// eslint-disable-next-line max-params
function UnitCell(a /*:number*/, b /*:number*/, c /*:number*/,
                  alpha /*:number*/, beta /*:number*/, gamma /*:number*/) {
  if (a <= 0 || b <= 0 || c <= 0 || alpha <= 0 || beta <= 0 || gamma <= 0) {
    throw Error('Zero or negative unit cell parameter(s).');
  }
  this.parameters = [a, b, c, alpha, beta, gamma];
  var deg2rad = Math.PI / 180.0;
  var cos_alpha = Math.cos(deg2rad * alpha);
  var cos_beta = Math.cos(deg2rad * beta);
  var cos_gamma = Math.cos(deg2rad * gamma);
  var sin_alpha = Math.sin(deg2rad * alpha);
  var sin_beta = Math.sin(deg2rad * beta);
  var sin_gamma = Math.sin(deg2rad * gamma);
  if (sin_alpha === 0 || sin_beta === 0 || sin_gamma === 0) {
    throw Error('Impossible angle - N*180deg.');
  }
  var cos_alpha_star_sin_beta = (cos_beta * cos_gamma - cos_alpha) / sin_gamma;
  var cos_alpha_star = cos_alpha_star_sin_beta / sin_beta;
  var s1rca2 = Math.sqrt(1.0 - cos_alpha_star * cos_alpha_star);
  // The orthogonalization matrix we use is described in ITfC B p.262:
  // "An alternative mode of orthogonalization, used by the Protein
  // Data Bank and most programs, is to align the a1 axis of the unit
  // cell with the Cartesian X_1 axis, and to align the a*_3 axis with the
  // Cartesian X_3 axis."
  /* eslint-disable no-multi-spaces, comma-spacing */
  var orth = [a,   b * cos_gamma,  c * cos_beta,
              0.0, b * sin_gamma, -c * cos_alpha_star_sin_beta,
              0.0, 0.0          ,  c * sin_beta * s1rca2];
  // based on xtal.js which is based on cctbx.uctbx
  var frac = [
    1.0 / a,
    -cos_gamma / (sin_gamma * a),
    -(cos_gamma * cos_alpha_star_sin_beta + cos_beta * sin_gamma) /
        (sin_beta * s1rca2 * sin_gamma * a),
    0.0,
    1.0 / (sin_gamma * b),
    cos_alpha_star / (s1rca2 * sin_gamma * b),
    0.0,
    0.0,
    1.0 / (sin_beta * s1rca2 * c)
  ];

  function multiply(xyz, mat) {
    var x = xyz[0], y = xyz[1], z = xyz[2];  // eslint-disable-line
    return [mat[0] * x + mat[1] * y + mat[2] * z,
                         mat[4] * y + mat[5] * z,
                                      mat[8] * z];
  }

  this.fractionalize = function (xyz) { return multiply(xyz, frac); };
  this.orthogonalize = function (xyz) { return multiply(xyz, orth); };
}

return UnitCell;
})();



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
      if (v < lower[j]) lower[j] = v;
      if (v > upper[j]) upper[j] = v;
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
    if (!a.is_main_conformer()) continue;
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
    if (atom.altloc !== '' && atom.altloc !== 'A') continue;
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
  this.resseq = null;
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
  this.resseq = parseInt(pdb_line.substring(22, 26), 10);
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

  if (!this.is_same_conformer(other)) return false;
  if (this.element === 'H' && other.element === 'H') return false;
  var dxyz2 = this.distance_sq(other);
  if (dxyz2 > MAX_DIST_SP_SQ) return false;
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
         a.xyz[1].toFixed(2) + ',' + a.xyz[2].toFixed(2) + ')';
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
    if (iu < 0 || iu >= this.xdim) continue;
    for (var iv = v-1; iv <= v+1; iv++) {
      if (iv < 0 || iv >= this.ydim) continue;
      for (var iw = w-1; iw <= w+1; iw++) {
        if (iw < 0 || iw >= this.zdim) continue;
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
  //var cnt = 0;
  for (var i = 0; i < cubes.boxes.length; i++) {
    var box = cubes.boxes[i];
    if (box.length === 0) continue;
    var nearby_atoms = cubes.get_nearby_atoms(i);
    for (var a = 0; a < box.length; a++) {
      var atom_id = box[a];
      for (var k = 0; k < nearby_atoms.length; k++) {
        var j = nearby_atoms[k];
        if (j > atom_id && atoms[atom_id].is_bonded_to(atoms[j])) {
          atoms[atom_id].bonds.push(j);
          atoms[j].bonds.push(atom_id);
          //cnt++;
        }
      }
    }
  }
  //console.log(atoms.length + ' atoms, ' + cnt + ' bonds.');
  this.cubes = cubes;
};

Model.prototype.get_nearest_atom = function (x, y, z, atom_name) {
  var box_id = this.cubes.find_box_id(x, y, z);
  var indices = this.cubes.get_nearby_atoms(box_id);
  var nearest = null;
  var min_d2 = Infinity;
  for (var i = 0; i < indices.length; i++) {
    var atom = this.atoms[indices[i]];
    if (atom_name !== undefined && atom_name !== null &&
        atom_name !== atom.name) continue;
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



var ElMap = (function () {
'use strict';

function GridArray(dim) {
  this.dim = dim; // dimensions of the grid for the entire unit cell
  this.values = new Float32Array(dim[0] * dim[1] * dim[2]);
}

function modulo(a, b) {
  var reminder = a % b;
  return reminder >= 0 ? reminder : reminder + b;
}

GridArray.prototype.grid2index = function (i, j, k) {
  i = modulo(i, this.dim[0]);
  j = modulo(j, this.dim[1]);
  k = modulo(k, this.dim[2]);
  return this.dim[2] * (this.dim[1] * i + j) + k;
};

GridArray.prototype.grid2frac = function (i, j, k) {
  return [i / this.dim[0], j / this.dim[1], k / this.dim[2]];
};

// return grid coordinates (rounded down) for the given fractional coordinates
GridArray.prototype.frac2grid = function (xyz) {
  // at one point "| 0" here made extract_block() 40% faster on V8 3.14,
  // but I don't see any effect now
  return [Math.floor(xyz[0] * this.dim[0]) | 0,
          Math.floor(xyz[1] * this.dim[1]) | 0,
          Math.floor(xyz[2] * this.dim[2]) | 0];
};

GridArray.prototype.set_grid_value = function (i, j, k, value) {
  var idx = this.grid2index(i, j, k);
  this.values[idx] = value;
};

GridArray.prototype.get_grid_value = function (i, j, k) {
  var idx = this.grid2index(i, j, k);
  return this.values[idx];
};

function ElMap() {
  this.unit_cell = null;
  this.grid = null;
  this.mean = 0.0;
  this.rms = 1.0;
}

ElMap.prototype.calculate_stddev = function (a, offset) {
  var sum = 0;
  var sq_sum = 0;
  var alen = a.length;
  for (var i = offset; i < alen; i++) {
    sum += a[i];
    sq_sum += a[i] * a[i];
  }
  var mean = sum / (alen - offset);
  var variance = sq_sum / (alen - offset) - mean * mean;
  this.mean = mean;
  this.rms = Math.sqrt(variance);
};

ElMap.prototype.abs_level = function (sigma) {
  return sigma * this.rms + this.mean;
};

// http://www.ccp4.ac.uk/html/maplib.html#description
// eslint-disable-next-line complexity
ElMap.prototype.from_ccp4 = function (buf) {
  if (buf.byteLength < 1024) throw Error('File shorter than 1024 bytes.');
  //console.log('buf type: ' + Object.prototype.toString.call(buf));
  // for now we assume both file and host are little endian
  var iview = new Int32Array(buf, 0, 256);
  // word 53 - character string 'MAP ' to identify file type
  if (iview[52] !== 0x2050414d) throw Error('not a CCP4 map');
  // map has 3 dimensions referred to as columns (fastest changing), rows
  // and sections (c-r-s)
  var n_crs = [iview[0], iview[1], iview[2]];
  this.mode = iview[3];
  var nb;
  if (this.mode === 2) nb = 4;
  else if (this.mode === 0) nb = 1;
  else throw Error('Only Mode 2 and Mode 0 of CCP4 map is supported.');
  var start = [iview[4], iview[5], iview[6]];
  var n_grid = [iview[7], iview[8], iview[9]];
  var nsymbt = iview[23]; // size of extended header in bytes
  if (1024 + nsymbt + nb*n_crs[0]*n_crs[1]*n_crs[2] !== buf.byteLength) {
    throw Error('ccp4 file too short or too long');
  }
  var fview = new Float32Array(buf, 0, buf.byteLength >> 2);
  this.unit_cell = new UnitCell(fview[10], fview[11], fview[12],
                                fview[13], fview[14], fview[15]);
  // MAPC, MAPR, MAPS - axis corresp to cols, rows, sections (1,2,3 for X,Y,Z)
  var map_crs = [iview[16], iview[17], iview[18]];
  var ax = map_crs.indexOf(1);
  var ay = map_crs.indexOf(2);
  var az = map_crs.indexOf(3);

  this.min = fview[19];
  this.max = fview[20];
  this.sg_number = iview[22];
  this.lskflg = iview[24];
  //console.log('map mean and rms:', this.mean.toFixed(4), this.rms.toFixed(4));
  this.grid = new GridArray(n_grid);
  if (nsymbt % 4 !== 0) {
    throw Error('CCP4 map with NSYMBT not divisible by 4 is not supported.');
  }
  var data_view;
  if (this.mode === 2) data_view = fview;
  else /* this.mode === 0 */ data_view = new Int8Array(buf);
  var idx = (1024 + nsymbt) / nb | 0;

  // We assume that if DMEAN and RMS from the header are not clearly wrong
  // they are what the user wants. Because the map can cover a small part
  // of the asu and its rmsd may be different than the total rmsd.
  this.mean = fview[21];
  this.rms = fview[54];
  if (this.mean < this.min || this.mean > this.max || this.rms <= 0) {
    this.calculate_stddev(data_view, idx);
  }
  var b1 = 1;
  var b0 = 0;
  // if the file was converted by mapmode2to0 - scale the data
  if (this.mode === 0 && iview[39] === -128 && iview[40] === 127) {
    // scaling f(x)=b1*x+b0 such that f(-128)=min and f(127)=max
    b1 = (this.max - this.min) / 255.0;
    b0 = 0.5 * (this.min + this.max + b1);
  }

  var end = [start[0] + n_crs[0], start[1] + n_crs[1], start[2] + n_crs[2]];
  var it = [0, 0, 0];
  for (it[2] = start[2]; it[2] < end[2]; it[2]++) { // sections
    for (it[1] = start[1]; it[1] < end[1]; it[1]++) { // rows
      for (it[0] = start[0]; it[0] < end[0]; it[0]++) { // cols
        this.grid.set_grid_value(it[ax], it[ay], it[az],
                                 b1 * data_view[idx] + b0);
        idx++;
      }
    }
  }
  if (nsymbt > 0) {
    var u8view = new Uint8Array(buf);
    for (var i = 0; i+80 <= nsymbt; i += 80) {
      var j;
      var symop = '';
      for (j = 0; j < 80; ++j) {
        symop += String.fromCharCode(u8view[1024 + i + j]);
      }
      if (/^\s*x\s*,\s*y\s*,\s*z\s*$/i.test(symop)) continue;  // skip x,y,z
      //console.log('sym ops', symop.trim());
      var mat = parse_symop(symop);
      // Note: we apply here symops to grid points instead of coordinates.
      // In the cases we came across it is equivalent, but in general not.
      for (j = 0; j < 3; ++j) {
        mat[j][3] = Math.round(mat[j][3] * n_grid[j]) | 0;
      }
      idx = (1024 + nsymbt) / nb | 0;
      var xyz = [0, 0, 0];
      for (it[2] = start[2]; it[2] < end[2]; it[2]++) { // sections
        for (it[1] = start[1]; it[1] < end[1]; it[1]++) { // rows
          for (it[0] = start[0]; it[0] < end[0]; it[0]++) { // cols
            for (j = 0; j < 3; ++j) {
              xyz[j] = it[ax] * mat[j][0] + it[ay] * mat[j][1] +
                       it[az] * mat[j][2] + mat[j][3];
            }
            this.grid.set_grid_value(xyz[0], xyz[1], xyz[2],
                                     b1 * data_view[idx] + b0);
            idx++;
          }
        }
      }
    }
  }
};

// symop -> matrix ([x,y,z] = matrix * [x,y,z,1])
function parse_symop(symop) {
  var ops = symop.toLowerCase().replace(/\s+/g, '').split(',');
  if (ops.length !== 3) throw Error('Unexpected symop: ' + symop);
  var mat = [];
  for (var i = 0; i < 3; i++) {
    var terms = ops[i].split(/(?=[+-])/);
    var row = [0, 0, 0, 0];
    for (var j = 0; j < terms.length; j++) {
      var term = terms[j];
      var sign = (term[0] === '-' ? -1 : 1);
      var m = terms[j].match(/^[+-]?([xyz])$/);
      if (m) {
        var pos = {x: 0, y: 1, z: 2}[m[1]];
        row[pos] = sign;
      } else {
        m = terms[j].match(/^[+-]?(\d)\/(\d)$/);
        if (!m) throw Error('What is ' + terms[j] + ' in ' + symop);
        row[3] = sign * Number(m[1]) / Number(m[2]);
      }
    }
    mat.push(row);
  }
  return mat;
}

// DSN6 MAP FORMAT
// http://www.uoxray.uoregon.edu/tnt/manual/node104.html
// Density values are stored as bytes.
ElMap.prototype.from_dsn6 = function (buf) {
  //console.log('buf type: ' + Object.prototype.toString.call(buf));
  var u8data = new Uint8Array(buf);
  var iview = new Int16Array(u8data.buffer);
  if (iview[18] !== 100) {
    var len = iview.length;  // or only header, 256?
    for (var n = 0; n < len; n++) {
      // swapping bytes with Uint8Array like this:
      // var tmp=u8data[n*2]; u8data[n*2]=u8data[n*2+1]; u8data[n*2+1]=tmp;
      // was slowing down this whole function 5x times (!?) on V8.
      var val = iview[n];
      iview[n] = ((val & 0xff) << 8) | ((val >> 8) & 0xff);
    }
  }
  if (iview[18] !== 100) {
    throw Error('Endian swap failed');
  }
  var origin = [iview[0], iview[1], iview[2]];
  var n_real = [iview[3], iview[4], iview[5]];
  var n_grid = [iview[6], iview[7], iview[8]];
  var cell_mult = 1.0 / iview[17];
  this.unit_cell = new UnitCell(cell_mult * iview[9],
                                cell_mult * iview[10],
                                cell_mult * iview[11],
                                cell_mult * iview[12],
                                cell_mult * iview[13],
                                cell_mult * iview[14]);
  this.grid = new GridArray(n_grid);
  var prod = iview[15] / 100;
  var plus = iview[16];
  //var data_scale_factor = iview[15] / iview[18] + iview[16];
  // bricks have 512 (8x8x8) values
  var offset = 512;
  var n_blocks = [Math.ceil(n_real[0] / 8),
                  Math.ceil(n_real[1] / 8),
                  Math.ceil(n_real[2] / 8)];
  for (var zz = 0; zz < n_blocks[2]; zz++) {
    for (var yy = 0; yy < n_blocks[1]; yy++) {
      for (var xx = 0; xx < n_blocks[0]; xx++) { // loop over bricks
        for (var k = 0; k < 8; k++) {
          var z = 8 * zz + k;
          for (var j = 0; j < 8; j++) {
            var y = 8 * yy + j;
            for (var i = 0; i < 8; i++) { // loop inside brick
              var x = 8 * xx + i;
              if (x < n_real[0] && y < n_real[1] && z < n_real[2]) {
                var density = (u8data[offset] - plus) / prod;
                offset++;
                this.grid.set_grid_value(origin[0] + x,
                                         origin[1] + y,
                                         origin[2] + z, density);
              } else {
                offset += 8 - i;
                break;
              }
            }
          }
        }
      }
    }
  }
  this.calculate_stddev(this.grid.values, 0);
  //this.show_debug_info();
};

ElMap.prototype.show_debug_info = function () {
  console.log('unit cell: ' + this.unit_cell.parameters.join(', '));
  console.log('grid: ' + this.grid.dim);
};

// Extract a block of density for calculating an isosurface using the
// separate marching cubes implementation.
ElMap.prototype.extract_block = function (radius, center) {
  var grid = this.grid;
  var unit_cell = this.unit_cell;
  var grid_min = [0, 0, 0];
  var grid_max = grid.dim;
  if (center) {
    var xyz_min = [center[0] - radius, center[1] - radius, center[2] - radius];
    var xyz_max = [center[0] + radius, center[1] + radius, center[2] + radius];
    var frac_min = unit_cell.fractionalize(xyz_min);
    var frac_max = unit_cell.fractionalize(xyz_max);
    grid_min = grid.frac2grid(frac_min);
    grid_max = grid.frac2grid(frac_max);
  }
  var nx = grid_max[0] - grid_min[0] + 1;
  var ny = grid_max[1] - grid_min[1] + 1;
  var nz = grid_max[2] - grid_min[2] + 1;
  var points = [];
  var values = [];
  for (var i = grid_min[0]; i <= grid_max[0]; i++) {
    for (var j = grid_min[1]; j <= grid_max[1]; j++) {
      for (var k = grid_min[2]; k <= grid_max[2]; k++) {
        var frac = grid.grid2frac(i, j, k);
        var orth = unit_cell.orthogonalize(frac);
        points.push(orth);
        var map_value = grid.get_grid_value(i, j, k);
        values.push(map_value);
      }
    }
  }
  this.block = {points: points, values: values, size: [nx, ny, nz]};
};

ElMap.prototype.isomesh_in_block = function (sigma, method) {
  var abs_level = this.abs_level(sigma);
  var bl = this.block;
  return isosurface(bl.size, bl.values, bl.points, abs_level, method);
};

return ElMap;
})();


// based on xtal.js/marchingcubes.js which is
// based on http://stemkoski.github.io/Three.js/Marching-Cubes.html

var isosurface = (function () {
'use strict';
/* eslint comma-spacing: 0, no-multi-spaces: 0 */

// Marching cubes lookup tables.
// These tables are straight from Paul Bourke's page:
// http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
// who in turn got them from Cory Gene Bloyd.
var edgeTable = new Int32Array([
  0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
  0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
  0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
  0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
  0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
  0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
  0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
  0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
  0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
  0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
  0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
  0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
  0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
  0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
  0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
  0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
  0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
  0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
  0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
  0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
  0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
  0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
  0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
  0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
  0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
  0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
  0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
  0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
  0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
  0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
  0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
  0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0]);

var triTable = [
  [],
  [0, 8, 3],
  [0, 1, 9],
  [1, 8, 3, 9, 8, 1],
  [1, 2, 10],
  [0, 8, 3, 1, 2, 10],
  [9, 2, 10, 0, 2, 9],
  [2, 8, 3, 2, 10, 8, 10, 9, 8],
  [3, 11, 2],
  [0, 11, 2, 8, 11, 0],
  [1, 9, 0, 2, 3, 11],
  [1, 11, 2, 1, 9, 11, 9, 8, 11],
  [3, 10, 1, 11, 10, 3],
  [0, 10, 1, 0, 8, 10, 8, 11, 10],
  [3, 9, 0, 3, 11, 9, 11, 10, 9],
  [9, 8, 10, 10, 8, 11],
  [4, 7, 8],
  [4, 3, 0, 7, 3, 4],
  [0, 1, 9, 8, 4, 7],
  [4, 1, 9, 4, 7, 1, 7, 3, 1],
  [1, 2, 10, 8, 4, 7],
  [3, 4, 7, 3, 0, 4, 1, 2, 10],
  [9, 2, 10, 9, 0, 2, 8, 4, 7],
  [2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4],
  [8, 4, 7, 3, 11, 2],
  [11, 4, 7, 11, 2, 4, 2, 0, 4],
  [9, 0, 1, 8, 4, 7, 2, 3, 11],
  [4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1],
  [3, 10, 1, 3, 11, 10, 7, 8, 4],
  [1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4],
  [4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3],
  [4, 7, 11, 4, 11, 9, 9, 11, 10],
  [9, 5, 4],
  [9, 5, 4, 0, 8, 3],
  [0, 5, 4, 1, 5, 0],
  [8, 5, 4, 8, 3, 5, 3, 1, 5],
  [1, 2, 10, 9, 5, 4],
  [3, 0, 8, 1, 2, 10, 4, 9, 5],
  [5, 2, 10, 5, 4, 2, 4, 0, 2],
  [2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8],
  [9, 5, 4, 2, 3, 11],
  [0, 11, 2, 0, 8, 11, 4, 9, 5],
  [0, 5, 4, 0, 1, 5, 2, 3, 11],
  [2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5],
  [10, 3, 11, 10, 1, 3, 9, 5, 4],
  [4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10],
  [5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3],
  [5, 4, 8, 5, 8, 10, 10, 8, 11],
  [9, 7, 8, 5, 7, 9],
  [9, 3, 0, 9, 5, 3, 5, 7, 3],
  [0, 7, 8, 0, 1, 7, 1, 5, 7],
  [1, 5, 3, 3, 5, 7],
  [9, 7, 8, 9, 5, 7, 10, 1, 2],
  [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3],
  [8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2],
  [2, 10, 5, 2, 5, 3, 3, 5, 7],
  [7, 9, 5, 7, 8, 9, 3, 11, 2],
  [9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11],
  [2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7],
  [11, 2, 1, 11, 1, 7, 7, 1, 5],
  [9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11],
  [5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0],
  [11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0],
  [11, 10, 5, 7, 11, 5],
  [10, 6, 5],
  [0, 8, 3, 5, 10, 6],
  [9, 0, 1, 5, 10, 6],
  [1, 8, 3, 1, 9, 8, 5, 10, 6],
  [1, 6, 5, 2, 6, 1],
  [1, 6, 5, 1, 2, 6, 3, 0, 8],
  [9, 6, 5, 9, 0, 6, 0, 2, 6],
  [5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8],
  [2, 3, 11, 10, 6, 5],
  [11, 0, 8, 11, 2, 0, 10, 6, 5],
  [0, 1, 9, 2, 3, 11, 5, 10, 6],
  [5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11],
  [6, 3, 11, 6, 5, 3, 5, 1, 3],
  [0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6],
  [3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9],
  [6, 5, 9, 6, 9, 11, 11, 9, 8],
  [5, 10, 6, 4, 7, 8],
  [4, 3, 0, 4, 7, 3, 6, 5, 10],
  [1, 9, 0, 5, 10, 6, 8, 4, 7],
  [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4],
  [6, 1, 2, 6, 5, 1, 4, 7, 8],
  [1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7],
  [8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6],
  [7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9],
  [3, 11, 2, 7, 8, 4, 10, 6, 5],
  [5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11],
  [0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6],
  [9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6],
  [8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6],
  [5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11],
  [0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7],
  [6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9],
  [10, 4, 9, 6, 4, 10],
  [4, 10, 6, 4, 9, 10, 0, 8, 3],
  [10, 0, 1, 10, 6, 0, 6, 4, 0],
  [8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10],
  [1, 4, 9, 1, 2, 4, 2, 6, 4],
  [3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4],
  [0, 2, 4, 4, 2, 6],
  [8, 3, 2, 8, 2, 4, 4, 2, 6],
  [10, 4, 9, 10, 6, 4, 11, 2, 3],
  [0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6],
  [3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10],
  [6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1],
  [9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3],
  [8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1],
  [3, 11, 6, 3, 6, 0, 0, 6, 4],
  [6, 4, 8, 11, 6, 8],
  [7, 10, 6, 7, 8, 10, 8, 9, 10],
  [0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10],
  [10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0],
  [10, 6, 7, 10, 7, 1, 1, 7, 3],
  [1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7],
  [2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9],
  [7, 8, 0, 7, 0, 6, 6, 0, 2],
  [7, 3, 2, 6, 7, 2],
  [2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7],
  [2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7],
  [1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11],
  [11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1],
  [8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6],
  [0, 9, 1, 11, 6, 7],
  [7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0],
  [7, 11, 6],
  [7, 6, 11],
  [3, 0, 8, 11, 7, 6],
  [0, 1, 9, 11, 7, 6],
  [8, 1, 9, 8, 3, 1, 11, 7, 6],
  [10, 1, 2, 6, 11, 7],
  [1, 2, 10, 3, 0, 8, 6, 11, 7],
  [2, 9, 0, 2, 10, 9, 6, 11, 7],
  [6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8],
  [7, 2, 3, 6, 2, 7],
  [7, 0, 8, 7, 6, 0, 6, 2, 0],
  [2, 7, 6, 2, 3, 7, 0, 1, 9],
  [1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6],
  [10, 7, 6, 10, 1, 7, 1, 3, 7],
  [10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8],
  [0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7],
  [7, 6, 10, 7, 10, 8, 8, 10, 9],
  [6, 8, 4, 11, 8, 6],
  [3, 6, 11, 3, 0, 6, 0, 4, 6],
  [8, 6, 11, 8, 4, 6, 9, 0, 1],
  [9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6],
  [6, 8, 4, 6, 11, 8, 2, 10, 1],
  [1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6],
  [4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9],
  [10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3],
  [8, 2, 3, 8, 4, 2, 4, 6, 2],
  [0, 4, 2, 4, 6, 2],
  [1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8],
  [1, 9, 4, 1, 4, 2, 2, 4, 6],
  [8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1],
  [10, 1, 0, 10, 0, 6, 6, 0, 4],
  [4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3],
  [10, 9, 4, 6, 10, 4],
  [4, 9, 5, 7, 6, 11],
  [0, 8, 3, 4, 9, 5, 11, 7, 6],
  [5, 0, 1, 5, 4, 0, 7, 6, 11],
  [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5],
  [9, 5, 4, 10, 1, 2, 7, 6, 11],
  [6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5],
  [7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2],
  [3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6],
  [7, 2, 3, 7, 6, 2, 5, 4, 9],
  [9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7],
  [3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0],
  [6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8],
  [9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7],
  [1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4],
  [4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10],
  [7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10],
  [6, 9, 5, 6, 11, 9, 11, 8, 9],
  [3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5],
  [0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11],
  [6, 11, 3, 6, 3, 5, 5, 3, 1],
  [1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6],
  [0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10],
  [11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5],
  [6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3],
  [5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2],
  [9, 5, 6, 9, 6, 0, 0, 6, 2],
  [1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8],
  [1, 5, 6, 2, 1, 6],
  [1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6],
  [10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0],
  [0, 3, 8, 5, 6, 10],
  [10, 5, 6],
  [11, 5, 10, 7, 5, 11],
  [11, 5, 10, 11, 7, 5, 8, 3, 0],
  [5, 11, 7, 5, 10, 11, 1, 9, 0],
  [10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1],
  [11, 1, 2, 11, 7, 1, 7, 5, 1],
  [0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11],
  [9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7],
  [7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2],
  [2, 5, 10, 2, 3, 5, 3, 7, 5],
  [8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5],
  [9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2],
  [9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2],
  [1, 3, 5, 3, 7, 5],
  [0, 8, 7, 0, 7, 1, 1, 7, 5],
  [9, 0, 3, 9, 3, 5, 5, 3, 7],
  [9, 8, 7, 5, 9, 7],
  [5, 8, 4, 5, 10, 8, 10, 11, 8],
  [5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0],
  [0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5],
  [10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4],
  [2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8],
  [0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11],
  [0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5],
  [9, 4, 5, 2, 11, 3],
  [2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4],
  [5, 10, 2, 5, 2, 4, 4, 2, 0],
  [3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9],
  [5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2],
  [8, 4, 5, 8, 5, 3, 3, 5, 1],
  [0, 4, 5, 1, 0, 5],
  [8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5],
  [9, 4, 5],
  [4, 11, 7, 4, 9, 11, 9, 10, 11],
  [0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11],
  [1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11],
  [3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4],
  [4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2],
  [9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3],
  [11, 7, 4, 11, 4, 2, 2, 4, 0],
  [11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4],
  [2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9],
  [9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7],
  [3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10],
  [1, 10, 2, 8, 7, 4],
  [4, 9, 1, 4, 1, 7, 7, 1, 3],
  [4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1],
  [4, 0, 3, 7, 4, 3],
  [4, 8, 7],
  [9, 10, 8, 10, 11, 8],
  [3, 0, 9, 3, 9, 11, 11, 9, 10],
  [0, 1, 10, 0, 10, 8, 8, 10, 11],
  [3, 1, 10, 11, 3, 10],
  [1, 2, 11, 1, 11, 9, 9, 11, 8],
  [3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9],
  [0, 2, 11, 8, 0, 11],
  [3, 2, 11],
  [2, 3, 8, 2, 8, 10, 10, 8, 9],
  [9, 10, 2, 0, 9, 2],
  [2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8],
  [1, 10, 2],
  [1, 3, 8, 9, 1, 8],
  [0, 9, 1],
  [0, 3, 8],
  []];

var cubeVerts = [[0,0,0], [1,0,0], [1,1,0], [0,1,0],
                 [0,0,1], [1,0,1], [1,1,1], [0,1,1]];
var edgeIndex = [[0,1], [1,2], [2,3], [3,0], [4,5], [5,6],
                 [6,7], [7,4], [0,4], [1,5], [2,6], [3,7]];


function check_input(dims, values, points) {
  if (dims[0] <= 0 || dims[1] <= 0 || dims[2] <= 0) {
    throw Error('Grid dimensions are zero along at least one edge');
  }
  var size_xyz = dims[0] * dims[1] * dims[2];
  if (values.length !== size_xyz || points.length !== size_xyz) {
    throw Error('isosurface: array size mismatch');
  }
}

// return offsets relative to vertex [0,0,0]
function calculate_vert_offsets(dims) {
  var vert_offsets = [];
  for (var i = 0; i < 8; ++i) {
    var v = cubeVerts[i];
    vert_offsets.push(v[0] + dims[2] * (v[1] + dims[1] * v[2]));
  }
  return vert_offsets;
}


function marching_cubes(dims, values, points, isolevel, snap) {
  var vlist = new Array(12);
  var vert_offsets = calculate_vert_offsets(dims);
  var vertex_values = new Float32Array(8);
  var vertex_points = [null, null, null, null, null, null, null, null];
  var size_x = dims[0];
  var size_y = dims[1];
  var size_z = dims[2];
  var vertices = [];
  var faces = [];
  var vertex_count = 0;
  for (var x = 0; x < size_x - 1; x++) {
    for (var y = 0; y < size_y - 1; y++) {
      for (var z = 0; z < size_z - 1; z++) {
        var offset0 = z + size_z * (y + size_y * x);
        var cubeindex = 0;
        var i, j;
        for (i = 0; i < 8; ++i) {
          j = offset0 + vert_offsets[i];
          cubeindex |= (values[j] < isolevel) ? 1 << i : 0;
        }
        if (cubeindex === 0 || cubeindex === 255) continue;
        for (i = 0; i < 8; ++i) {
          j = offset0 + vert_offsets[i];
          vertex_values[i] = values[j];
          vertex_points[i] = points[j];
        }

        // 12 bit number, indicates which edges are crossed by the isosurface
        var edge_mask = edgeTable[cubeindex];

        // check which edges are crossed, and estimate the point location
        // using a weighted average of scalar values at edge endpoints.
        for (i = 0; i < 12; ++i) {
          if ((edge_mask & (1 << i)) !== 0) {
            var e = edgeIndex[i];
            var mu = (isolevel - vertex_values[e[0]]) /
                     (vertex_values[e[1]] - vertex_values[e[0]]);
            if (snap === true) {
              if (mu > 0.85) mu = 1;
              else if (mu < 0.15) mu = 0;
            }
            var p1 = vertex_points[e[0]];
            var p2 = vertex_points[e[1]];
            // TODO: avoid duplicated vertices among neighbouring cells
            vertices.push(p1[0] + (p2[0] - p1[0]) * mu,
                          p1[1] + (p2[1] - p1[1]) * mu,
                          p1[2] + (p2[2] - p1[2]) * mu);
            vlist[i] = vertex_count++;
          }
        }
        // construct triangles -- get correct vertices from triTable.
        var t = triTable[cubeindex];
        for (i = 0; i < t.length; i++) {
          faces.push(vlist[t[i]]);
        }
      }
    }
  }
  return { vertices: vertices, faces: faces };
}

function isosurface(dims, values, points, isolevel, method) {
  check_input(dims, values, points);
  //if (method === 'marching tetrahedra') {
  //  return marching_tetrahedra(dims, values, points, isolevel);
  //}
  return marching_cubes(dims, values, points, isolevel,
                        method === 'snapped MC');
}

return isosurface;
})();



var Viewer = (function () {
'use strict';

var ColorSchemes = { // accessible as Viewer.ColorSchemes
  dark: {
    bg: 0x000000,
    map_den: 0x3362B2,
    map_pos: 0x298029,
    map_neg: 0x8B2E2E,
    center: 0xC997B0,
    cell_box: 0xFFFFFF,
    // atoms
    H: 0x858585, // H is normally invisible
    // C, N and O are taken approximately (by color-picker) from coot
    C: 0xb3b300,
    N: 0x7EAAFB,
    O: 0xF24984,
    S: 0x40ff40, // S in coot is too similar to C, here it is greener
    // Coot doesn't define other colors (?)
    MG: 0xc0c0c0,
    P: 0xffc040,
    CL: 0xa0ff60,
    CA: 0xffffff,
    MN: 0xff90c0,
    FE: 0xa03000,
    NI: 0x00ff80,
    def: 0xa0a0a0 // default atom color
  },
  light: { // like in Coot after Edit > Background Color > White
    bg: 0xFFFFFF,
    map_den: 0x3362B2,
    map_pos: 0x298029,
    map_neg: 0x8B2E2E,
    center: 0xC7C769,
    cell_box: 0x000000,
    H: 0x999999,
    C: 0xA96464,
    N: 0x1C51B3,
    O: 0xC33869,
    S: 0x9E7B3D,
    def: 0x808080
  }
};

var auto_speed = 1.0;  // accessible as Viewer.auto_speed

// relative position on canvas in normalized device coordinates [-1, +1]
function relX(evt) { return 2 * evt.pageX / window.innerWidth - 1; }
function relY(evt) { return 1 - 2 * evt.pageY / window.innerHeight; }

// map 2d position to sphere with radius 1.
function project_on_ball(x, y) {
  var z = 0;
  var length_sq = x * x + y * y;
  if (length_sq < 1) {  // in ellipse
    z = Math.sqrt(1.0 - length_sq);
  } else {  // in a corner
    var length = Math.sqrt(length_sq);
    x /= length;
    y /= length;
  }
  return [x, y, z];  // guaranteed to be normalized
}

var _raycaster;
function get_raycaster(coords, camera) {
  if (_raycaster === undefined) _raycaster = new THREE.Raycaster();
  _raycaster.setFromCamera(coords, camera);
  _raycaster.near = camera.near;
  _raycaster.far = camera.far - 0.2 * (camera.far - camera.near); // 20% in fog
  _raycaster.linePrecision = 0.2;
  return _raycaster;
}

var STATE = {NONE: -1, ROTATE: 0, PAN: 1, ZOOM: 2, PAN_ZOOM: 3, SLAB: 4,
             ROLL: 5, AUTO_ROTATE: 6, GO: 7};


// based on three.js/examples/js/controls/OrthographicTrackballControls.js
var Controls = function (camera, target) {
  var _state = STATE.NONE;
  var _rotate_start = new THREE.Vector3();
  var _rotate_end = new THREE.Vector3();
  var _zoom_start = new THREE.Vector2();
  var _zoom_end = new THREE.Vector2();
  var _pinch_start = 0;
  var _pinch_end = 0;
  var _pan_start = new THREE.Vector2();
  var _pan_end = new THREE.Vector2();
  var _panned = true;
  var _slab_width = 10.0;
  var _rock_state = 0.0;
  var _go_func = null;

  function change_slab_width(delta) {
    _slab_width = Math.max(_slab_width + delta, 0.01);
  }

  function rotate_camera(eye) {
    var quat = new THREE.Quaternion();
    quat.setFromUnitVectors(_rotate_end, _rotate_start);
    eye.applyQuaternion(quat);
    camera.up.applyQuaternion(quat);
    _rotate_end.applyQuaternion(quat);
    _rotate_start.copy(_rotate_end);
  }

  function zoom_camera(eye) {
    var dx = _zoom_end.x - _zoom_start.x;
    var dy = _zoom_end.y - _zoom_start.y;
    if (_state === STATE.ZOOM) {
      camera.zoom /= (1 - dx + dy);
    } else if (_state === STATE.SLAB) {
      change_slab_width(10.0 * dx);
      target.addScaledVector(eye, -5.0 / eye.length() * dy);
    } else if (_state === STATE.ROLL) {
      camera.up.applyAxisAngle(eye, 0.05 * (dx - dy));
    }
    _zoom_start.copy(_zoom_end);
  }

  function pan_camera(eye) {
    var dx = _pan_end.x - _pan_start.x;
    var dy = _pan_end.y - _pan_start.y;
    dx *= 0.5 * (camera.right - camera.left) / camera.zoom;
    dy *= 0.5 * (camera.bottom - camera.top) / camera.zoom;
    var pan = eye.clone().cross(camera.up).setLength(dx);
    pan.addScaledVector(camera.up, dy / camera.up.length());
    camera.position.add(pan);
    target.add(pan);
    _pan_start.copy(_pan_end);
  }

  this.toggle_state = function (toggled, params) {
    _state = (_state === toggled ? STATE.NONE : toggled);
    _rock_state = params.rock ? 0.0 : null;
  };

  this.is_going = function () { return _state === STATE.GO; };

  function auto_rotate(eye) {
    _rotate_start.copy(eye).normalize();
    var speed = 3e-4 * auto_speed;
    if (_rock_state !== null) {
      _rock_state += 0.02;
      speed = 4e-5 * auto_speed * Math.cos(_rock_state);
    }
    _rotate_end.crossVectors(camera.up, eye).multiplyScalar(speed)
      .add(_rotate_start);
  }

  this.update = function () {
    var changed = false;
    var eye = camera.position.clone().sub(target);
    if (_state === STATE.AUTO_ROTATE) {
      auto_rotate(eye);
    }
    if (!_rotate_start.equals(_rotate_end)) {
      rotate_camera(eye);
      changed = true;
    }
    if (_pinch_end !== _pinch_start) {
      camera.zoom *= _pinch_end / _pinch_start;
      _pinch_start = _pinch_end;
      changed = true;
    }
    if (!_zoom_end.equals(_zoom_start)) {
      zoom_camera(eye);
      changed = true;
    }
    if (!_pan_end.equals(_pan_start)) {
      pan_camera(eye);
      _panned = true;
      changed = true;
    }
    camera.position.addVectors(target, eye);
    if (_state === STATE.GO) {
      _go_func();
      changed = true;
    }
    camera.lookAt(target);
    return changed;
  };

  this.start = function (new_state, x, y, dist) {
    if (_state === STATE.NONE || _state === STATE.AUTO_ROTATE) {
      _state = new_state;
    }
    this.move(x, y, dist);
    switch (_state) {
      case STATE.ROTATE:
        _rotate_start.copy(_rotate_end);
        break;
      case STATE.ZOOM:
      case STATE.SLAB:
      case STATE.ROLL:
        _zoom_start.copy(_zoom_end);
        break;
      case STATE.PAN:
        _pan_start.copy(_pan_end);
        _panned = false;
        break;
      case STATE.PAN_ZOOM:
        _pinch_start = _pinch_end;
        _pan_start.copy(_pan_end);
        break;
    }
  };

  this.move = function (x, y, dist) {
    switch (_state) {
      case STATE.ROTATE:
        var xyz = project_on_ball(x, y);
        //console.log(this.camera.projectionMatrix);
        //console.log(this.camera.matrixWorld);
        // TODO maybe use project()/unproject()/applyProjection()
        var eye = camera.position.clone().sub(target);
        _rotate_end.crossVectors(camera.up, eye).setLength(xyz[0]);
        _rotate_end.addScaledVector(camera.up, xyz[1] / camera.up.length());
        _rotate_end.addScaledVector(eye, xyz[2] / eye.length());
        break;
      case STATE.ZOOM:
      case STATE.SLAB:
      case STATE.ROLL:
        _zoom_end.set(x, y);
        break;
      case STATE.PAN:
        _pan_end.set(x, y);
        break;
      case STATE.PAN_ZOOM:
        _pan_end.set(x, y);
        _pinch_end = dist;
        break;
    }
  };

  this.stop = function (model_bag) {
    var atom = null;
    if (_state === STATE.PAN && !_panned && model_bag) {
      atom = model_bag.pick_atom(get_raycaster(_pan_start, camera));
    }
    _state = STATE.NONE;
    _rotate_start.copy(_rotate_end);
    _pinch_start = _pinch_end;
    _pan_start.copy(_pan_end);
    if (atom !== null) { // center on atom
      this.go_to(new THREE.Vector3(atom.xyz[0], atom.xyz[1], atom.xyz[2]));
    }
  };

  this.slab_width = function () { return _slab_width; };
  this.change_slab_width = change_slab_width;

  this.go_to = function (targ, cam_pos, cam_up, steps) {
    if ((!targ || targ.distanceToSquared(target) < 0.1) &&
        (!cam_pos || cam_pos.distanceToSquared(camera.position) < 0.1) &&
        (!cam_up || cam_up.distanceToSquared(camera.up) < 0.1)) {
      return;
    }
    _state = STATE.GO;
    steps = steps || (60 / auto_speed);
    var alphas = [];
    var prev_pos = 0;
    for (var i = 1; i <= steps; ++i) {
      var pos = i / steps;
      // quadratic easing
      pos = pos < 0.5 ? 2 * pos * pos : -2 * pos * (pos-2) - 1;
      alphas.push((pos - prev_pos) / (1 - prev_pos));
      prev_pos = pos;
    }
    _go_func = function () {
      var a = alphas.shift();
      if (targ) {
        // unspecified cam_pos - camera stays in the same distance to target
        if (!cam_pos) camera.position.sub(target);
        target.lerp(targ, a);
        if (!cam_pos) camera.position.add(target);
      }
      if (cam_pos) camera.position.lerp(cam_pos, a);
      if (cam_up) camera.up.lerp(cam_up, a);
      if (alphas.length === 0) {
        _state = STATE.NONE;
        _go_func = null;
      }
    };
  };
};


// constants

var CUBE_EDGES = [[0, 0, 0], [1, 0, 0],
                  [0, 0, 0], [0, 1, 0],
                  [0, 0, 0], [0, 0, 1],
                  [1, 0, 0], [1, 1, 0],
                  [1, 0, 0], [1, 0, 1],
                  [0, 1, 0], [1, 1, 0],
                  [0, 1, 0], [0, 1, 1],
                  [0, 0, 1], [1, 0, 1],
                  [0, 0, 1], [0, 1, 1],
                  [1, 0, 1], [1, 1, 1],
                  [1, 1, 0], [1, 1, 1],
                  [0, 1, 1], [1, 1, 1]];

var COLOR_AIMS = ['element', 'B-factor', 'occupancy', 'index', 'chain'];
var RENDER_STYLES = ['lines', 'trace', 'ribbon', 'lines+balls'];
var MAP_STYLES = ['marching cubes', 'snapped MC'];

function make_center_cube(size, ctr, color) {
  var geometry = new THREE.Geometry();
  for (var i = 0; i < CUBE_EDGES.length; i++) {
    var a = CUBE_EDGES[i];
    var x = ctr.x + size * (a[0] - 0.5);
    var y = ctr.y + size * (a[1] - 0.5);
    var z = ctr.z + size * (a[2] - 0.5);
    geometry.vertices.push(new THREE.Vector3(x, y, z));
  }
  var material = new THREE.LineBasicMaterial({color: color, linewidth: 2});
  return new THREE.LineSegments(geometry, material);
}

function make_unitcell_box(uc, color) {
  if (!uc) {
    throw Error('Unit cell not defined!');
  }
  var geometry = new THREE.Geometry();
  for (var i = 0; i < CUBE_EDGES.length; i++) {
    var xyz = uc.orthogonalize(CUBE_EDGES[i]);
    geometry.vertices.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]));
  }
  geometry.colors.push(
    new THREE.Color(0xff0000), new THREE.Color(0xffaa00),
    new THREE.Color(0x00ff00), new THREE.Color(0xaaff00),
    new THREE.Color(0x0000ff), new THREE.Color(0x00aaff)
  );
  for (var j = 6; j < CUBE_EDGES.length; j++) {
    geometry.colors.push(color);
  }
  var material = new THREE.LineBasicMaterial({vertexColors:
                                                THREE.VertexColors});
  return new THREE.LineSegments(geometry, material);
}

function rainbow_value(v, vmin, vmax) {
  var c = new THREE.Color(0xe0e0e0);
  if (vmin < vmax) {
    var ratio = (v - vmin) / (vmax - vmin);
    var hue = (240 - (240 * ratio)) / 360;
    c.setHSL(hue, 1.0, 0.5);
  }
  return c;
}

function color_by(style, atoms, elem_colors) {
  var color_func;
  var i;
  var last_atom = atoms[atoms.length-1];
  if (style === 'index') {
    color_func = function (atom) {
      return rainbow_value(atom.i_seq, 0, last_atom.i_seq);
    };
  } else if (style === 'B-factor') {
    var vmin = Infinity;
    var vmax = -Infinity;
    for (i = 0; i < atoms.length; i++) {
      var v = atoms[i].b;
      if (v > vmax) vmax = v;
      if (v < vmin) vmin = v;
    }
    console.log('B-factors in [' + vmin + ', ' + vmax + ']');
    color_func = function (atom) {
      return rainbow_value(atom.b, vmin, vmax);
    };
  } else if (style === 'occupancy') {
    color_func = function (atom) {
      return rainbow_value(atom.occ, 0, 1);
    };
  } else if (style === 'chain') {
    color_func = function (atom) {
      return rainbow_value(atom.chain_index, 0, last_atom.chain_index);
    };
  } else { // element
    color_func = function (atom) {
      return elem_colors[atom.element] || elem_colors.def;
    };
  }
  var colors = [];
  for (i = 0; i < atoms.length; i++) {
    colors.push(color_func(atoms[i]));
  }
  return colors;
}

function make_balls(visible_atoms, colors, ball_size) {
  var ball_texture = new THREE.TextureLoader().load('src/ball.png');
  var pt_geometry = new THREE.Geometry();
  for (var i = 0; i < visible_atoms.length; i++) {
    var xyz = visible_atoms[i].xyz;
    pt_geometry.vertices.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]));
    pt_geometry.colors.push(colors[i]);
  }
  var pt_material = new THREE.PointsMaterial({
    vertexColors: THREE.VertexColors,
    map: ball_texture,
    size: ball_size,
    alphaTest: 0.5
  });
  return new THREE.Points(pt_geometry, pt_material);
}

function make_segment_geometry(segment, colors, c_offset, smooth) {
  var geometry = new THREE.Geometry();
  var i, xyz;
  if (!smooth || smooth < 2) {
    for (i = 0; i < segment.length; i++) {
      xyz = segment[i].xyz;
      geometry.vertices.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]));
      geometry.colors.push(colors[c_offset+i]);
    }
  } else {
    var points = [];
    for (i = 0; i < segment.length; i++) {
      xyz = segment[i].xyz;
      points.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]));
    }
    for (i = 0; i < segment.length - 1; i++) {
      for (var j = 0; j < smooth; ++j) {
        geometry.colors.push(colors[c_offset+i]);
      }
    }
    geometry.colors.push(colors[c_offset+i]);
    var curve = new THREE.CatmullRomCurve3(points);
    geometry.vertices = curve.getPoints(geometry.colors.length - 1);
  }
  return geometry;
}

// Add a representation of an unbonded atom as a cross to geometry
function add_isolated_atom(geometry, atom, color) {
  var c = atom.xyz;
  var R = 0.7;
  geometry.vertices.push(new THREE.Vector3(c[0]-R, c[1], c[2]));
  geometry.vertices.push(new THREE.Vector3(c[0]+R, c[1], c[2]));
  geometry.vertices.push(new THREE.Vector3(c[0], c[1]-R, c[2]));
  geometry.vertices.push(new THREE.Vector3(c[0], c[1]+R, c[2]));
  geometry.vertices.push(new THREE.Vector3(c[0], c[1], c[2]-R));
  geometry.vertices.push(new THREE.Vector3(c[0], c[1], c[2]+R));
  for (var i = 0; i < 6; i++) {
    geometry.colors.push(color);
  }
}

function set_colors(palette, o) {
  var scheme = ColorSchemes[palette];
  for (var key in scheme) {
    if (o[key]) {
      o[key].set(scheme[key]);
    } else {
      o[key] = new THREE.Color(scheme[key]);
    }
  }
  o.name = palette;
  return o;
}


function MapBag(map, is_diff_map) {
  this.map = map;
  this.name = '';
  this.isolevel = is_diff_map ? 3.0 : 1.5;
  this.visible = true;
  this.types = is_diff_map ? ['map_pos', 'map_neg'] : ['map_den'];
  this.block_ctr = new THREE.Vector3(Infinity, 0, 0);
  this.el_objects = []; // three.js objects
}


function ModelBag(model, config) {
  this.model = model;
  this.name = '';
  this.visible = true;
  this.conf = config;
  this.atomic_objects = null; // list of three.js objects
}

ModelBag.prototype.pick_atom = function (raycaster) {
  var intersects = raycaster.intersectObjects(this.atomic_objects);
  if (intersects.length < 1) return null;
  var p = intersects[0].point;
  return this.model.get_nearest_atom(p.x, p.y, p.z);
};

ModelBag.prototype.get_visible_atoms = function () {
  var atoms = this.model.atoms;
  if (this.conf.hydrogens || !this.model.has_hydrogens) {
    return atoms;
  }
  var non_h = [];
  for (var i = 0; i < atoms.length; i++) {
    if (atoms[i].element !== 'H') {
      non_h.push(atoms[i]);
    }
  }
  return non_h;
};

ModelBag.prototype.add_bonds = function (ligands_only, ball_size) {
  var visible_atoms = this.get_visible_atoms();
  var color_style = ligands_only ? 'element' : this.conf.color_aim;
  var colors = color_by(color_style, visible_atoms, this.conf.colors);
  var geometry = new THREE.Geometry();
  var opt = { hydrogens: this.conf.hydrogens,
              ligands_only: ligands_only,
              balls: this.conf.render_style === 'lines+balls' };
  for (var i = 0; i < visible_atoms.length; i++) {
    var atom = visible_atoms[i];
    var color = colors[i];
    if (ligands_only && !atom.is_ligand) continue;
    if (atom.bonds.length === 0 && !opt.balls) { // nonbonded, draw star
      add_isolated_atom(geometry, atom, color);
    } else { // bonded, draw lines
      for (var j = 0; j < atom.bonds.length; j++) {
        // TODO: one line per bond (not trivial, because coloring)
        var other = this.model.atoms[atom.bonds[j]];
        if (!opt.hydrogens && other.element === 'H') continue;
        // Coot show X-H bonds as thinner lines in a single color.
        // Here we keep it simple and render such bonds like all others.
        if (opt.ligands_only && !other.is_ligand) continue;
        var mid = atom.midpoint(other);
        var vmid = new THREE.Vector3(mid[0], mid[1], mid[2]);
        var vatom = new THREE.Vector3(atom.xyz[0], atom.xyz[1], atom.xyz[2]);
        if (opt.balls) {
          vatom.lerp(vmid, 0.3); // TODO: use ball_size
        }
        geometry.vertices.push(vatom, vmid);
        geometry.colors.push(color, color);
      }
    }
  }
  var material = new THREE.LineBasicMaterial({
    vertexColors: THREE.VertexColors,
    linewidth: this.conf.line_width
  });
  //console.log('make_bonds() vertex count: ' + geometry.vertices.length);
  this.atomic_objects.push(new THREE.LineSegments(geometry, material));
  if (opt.balls) {
    this.atomic_objects.push(make_balls(visible_atoms, colors, ball_size));
  }
};

ModelBag.prototype.add_trace = function (smoothness) {
  var segments = this.model.extract_trace();
  var visible_atoms = [].concat.apply([], segments);
  var colors = color_by(this.conf.color_aim, visible_atoms, this.conf.colors);
  var k = 0;
  var material = new THREE.LineBasicMaterial({
    vertexColors: THREE.VertexColors,
    linewidth: this.conf.line_width
  });
  for (var i = 0; i < segments.length; i++) {
    var geom = make_segment_geometry(segments[i], colors, k, smoothness);
    k += segments[i].length;
    var line = new THREE.Line(geom, material);
    this.atomic_objects.push(line);
  }
};

ModelBag.prototype.add_ribbon = function () {
  // for now it's not really a ribbon
  this.add_trace(8);
};


function Viewer(element_id) {
  this.config = {
    bond_line: 4.0, // for 700px height (in Coot it also depends on height)
    map_line: 1.25,  // for any height
    map_radius: 10.0,
    map_style: MAP_STYLES[0],
    render_style: RENDER_STYLES[0],
    color_aim: COLOR_AIMS[0],
    colors: set_colors('dark', {}),
    hydrogens: false,
    line_width: 0 // it will be set in resize()
  };

  // rendered objects
  this.model_bags = [];
  this.map_bags = [];
  this.decor = {cell_box: null, selection: null};
  this.nav = null;

  this.last_ctr = new THREE.Vector3(Infinity, 0, 0);
  this.initial_hud_text = null;
  this.initial_hud_bg = '';
  this.selected_atom = null;
  this.active_model_bag = null;
  this.scene = new THREE.Scene();
  this.target = new THREE.Vector3();
  this.camera = new THREE.OrthographicCamera();
  this.scene.add(this.camera);
  this.scene.fog = new THREE.Fog(this.config.colors.bg, 0, 1);
  this.light = new THREE.AmbientLight(0xffffff);
  this.scene.add(this.light);
  this.controls = new Controls(this.camera, this.target);

  if (typeof document === 'undefined') return;  // for testing on node

  this.renderer = new THREE.WebGLRenderer({antialias: true});
  this.renderer.setClearColor(this.config.colors.bg, 1);
  this.renderer.setPixelRatio(window.devicePixelRatio);
  this.resize();
  this.camera.zoom = this.camera.right / 35.0;
  var container = document.getElementById(element_id);
  if (container === null) { // for testing
    return;
  }
  container.appendChild(this.renderer.domElement);
  if (window.Stats) {
    this.stats = new window.Stats();
    container.appendChild(this.stats.dom);
  }

  window.addEventListener('resize', this.resize.bind(this));
  window.addEventListener('keydown', this.keydown.bind(this));
  window.addEventListener('contextmenu', function (e) { e.preventDefault(); });
  window.addEventListener('mousewheel', this.mousewheel.bind(this));
  window.addEventListener('MozMousePixelScroll', this.mousewheel.bind(this));
  window.addEventListener('mousedown', this.mousedown.bind(this));
  window.addEventListener('touchstart', this.touchstart.bind(this));
  window.addEventListener('touchmove', this.touchmove.bind(this));
  window.addEventListener('touchend', this.touchend.bind(this));
  window.addEventListener('touchcancel', this.touchend.bind(this));
  window.addEventListener('dblclick', this.dblclick.bind(this));

  this.controls.update();

  var self = this;

  this.mousemove = function (event) {
    event.preventDefault();
    //event.stopPropagation();
    self.controls.move(relX(event), relY(event));
  };

  this.mouseup = function (event) {
    event.preventDefault();
    event.stopPropagation();
    self.controls.stop(self.active_model_bag);
    document.removeEventListener('mousemove', self.mousemove);
    document.removeEventListener('mouseup', self.mouseup);
    self.redraw_maps();
  };
}

Viewer.prototype.hud = function (text, type) {
  if (typeof document === 'undefined') return;  // for testing on node
  var el = document && document.getElementById('hud');
  if (el) {
    if (this.initial_hud_text === null) {
      this.initial_hud_text = el.textContent;
      this.initial_hud_bg = el.style['background-color'];
    }
    el.textContent = (text !== undefined ? text : this.initial_hud_text);
    el.style['background-color'] = (type !== 'ERR' ? this.initial_hud_bg
                                                   : '#b00');
    if (type === 'ERR') console.log('ERR: ' + text);
  } else {
    console.log('hud: ' + text);
  }
};

Viewer.prototype.toggle_help = function () {
  var el = document.getElementById('help');
  if (el) {
    el.style.display = el.style.display === 'block' ? 'none' : 'block';
  }
};

Viewer.prototype.redraw_center = function () {
  if (this.target.distanceToSquared(this.last_ctr) > 0.0001) {
    this.last_ctr.copy(this.target);
    if (this.mark) {
      this.scene.remove(this.mark);
    }
    this.mark = make_center_cube(0.1, this.target, this.config.colors.center);
    this.scene.add(this.mark);
  }
};

Viewer.prototype.redraw_maps = function (force) {
  this.redraw_center();
  for (var i = 0; i < this.map_bags.length; i++) {
    var map_bag = this.map_bags[i];
    if (force || this.target.distanceToSquared(map_bag.block_ctr) > 0.01) {
      this.redraw_map(map_bag);
    }
  }
};

Viewer.prototype.clear_el_objects = function (map_bag) {
  for (var i = 0; i < map_bag.el_objects.length; i++) {
    this.scene.remove(map_bag.el_objects[i]);
    map_bag.el_objects[i].geometry.dispose();
  }
  map_bag.el_objects = [];
};

Viewer.prototype.clear_atomic_objects = function (model) {
  if (model.atomic_objects) {
    for (var i = 0; i < model.atomic_objects.length; i++) {
      this.scene.remove(model.atomic_objects[i]);
    }
  }
  model.atomic_objects = null;
};

Viewer.prototype.set_atomic_objects = function (model_bag) {
  model_bag.atomic_objects = [];
  switch (model_bag.conf.render_style) {
    case 'lines':
      model_bag.add_bonds();
      break;
    case 'lines+balls':
      var h_scale = this.camera.projectionMatrix.elements[5];
      var ball_size = Math.max(1, 80 * h_scale);
      model_bag.add_bonds(false, ball_size);
      break;
    case 'trace':  // + lines for ligands
      model_bag.add_trace();
      model_bag.add_bonds(true);
      break;
    case 'ribbon':
      model_bag.add_ribbon();
      model_bag.add_bonds(true);
      break;
  }
  for (var i = 0; i < model_bag.atomic_objects.length; i++) {
    this.scene.add(model_bag.atomic_objects[i]);
  }
};

Viewer.prototype.toggle_map_visibility = function (map_bag) {
  if (typeof map_bag === 'number') {
    map_bag = this.map_bags[map_bag];
  }
  map_bag.visible = !map_bag.visible;
  this.redraw_map(map_bag);
};

Viewer.prototype.redraw_map = function (map_bag) {
  this.clear_el_objects(map_bag);
  if (map_bag.visible) {
    map_bag.map.block = null;
    this.add_el_objects(map_bag);
  }
};

Viewer.prototype.toggle_model_visibility = function (model_bag) {
  model_bag = model_bag || this.active_model_bag;
  model_bag.visible = !model_bag.visible;
  this.redraw_model(model_bag);
};

Viewer.prototype.redraw_model = function (model_bag) {
  this.clear_atomic_objects(model_bag);
  if (model_bag.visible) {
    this.set_atomic_objects(model_bag);
  }
};

Viewer.prototype.redraw_models = function () {
  for (var i = 0; i < this.model_bags.length; i++) {
    this.redraw_model(this.model_bags[i]);
  }
};

Viewer.prototype.add_el_objects = function (map_bag) {
  if (!map_bag.visible) return;
  if (!map_bag.map.block) {
    map_bag.block_ctr.copy(this.target);
    map_bag.map.extract_block(this.config.map_radius,
                              [this.target.x, this.target.y, this.target.z]);
  }
  for (var i = 0; i < map_bag.types.length; i++) {
    var mtype = map_bag.types[i];
    var isolevel = (mtype === 'map_neg' ? -1 : 1) * map_bag.isolevel;
    var iso = map_bag.map.isomesh_in_block(isolevel, this.config.map_style);
    var geom = new THREE.BufferGeometry();
    geom.addAttribute('position',
                 new THREE.BufferAttribute(new Float32Array(iso.vertices), 3));
    /* old version - mesh instead of lines
    geom.setIndex(new THREE.BufferAttribute(new Uint32Array(iso.faces), 1));
    var material = new THREE.MeshBasicMaterial({
      color: this.config.colors[mtype],
      wireframe: true,
      wireframeLinewidth: this.config.map_line
    });
    var obj = new THREE.Mesh(geom, material);
    */

    var faces = iso.faces;
    var arr = new Uint32Array(faces.length * 2);
    for (var j = 0; j < faces.length; j += 3) {
      arr[2*j] = faces[j];
      arr[2*j+1] = faces[j+1];
      arr[2*j+2] = faces[j+1];
      arr[2*j+3] = faces[j+2];
      arr[2*j+4] = faces[j+2];
      arr[2*j+5] = faces[j];
    }
    geom.setIndex(new THREE.BufferAttribute(arr, 1));
    var material = new THREE.LineBasicMaterial({
      color: this.config.colors[mtype],
      linewidth: this.config.map_line
    });
    var obj = new THREE.LineSegments(geom, material);

    map_bag.el_objects.push(obj);
    this.scene.add(obj);
  }
};

Viewer.prototype.change_isolevel_by = function (map_idx, delta) {
  if (map_idx >= this.map_bags.length) return;
  var map_bag = this.map_bags[map_idx];
  map_bag.isolevel += delta;
  var abs_level = map_bag.map.abs_level(map_bag.isolevel);
  this.hud('map ' + (map_idx+1) + ' level =  ' + abs_level.toFixed(4) +
           'e/A^3 (' + map_bag.isolevel.toFixed(2) + ' rmsd)');
  //TODO: move slow part into update()
  this.clear_el_objects(map_bag);
  this.add_el_objects(map_bag);
};

Viewer.prototype.change_map_radius = function (delta) {
  var RMIN = 2;
  var RMAX = 40;
  var cf = this.config;
  cf.map_radius = Math.min(Math.max(cf.map_radius + delta, RMIN), RMAX);
  var info = 'map "radius": ' + cf.map_radius;
  if (cf.map_radius === RMAX) info += ' (max)';
  else if (cf.map_radius === RMIN) info += ' (min)';
  this.hud(info);
  this.redraw_maps(true); //TODO: move slow part into update()
};

Viewer.prototype.toggle_cell_box = function () {
  if (this.decor.cell_box) {
    this.scene.remove(this.decor.cell_box);
    this.decor.cell_box = null;
  } else {
    var uc = null;
    if (this.model_bags.length > 0) {
      uc = this.model_bags[0].model.unit_cell;
    }
    // model may not have unit cell
    if (!uc && this.map_bags.length > 0) {
      uc = this.map_bags[0].map.unit_cell;
    }
    if (uc) {
      this.decor.cell_box = make_unitcell_box(uc, this.config.colors.cell_box);
      this.scene.add(this.decor.cell_box);
    }
  }
};

Viewer.prototype.shift_clip = function (away) {
  var eye = this.camera.position.clone().sub(this.target).setLength(1);
  if (!away) {
    eye.negate();
  }
  this.target.add(eye);
  this.camera.position.add(eye);
  this.update_camera();
  this.redraw_maps();
  this.hud('clip shifted by [' + vec3_to_str(eye, 2, ' ') + ']');
};

Viewer.prototype.go_to_nearest_Ca = function () {
  var t = this.target;
  if (this.active_model_bag === null) return;
  var a = this.active_model_bag.model.get_nearest_atom(t.x, t.y, t.z, 'CA');
  if (a) {
    this.hud(a.long_label());
    //this.set_selection(a);
    this.controls.go_to(new THREE.Vector3(a.xyz[0], a.xyz[1], a.xyz[2]),
                        null, null, 30 / auto_speed);
    this.selected_atom = a;
  } else {
    this.hud('no nearby CA');
  }
};

Viewer.prototype.redraw_all = function () {
  this.scene.fog.color = this.config.colors.bg;
  if (this.renderer) this.renderer.setClearColor(this.config.colors.bg, 1);
  this.redraw_models();
  this.redraw_maps(true);
};

function next(elem, arr) {
  return arr[(arr.indexOf(elem) + 1) % arr.length];
}

function vec3_to_str(vec, n, sep) {
  return vec.x.toFixed(n) + sep + vec.y.toFixed(n) + sep + vec.z.toFixed(n);
}

Viewer.prototype.keydown = function (evt) {  // eslint-disable-line complexity
  var key = evt.keyCode;
  switch (key) {
    case 84:  // t
      this.config.render_style = next(this.config.render_style, RENDER_STYLES);
      this.hud('rendering as ' + this.config.render_style);
      this.redraw_models();
      break;
    case 67:  // c
      if (evt.shiftKey) {
        set_colors(next(this.config.colors.name, Object.keys(ColorSchemes)),
                   this.config.colors);
        this.hud('color scheme: ' + this.config.colors.name);
        this.redraw_all();
      } else { // color-by
        this.config.color_aim = next(this.config.color_aim, COLOR_AIMS);
        this.hud('coloring by ' + this.config.color_aim);
        this.redraw_models();
      }
      break;
    case 87:  // w
      this.config.map_style = next(this.config.map_style, MAP_STYLES);
      this.hud('map style: ' + this.config.map_style);
      this.redraw_maps(true);
      break;
    case 89:  // y
      this.config.hydrogens = !this.config.hydrogens;
      this.hud((this.config.hydrogens ? 'show' : 'hide') +
               ' hydrogens (if any)');
      this.redraw_models();
      break;
    case 107:  // add
    case 61:  // equals/firefox
    case 187:  // equal sign
      this.change_isolevel_by(evt.shiftKey ? 1 : 0, 0.1);
      break;
    case 109:  // subtract
    case 173:  // minus/firefox
    case 189:  // dash
      this.change_isolevel_by(evt.shiftKey ? 1 : 0, -0.1);
      break;
    case 219:  // [
      this.change_map_radius(-2);
      break;
    case 221:  // ]
      this.change_map_radius(2);
      break;
    case 68:  // d
    case 70:  // f
      this.controls.change_slab_width(key === 68 ? -0.1 : +0.1);
      this.update_camera();
      this.hud('clip width: ' + (this.camera.far-this.camera.near).toFixed(1));
      break;
    case 77:  // m
    case 78:  // n
      this.camera.zoom *= (key === 77 ? 1.03 : (1 / 1.03));
      this.update_camera();
      this.hud('zoom: ' + this.camera.zoom.toFixed(2));
      break;
    case 80:  // p
      if (evt.shiftKey) {
        window.location.hash = '#xyz=' + vec3_to_str(this.target, 1, ',') +
          '&eye=' + vec3_to_str(this.camera.position, 1, ',') +
          '&zoom=' + this.camera.zoom.toFixed(0);
        this.hud('copy URL from the location bar');
      } else {
        this.go_to_nearest_Ca();
      }
      break;
    case 51:  // 3
    case 99:  // numpad 3
      this.shift_clip(true);
      break;
    case 108:  // numpad period (Linux)
    case 110:  // decimal point (Mac)
      this.shift_clip(false);
      break;
    case 85:  // u
      this.hud('toggled unit cell box');
      this.toggle_cell_box();
      break;
    case 73:  // i
      this.hud('toggled camera movement');
      this.controls.toggle_state(STATE.AUTO_ROTATE, {rock: evt.shiftKey});
      break;
    case 82:  // r
      if (evt.shiftKey) {
        this.hud('redraw!');
        this.redraw_all();
      } else {
        this.hud('model recentered');
        this.recenter();
      }
      break;
    case 72:  // h
      this.toggle_help();
      break;
    case 36: // Home
    case 35: // End
      this.config.bond_line += (key === 36 ? 0.2 : -0.2);
      this.config.bond_line = Math.max(this.config.bond_line, 0.1);
      this.resize(); // overkill
      this.hud('bond width: ' + this.config.line_width.toFixed(1));
      break;
    case 16: // shift
    case 17: // ctrl
    case 18: // alt
    case 225: // altgr
      break;
    case 32: // Space
      this.center_next_residue(evt.shiftKey);
      break;
    default:
      this.hud('Nothing here. Press H for help.');
      break;
  }
};

Viewer.prototype.mousedown = function (event) {
  event.preventDefault();
  event.stopPropagation();
  var state = STATE.NONE;
  if (event.button === 1 || (event.button === 0 && event.ctrlKey)) {
    state = STATE.PAN;
  } else if (event.button === 0) {
    // in Coot shift+Left is labeling atoms like dblclick, + rotation
    if (event.shiftKey) {
      this.dblclick(event);
    }
    state = STATE.ROTATE;
  } else if (event.button === 2) {
    if (event.ctrlKey) {
      state = event.shiftKey ? STATE.ROLL : STATE.SLAB;
    } else {
      state = STATE.ZOOM;
    }
  }
  this.controls.start(state, relX(event), relY(event));
  document.addEventListener('mousemove', this.mousemove);
  document.addEventListener('mouseup', this.mouseup);
};

Viewer.prototype.dblclick = function (event) {
  if (event.button !== 0) return;
  if (this.decor.selection) {
    this.scene.remove(this.decor.selection);
    this.decor.selection = null;
  }
  var mouse = new THREE.Vector2(relX(event), relY(event));
  var atom;
  if (this.active_model_bag !== null) {
    atom = this.active_model_bag.pick_atom(get_raycaster(mouse, this.camera));
  }
  if (atom) {
    this.hud(atom.long_label());
    this.set_selection(atom);
  } else {
    this.hud();
  }
};

Viewer.prototype.set_selection = function (atom) {
  var geometry = new THREE.Geometry();
  geometry.vertices.push(new THREE.Vector3(atom.xyz[0], atom.xyz[1],
                                           atom.xyz[2]));
  var color = this.config.colors[atom.element] || this.config.colors.def;
  var material = new THREE.PointsMaterial({size: 3, color: color});
  this.decor.selection = new THREE.Points(geometry, material);
  this.scene.add(this.decor.selection);
};

// for two-finger touch events
function touch_info(evt) {
  var touches = evt.touches;
  var dx = touches[0].pageX - touches[1].pageX;
  var dy = touches[0].pageY - touches[1].pageY;
  return {pageX: (touches[0].pageX + touches[1].pageX) / 2,
          pageY: (touches[0].pageY + touches[1].pageY) / 2,
          dist: Math.sqrt(dx * dx + dy * dy)};
}

Viewer.prototype.touchstart = function (event) {
  var touches = event.touches;
  if (touches.length === 1) {
    this.controls.start(STATE.ROTATE, relX(touches[0]), relY(touches[0]));
  } else { // for now using only two touches
    var info = touch_info(event);
    this.controls.start(STATE.PAN_ZOOM, relX(info), relY(info), info.dist);
  }
};

Viewer.prototype.touchmove = function (event) {
  event.preventDefault();
  event.stopPropagation();
  var touches = event.touches;
  if (touches.length === 1) {
    this.controls.move(relX(touches[0]), relY(touches[0]));
  } else { // for now using only two touches
    var info = touch_info(event);
    this.controls.move(relX(info), relY(info), info.dist);
  }
};

Viewer.prototype.touchend = function (/*event*/) {
  this.controls.stop();
  this.redraw_maps();
};

Viewer.prototype.mousewheel = function (evt) {
  evt.preventDefault();
  evt.stopPropagation();
  // evt.wheelDelta for WebKit, evt.detail for Firefox
  var delta = evt.wheelDelta ? evt.wheelDelta / 2000
                             : (evt.detail || 0) / -1000;
  this.change_isolevel_by(evt.shiftKey ? 1 : 0, delta);
};

Viewer.prototype.resize = function (/*evt*/) {
  var width = window.innerWidth;
  var height = window.innerHeight;
  this.camera.left = -width;
  this.camera.right = width;
  this.camera.top = height;
  this.camera.bottom = -height;
  this.camera.updateProjectionMatrix();
  this.renderer.setSize(width, height);
  var line_width = this.config.bond_line * height / 700;
  if (line_width !== this.config.line_width) {
    this.config.line_width = line_width;
    this.redraw_models();
  }
};

// makes sense only for full-window viewer
function parse_fragment() {
  var ret = {};
  if (typeof window === 'undefined') return ret;
  var params = window.location.hash.substr(1).split('&');
  for (var i = 0; i < params.length; i++) {
    var kv = params[i].split('=');
    var val = kv[1];
    if (kv[0] === 'xyz' || kv[0] === 'eye') {
      val = val.split(',').map(Number);
    } else if (kv[0] === 'zoom') {
      val = Number(val);
    }
    ret[kv[0]] = val;
  }
  return ret;
}

// If xyz set recenter on it looking toward the model center.
// Otherwise recenter on the model center looking along the z axis.
Viewer.prototype.recenter = function (xyz, eye, steps) {
  var new_up = null;
  var ctr;
  if (xyz == null || eye == null) {
    ctr = this.active_model_bag.model.get_center();
  }
  if (eye) {
    eye = new THREE.Vector3(eye[0], eye[1], eye[2]);
  }
  if (xyz == null) { // center on the molecule
    if (this.active_model_bag === null) return;
    xyz = new THREE.Vector3(ctr[0], ctr[1], ctr[2]);
    if (!eye) {
      eye = xyz.clone();
      eye.z += 100;
      new_up = THREE.Object3D.DefaultUp; // Vector3(0, 1, 0)
    }
  } else {
    xyz = new THREE.Vector3(xyz[0], xyz[1], xyz[2]);
    if (eye == null && this.active_model_bag !== null) {
      // look toward the center of the molecule
      eye = new THREE.Vector3(ctr[0], ctr[1], ctr[2]);
      eye.sub(xyz).negate().setLength(100); // we store now (eye - xyz)
      new_up = new THREE.Vector3(0, 1, 0).projectOnPlane(eye);
      var len = new_up.length();
      if (len < 0.1) { // the center is in [0,1,0] direction
        new_up.set(1, 0, 0).projectOnPlane(eye);
        len = new_up.length();
      }
      new_up.divideScalar(len); // normalizes
      eye.add(xyz);
    }
  }
  this.controls.go_to(xyz, eye, new_up, steps);
};

Viewer.prototype.center_next_residue = function (back) {
  if (!this.active_model_bag) return;
  var a = this.active_model_bag.model.next_residue(this.selected_atom, back);
  if (a) {
    this.hud('-> ' + a.long_label());
    this.controls.go_to(new THREE.Vector3(a.xyz[0], a.xyz[1], a.xyz[2]),
                        null, null, 30 / auto_speed);
    this.selected_atom = a;
  }
};

Viewer.prototype.update_camera = function () {
  var dxyz = this.camera.position.distanceTo(this.target);
  // the far plane is more distant from the target than the near plane (3:1)
  var w = 0.25 * this.controls.slab_width() / this.camera.zoom;
  this.camera.near = dxyz * (1 - w);
  this.camera.far = dxyz * (1 + 3 * w);
  //this.light.position.copy(this.camera.position);
  var h_scale = this.camera.projectionMatrix.elements[5];
  this.camera.updateProjectionMatrix();
  // temporary hack - scaling balls
  if (h_scale !== this.camera.projectionMatrix.elements[5]) {
    var ball_size = Math.max(1, 80 * this.camera.projectionMatrix.elements[5]);
    for (var i = 0; i < this.model_bags.length; i++) {
      var obj = this.model_bags[i].atomic_objects;
      if (obj.length === 2 && obj[1].material.size) {
        obj[1].material.size = ball_size;
      }
    }
  }
};

Viewer.prototype.render = function render() {
  if (this.controls.update()) {
    this.update_camera();
  }
  if (!this.controls.is_going()) {
    this.redraw_maps();
  }
  this.renderer.render(this.scene, this.camera);
  if (this.nav) {
    this.nav.renderer.render(this.nav.scene, this.camera);
  }
  if (true) { // TODO
    window.requestAnimationFrame(render.bind(this));
  }
  if (this.stats) {
    this.stats.update();
  }
};

Viewer.prototype.set_model = function (model) {
  var model_bag = new ModelBag(model, this.config);
  this.model_bags.push(model_bag);
  this.set_atomic_objects(model_bag);
  this.active_model_bag = model_bag;
};

Viewer.prototype.add_map = function (map, is_diff_map) {
  //map.show_debug_info();
  var map_bag = new MapBag(map, is_diff_map);
  this.map_bags.push(map_bag);
  this.add_el_objects(map_bag);
};

Viewer.prototype.load_file = function (url, binary, callback) {
  var req = new XMLHttpRequest();
  req.open('GET', url, true);
  if (binary) {
    req.responseType = 'arraybuffer';
  } else {
    // http://stackoverflow.com/questions/7374911/
    req.overrideMimeType('text/plain');
  }
  var self = this;
  req.onreadystatechange = function () {
    if (req.readyState === 4) {
      // chrome --allow-file-access-from-files gives status 0
      if (req.status === 200 || (req.status === 0 && req.response !== null)) {
        try {
          callback(req);
        } catch (e) {
          self.hud('Error: ' + e.message + '\nin ' + url, 'ERR');
        }
      } else {
        self.hud('Failed to fetch ' + url, 'ERR');
      }
    }
  };
  req.send(null);
};

// Load molecular model from PDB file and centers the view
Viewer.prototype.load_pdb = function (url, options) {
  options = options || {};
  var self = this;
  this.load_file(url, false, function (req) {
    var model = new Model();
    model.from_pdb(req.responseText);
    self.set_model(model);
    var frag = parse_fragment();
    if (frag.zoom) self.camera.zoom = frag.zoom;
    self.recenter(options.center || frag.xyz, frag.eye, 1);
    if (options.callback) options.callback();
  });
};

Viewer.prototype.load_map = function (url, is_diff_map, filetype, callback) {
  if (filetype !== 'ccp4' && filetype !== 'dsn6') {
    throw Error('Unknown map filetype.');
  }
  var self = this;
  this.load_file(url, true, function (req) {
    var map = new ElMap();
    if (filetype === 'ccp4') map.from_ccp4(req.response);
    else /* === 'dsn6'*/ map.from_dsn6(req.response);
    self.add_map(map, is_diff_map);
    if (callback) callback();
  });
};

// Load a normal map and a difference map.
// To show the first map ASAP we do not download both maps in parallel.
Viewer.prototype.load_ccp4_maps = function (url1, url2, callback) {
  var self = this;
  this.load_map(url1, false, 'ccp4', function () {
    self.load_map(url2, true, 'ccp4', callback);
  });
};

// TODO: navigation window like in gimp and mifit
/*
Viewer.prototype.show_nav = function (inset_id) {
  var inset = document.getElementById(inset_id);
  if (!inset) return;
  inset.style.display = 'block';
  var nav = {};
  nav.renderer = new THREE.WebGLRenderer();
  nav.renderer.setClearColor(0x555555, 1);
  nav.renderer.setSize(200, 200);
  inset.appendChild(nav.renderer.domElement);
  //nav.scene = new THREE.Scene();
  nav.scene = this.scene;
  //var light = new THREE.AmbientLight(0xffffff);
  //nav.scene.add(light);
  this.nav = nav;
};
*/

Viewer.ColorSchemes = ColorSchemes;
Viewer.auto_speed = auto_speed;

return Viewer;
})();

if (typeof module !== 'undefined') module.exports = {
  VERSION: UGLYMOL_VERSION,
  UnitCell: UnitCell,
  Model: Model,
  ElMap: ElMap,
  isosurface: isosurface,
  Viewer: Viewer
};
