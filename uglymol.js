/*!
 * UglyMol v0.7.0. Macromolecular Viewer for Crystallographers.
 * Copyright 2014 Nat Echols
 * Copyright 2016 Diamond Light Source Ltd
 * Copyright 2016 Marcin Wojdyr
 * Released under the MIT License.
 */
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
typeof define === 'function' && define.amd ? define(['exports'], factory) :
(global = global || self, factory(global.UM = {}));
}(this, (function (exports) { 'use strict';

var VERSION = exports.VERSION = '0.7.0';


// @flow

var UnitCell = function UnitCell(a /*:number*/, b /*:number*/, c /*:number*/,
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
  var cos_alpha_star_sin_beta = (cos_beta * cos_gamma - cos_alpha) /
                                  sin_gamma;
  var cos_alpha_star = cos_alpha_star_sin_beta / sin_beta;
  var s1rca2 = Math.sqrt(1.0 - cos_alpha_star * cos_alpha_star);
  // The orthogonalization matrix we use is described in ITfC B p.262:
  // "An alternative mode of orthogonalization, used by the Protein
  // Data Bank and most programs, is to align the a1 axis of the unit
  // cell with the Cartesian X_1 axis, and to align the a*_3 axis with the
  // Cartesian X_3 axis."
  //
  // Zeros in the matrices below are kept to make matrix multiplication
  // faster: they make extract_block() 2x (!) faster on V8 4.5.103,
  // no difference on FF 50.
  /* eslint-disable no-multi-spaces, comma-spacing */
  this.orth = [a, b * cos_gamma,c * cos_beta,
               0.0, b * sin_gamma, -c * cos_alpha_star_sin_beta,
               0.0, 0.0        ,c * sin_beta * s1rca2];
  // based on xtal.js which is based on cctbx.uctbx
  this.frac = [
    1.0 / a,
    -cos_gamma / (sin_gamma * a),
    -(cos_gamma * cos_alpha_star_sin_beta + cos_beta * sin_gamma) /
        (sin_beta * s1rca2 * sin_gamma * a),
    0.0,
    1.0 / (sin_gamma * b),
    cos_alpha_star / (s1rca2 * sin_gamma * b),
    0.0,
    0.0,
    1.0 / (sin_beta * s1rca2 * c) ];
};

UnitCell.prototype.fractionalize = function fractionalize (xyz /*:[number,number,number]*/) {
  return multiply(xyz, this.frac);
};

UnitCell.prototype.orthogonalize = function orthogonalize (xyz /*:[number,number,number]*/) {
  return multiply(xyz, this.orth);
};

// This function is only used with matrices frac and orth, which have 3 zeros.
// We skip these elements, but it doesn't affect performance (on FF50 and V8).
function multiply(xyz, mat) {
  /* eslint-disable indent */
  return [mat[0] * xyz[0]  + mat[1] * xyz[1]  + mat[2] * xyz[2],
        /*mat[3] * xyz[0]*/+ mat[4] * xyz[1]  + mat[5] * xyz[2],
        /*mat[6] * xyz[0]  + mat[7] * xyz[1]*/+ mat[8] * xyz[2]];
}

// @flow

/*::
 type Num3 = [number, number, number];
 */

var AMINO_ACIDS = [
  'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
  'LYS', 'MET', 'MSE', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'UNK' ];
var NUCLEIC_ACIDS = [
  'DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U', 'rA', 'rC', 'rG', 'rU',
  'Ar', 'Cr', 'Gr', 'Ur' ];

var NOT_LIGANDS = ['HOH'].concat(AMINO_ACIDS, NUCLEIC_ACIDS);

function modelsFromPDB(pdb_string/*:string*/) {
  var models = [new Model()];
  var pdb_tail = models[0].from_pdb(pdb_string.split('\n'));
  while (pdb_tail != null) {
    var model = new Model();
    pdb_tail = model.from_pdb(pdb_tail);
    if (model.atoms.length > 0) { models.push(model); }
  }
  return models;
}

var Model = function Model() {
  this.atoms = [];
  this.unit_cell = null;
  this.space_group = null;
  this.has_hydrogens = false;
  this.lower_bound = [0, 0, 0];
  this.upper_bound = [0, 0, 0];
  this.residue_map = null;
  this.cubes = null;
};

Model.prototype.from_pdb = function from_pdb (pdb_lines /*:string[]*/) /*:?string[]*/ {
  var chain_index = 0;// will be ++'ed for the first atom
  var last_chain = null;
  var atom_i_seq = 0;
  var continuation = null;
  for (var i = 0; i < pdb_lines.length; i++) {
    var line = pdb_lines[i];
    var rec_type = line.substring(0, 6).toUpperCase();
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
      this.atoms.push(new_atom);
    } else if (rec_type === 'ANISOU') ; else if (rec_type === 'CRYST1') {
      var a = parseFloat(line.substring(6, 15));
      var b = parseFloat(line.substring(15, 24));
      var c = parseFloat(line.substring(24, 33));
      var alpha = parseFloat(line.substring(33, 40));
      var beta = parseFloat(line.substring(40, 47));
      var gamma = parseFloat(line.substring(47, 54));
      //const sg_symbol = line.substring(55, 66);
      this.unit_cell = new UnitCell(a, b, c, alpha, beta, gamma);
    } else if (rec_type.substring(0, 3) === 'TER') {
      last_chain = null;
    } else if (rec_type === 'ENDMDL') {
      for (; i < pdb_lines.length; i++) {
        if (pdb_lines[i].substring(0, 6).toUpperCase() === 'MODEL ') {
          continuation = pdb_lines.slice(i);
          break;
        }
      }
      break;
    }
  }
  if (this.atoms.length === 0) { throw Error('No atom records found.'); }
  this.calculate_bounds();
  this.calculate_connectivity();
  return continuation;
};

Model.prototype.calculate_bounds = function calculate_bounds () {
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

Model.prototype.next_residue = function next_residue (atom /*:?Atom*/, backward /*:?boolean*/) {
  var len = this.atoms.length;
  var start = (atom ? atom.i_seq : 0) + len;// +len to avoid idx<0 below
  for (var i = (atom ? 1 : 0); i < len; i++) {
    var idx = (start + (backward ? -i : i)) % len;
    var a = this.atoms[idx];
    if (!a.is_main_conformer()) { continue; }
    if ((a.name === 'CA' && a.element === 'C') || a.name === 'P') {
      return a;
    }
  }
};

Model.prototype.extract_trace = function extract_trace () {
  var segments = [];
  var current_segment = [];
  var last_atom = null;
  for (var i = 0; i < this.atoms.length; i++) {
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

Model.prototype.get_residues = function get_residues () {
  if (this.residue_map != null) { return this.residue_map; }
  var residues = {};
  for (var i = 0; i < this.atoms.length; i++) {
    var atom = this.atoms[i];
    var resid = atom.resid();
    var reslist = residues[resid];
    if (reslist === undefined) {
      residues[resid] = [atom];
    } else {
      reslist.push(atom);
    }
  }
  this.residue_map = residues;
  return residues;
};

// tangent vector to the ribbon representation
Model.prototype.calculate_tangent_vector = function calculate_tangent_vector (residue /*:Atom[]*/) {
  var a1 = null;
  var a2 = null;
  // it may be too simplistic
  var peptide = (residue[0].resname.length === 3);
  var name1 = peptide ? 'C' : 'C2\'';
  var name2 = peptide ? 'O' : 'O4\'';
  for (var i = 0; i < residue.length; i++) {
    var atom = residue[i];
    if (!atom.is_main_conformer()) { continue; }
    if (atom.name === name1) {
      a1 = atom.xyz;
    } else if (atom.name === name2) {
      a2 = atom.xyz;
    }
  }
  if (a1 === null || a2 === null) { return [0, 0, 1]; } // arbitrary value
  var d = [a1[0]-a2[0], a1[1]-a2[1], a1[2]-a2[2]];
  var len = Math.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
  return [d[0]/len, d[1]/len, d[2]/len];
};

Model.prototype.get_center = function get_center () {
  var xsum = 0, ysum = 0, zsum = 0;// eslint-disable-line
  var n_atoms = this.atoms.length;
  for (var i = 0; i < n_atoms; i++) {
    var xyz = this.atoms[i].xyz;
    xsum += xyz[0];
    ysum += xyz[1];
    zsum += xyz[2];
  }
  return [xsum / n_atoms, ysum / n_atoms, zsum / n_atoms];
};

Model.prototype.calculate_connectivity = function calculate_connectivity () {
  var atoms = this.atoms;
  var cubes = new Cubicles(atoms, 3.0, this.lower_bound, this.upper_bound);
  //let cnt = 0;
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
          //cnt++;
        }
      }
    }
  }
  //console.log(atoms.length + ' atoms, ' + cnt + ' bonds.');
  this.cubes = cubes;
};

Model.prototype.get_nearest_atom = function get_nearest_atom (x /*:number*/, y /*:number*/, z /*:number*/,
                 atom_name /*:?string*/) {
  var cubes = this.cubes;
  if (cubes == null) { throw Error('Missing Cubicles'); }
  var box_id = cubes.find_box_id(x, y, z);
  var indices = cubes.get_nearby_atoms(box_id);
  var nearest = null;
  var min_d2 = Infinity;
  for (var i = 0; i < indices.length; i++) {
    var atom = this.atoms[indices[i]];
    if (atom_name != null && atom_name !== atom.name) { continue; }
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

// Single atom and associated labels
var Atom = function Atom() {
  this.hetero = false;
  this.name = '';
  this.altloc = '';
  this.resname = '';
  this.chain = '';
  this.chain_index = -1;
  this.resseq = -1;
  this.icode = null;
  this.xyz = [0, 0, 0];
  this.occ = 1.0;
  this.b = 0;
  this.element = '';
  this.i_seq = -1;
  this.is_ligand = null;
  this.bonds = [];
};

// http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
Atom.prototype.from_pdb_line = function from_pdb_line (pdb_line /*:string*/) {
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
  //if (pdb_line.length >= 80) {
  //this.charge = pdb_line.substring(78, 80).trim();
  //}
  this.is_ligand = (NOT_LIGANDS.indexOf(this.resname) === -1);
};

Atom.prototype.b_as_u = function b_as_u () {
  // B = 8 * pi^2 * u^2
  return Math.sqrt(this.b / (8 * 3.14159 * 3.14159));
};

Atom.prototype.distance_sq = function distance_sq (other /*:Atom*/) {
  var dx = this.xyz[0] - other.xyz[0];
  var dy = this.xyz[1] - other.xyz[1];
  var dz = this.xyz[2] - other.xyz[2];
  return dx*dx + dy*dy + dz*dz;
};

Atom.prototype.distance = function distance (other /*:Atom*/) {
  return Math.sqrt(this.distance_sq(other));
};

Atom.prototype.midpoint = function midpoint (other /*:Atom*/) {
  return [(this.xyz[0] + other.xyz[0]) / 2,
          (this.xyz[1] + other.xyz[1]) / 2,
          (this.xyz[2] + other.xyz[2]) / 2];
};

Atom.prototype.is_hydrogen = function is_hydrogen () {
  return this.element === 'H' || this.element === 'D';
};

Atom.prototype.is_ion = function is_ion () {
  return this.element === this.resname;
};

Atom.prototype.is_water = function is_water () {
  return this.resname === 'HOH';
};

Atom.prototype.is_same_residue = function is_same_residue (other /*:Atom*/, ignore_altloc /*:?boolean*/) {
  return other.resseq === this.resseq && other.icode === this.icode &&
         other.chain === this.chain && other.resname === this.resname &&
         (ignore_altloc || other.altloc === this.altloc);
};

Atom.prototype.is_same_conformer = function is_same_conformer (other /*:Atom*/) {
  return this.altloc === '' || other.altloc === '' ||
         this.altloc === other.altloc;
};

Atom.prototype.is_main_conformer = function is_main_conformer () {
  return this.altloc === '' || this.altloc === 'A';
};

Atom.prototype.bond_radius = function bond_radius () { // rather crude
  if (this.element === 'H') { return 1.3; }
  if (this.element === 'S' || this.element === 'P') { return 2.43; }
  return 1.99;
};

Atom.prototype.is_bonded_to = function is_bonded_to (other /*:Atom*/) {
  var MAX_DIST = 2.2 * 2.2;
  if (!this.is_same_conformer(other)) { return false; }
  var dxyz2 = this.distance_sq(other);
  if (dxyz2 > MAX_DIST) { return false; }
  if (this.element === 'H' && other.element === 'H') { return false; }
  return dxyz2 <= this.bond_radius() * other.bond_radius();
};

Atom.prototype.resid = function resid () {
  return this.resseq + '/' + this.chain;
};

Atom.prototype.long_label = function long_label () {
  var a = this;
  return a.name + ' /' + a.resseq + ' ' + a.resname + '/' + a.chain +
         ' - occ: ' + a.occ.toFixed(2) + ' bf: ' + a.b.toFixed(2) +
         ' ele: ' + a.element + ' pos: (' + a.xyz[0].toFixed(2) + ',' +
         a.xyz[1].toFixed(2) + ',' + a.xyz[2].toFixed(2) + ')';
};

Atom.prototype.short_label = function short_label () {
  var a = this;
  return a.name + ' /' + a.resseq + ' ' + a.resname + '/' + a.chain;
};
/*:: export type AtomT = Atom; */


// Partition atoms into boxes for quick neighbor searching.
var Cubicles = function Cubicles(atoms/*:Atom[]*/, box_length/*:number*/,
            lower_bound/*:Num3*/, upper_bound/*:Num3*/) {
  this.boxes = [];
  this.box_length = box_length;
  this.lower_bound = lower_bound;
  this.upper_bound = upper_bound;
  this.xdim = Math.ceil((upper_bound[0] - lower_bound[0]) / box_length);
  this.ydim = Math.ceil((upper_bound[1] - lower_bound[1]) / box_length);
  this.zdim = Math.ceil((upper_bound[2] - lower_bound[2]) / box_length);
  //console.log("Cubicles: " + this.xdim + "x" + this.ydim + "x" + this.zdim);
  var nxyz = this.xdim * this.ydim * this.zdim;
  for (var j = 0; j < nxyz; j++) {
    this.boxes.push([]);
  }
  for (var i = 0; i < atoms.length; i++) {
    var xyz = atoms[i].xyz;
    var box_id = this.find_box_id(xyz[0], xyz[1], xyz[2]);
    if (box_id === null) {
      throw Error('wrong cubicle');
    }
    this.boxes[box_id].push(i);
  }
};

Cubicles.prototype.find_box_id = function find_box_id (x/*:number*/, y/*:number*/, z/*:number*/) {
  var xstep = Math.floor((x - this.lower_bound[0]) / this.box_length);
  var ystep = Math.floor((y - this.lower_bound[1]) / this.box_length);
  var zstep = Math.floor((z - this.lower_bound[2]) / this.box_length);
  var box_id = (zstep * this.ydim + ystep) * this.xdim + xstep;
  if (box_id < 0 || box_id >= this.boxes.length) { throw Error('Ups!'); }
  return box_id;
};

Cubicles.prototype.get_nearby_atoms = function get_nearby_atoms (box_id/*:number*/) {
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

// @flow
/*:: type Num3 = [number, number, number] */

var Block = function Block() {
  this._points = null;
  this._values = null;
  this._size = [0, 0, 0];
};

Block.prototype.set = function set (points /*:Num3[]*/, values/*:number[]*/, size/*:Num3*/) {
  if (size[0] <= 0 || size[1] <= 0 || size[2] <= 0) {
    throw Error('Grid dimensions are zero along at least one edge');
  }
  var len = size[0] * size[1] * size[2];
  if (values.length !== len || points.length !== len) {
    throw Error('isosurface: array size mismatch');
  }

  this._points = points;
  this._values = values;
  this._size = size;
};

Block.prototype.clear = function clear () {
  this._points = null;
  this._values = null;
};

Block.prototype.empty = function empty () /*:boolean*/ {
  return this._values === null;
};

Block.prototype.isosurface = function isosurface (isolevel /*: number*/, method /*: string*/) {
  //if (method === 'marching tetrahedra') {
  //return marchingTetrahedra(block, isolevel);
  //}
  return marchingCubes(this._size, this._values, this._points,
                       isolevel, method);
};

/* eslint comma-spacing: 0, no-multi-spaces: 0 */

var edgeTable = new Int32Array([
  0x0  , 0x0  , 0x202, 0x302, 0x406, 0x406, 0x604, 0x704,
  0x804, 0x805, 0xa06, 0xa06, 0xc0a, 0xd03, 0xe08, 0xf00,
  0x90 , 0x98 , 0x292, 0x292, 0x496, 0x49e, 0x694, 0x694,
  0x894, 0x894, 0xa96, 0xa96, 0xc9a, 0xc92, 0xe91, 0xe90,
  0x230, 0x230, 0x33 , 0x13a, 0x636, 0x636, 0x434, 0x43c,
  0xa34, 0xa35, 0x837, 0x936, 0xe3a, 0xf32, 0xc31, 0xd30,
  0x2a0, 0x2a8, 0xa3 , 0xaa , 0x6a6, 0x6af, 0x5a4, 0x4ac,
  0xaa4, 0xaa4, 0x9a6, 0x8a6, 0xfaa, 0xea3, 0xca1, 0xca0,
  0x460, 0x460, 0x662, 0x762, 0x66 , 0x66 , 0x265, 0x364,
  0xc64, 0xc65, 0xe66, 0xe66, 0x86a, 0x863, 0xa69, 0xa60,
  0x4f0, 0x4f8, 0x6f2, 0x6f2, 0xf6 , 0xfe , 0x2f5, 0x2fc,
  0xcf4, 0xcf4, 0xef6, 0xef6, 0x8fa, 0x8f3, 0xaf9, 0xaf0,
  0x650, 0x650, 0x453, 0x552, 0x256, 0x256, 0x54 , 0x154,
  0xe54, 0xf54, 0xc57, 0xd56, 0xa5a, 0xb52, 0x859, 0x950,
  0x7c0, 0x6c1, 0x5c2, 0x4c2, 0x3c6, 0x2ce, 0xc5 , 0xc4 ,
  0xfc4, 0xec5, 0xdc6, 0xcc6, 0xbca, 0xac2, 0x8c1, 0x8c0,
  0x8c0, 0x8c0, 0xac2, 0xbc2, 0xcc6, 0xcc6, 0xec4, 0xfcc,
  0xc4 , 0xc5 , 0x2c6, 0x3c6, 0x4c2, 0x5c2, 0x6c1, 0x7c0,
  0x950, 0x859, 0xb52, 0xa5a, 0xd56, 0xc57, 0xe54, 0xe5c,
  0x154, 0x54 , 0x25e, 0x256, 0x552, 0x453, 0x658, 0x650,
  0xaf0, 0xaf0, 0x8f3, 0x8fa, 0xef6, 0xef6, 0xcf4, 0xcfc,
  0x2f4, 0x3f5, 0xff , 0x1f6, 0x6f2, 0x6f3, 0x4f9, 0x5f0,
  0xa60, 0xa69, 0x863, 0x86a, 0xe66, 0xe67, 0xd65, 0xc6c,
  0x364, 0x265, 0x166, 0x66 , 0x76a, 0x663, 0x460, 0x460,
  0xca0, 0xca0, 0xea2, 0xfa2, 0x8a6, 0x8a6, 0xaa4, 0xba4,
  0x4ac, 0x5a4, 0x6ae, 0x7a6, 0xaa , 0xa3 , 0x2a8, 0x2a0,
  0xd30, 0xc31, 0xf32, 0xe3a, 0x936, 0x837, 0xb35, 0xa34,
  0x43c, 0x434, 0x73e, 0x636, 0x13a, 0x33 , 0x339, 0x230,
  0xe90, 0xe90, 0xc92, 0xc9a, 0xa96, 0xa96, 0x894, 0x89c,
  0x694, 0x695, 0x49f, 0x496, 0x292, 0x392, 0x98 , 0x90 ,
  0xf00, 0xe08, 0xd03, 0xc0a, 0xa06, 0xa0e, 0x805, 0x804,
  0x704, 0x604, 0x506, 0x406, 0x302, 0x202, 0x0  , 0x0]);


// generated from classical triTable by tools/isolut.py
var segTable = [
  [],
  [],
  [1, 9],
  [1, 8, 1, 9],
  [2, 10, 10, 1],
  [2, 10, 10, 1],
  [9, 2, 2, 10, 10, 9],
  [2, 8, 2, 10, 10, 8, 10, 9],
  [11, 2],
  [0, 11, 11, 2],
  [1, 9, 11, 2],
  [1, 11, 11, 2, 1, 9, 9, 11],
  [3, 10, 10, 1, 11, 10],
  [0, 10, 10, 1, 8, 10, 11, 10],
  [3, 9, 11, 9, 11, 10, 10, 9],
  [8, 10, 10, 9, 11, 10],
  [4, 7],
  [4, 3, 4, 7],
  [1, 9, 4, 7],
  [4, 1, 1, 9, 4, 7, 7, 1],
  [2, 10, 10, 1, 4, 7],
  [3, 4, 4, 7, 2, 10, 10, 1],
  [9, 2, 2, 10, 10, 9, 4, 7],
  [2, 10, 10, 9, 9, 2, 9, 7, 7, 2, 4, 7],
  [4, 7, 11, 2],
  [11, 4, 4, 7, 11, 2, 2, 4],
  [1, 9, 4, 7, 11, 2],
  [4, 7, 11, 4, 11, 9, 11, 2, 2, 9, 1, 9],
  [3, 10, 10, 1, 11, 10, 4, 7],
  [1, 11, 11, 10, 10, 1, 1, 4, 4, 11, 4, 7],
  [4, 7, 0, 11, 11, 9, 11, 10, 10, 9],
  [4, 7, 11, 4, 11, 9, 11, 10, 10, 9],
  [9, 5, 5, 4],
  [9, 5, 5, 4],
  [0, 5, 5, 4, 1, 5],
  [8, 5, 5, 4, 3, 5, 1, 5],
  [2, 10, 10, 1, 9, 5, 5, 4],
  [2, 10, 10, 1, 9, 5, 5, 4],
  [5, 2, 2, 10, 10, 5, 5, 4, 4, 2],
  [2, 10, 10, 5, 5, 2, 5, 3, 5, 4, 4, 3],
  [9, 5, 5, 4, 11, 2],
  [0, 11, 11, 2, 9, 5, 5, 4],
  [0, 5, 5, 4, 1, 5, 11, 2],
  [1, 5, 5, 2, 5, 8, 8, 2, 11, 2, 5, 4],
  [10, 3, 11, 10, 10, 1, 9, 5, 5, 4],
  [9, 5, 5, 4, 8, 1, 8, 10, 10, 1, 11, 10],
  [5, 4, 0, 5, 0, 11, 11, 5, 11, 10, 10, 5],
  [5, 4, 8, 5, 8, 10, 10, 5, 11, 10],
  [9, 7, 5, 7, 9, 5],
  [9, 3, 9, 5, 5, 3, 5, 7],
  [0, 7, 1, 7, 1, 5, 5, 7],
  [1, 5, 5, 3, 5, 7],
  [9, 7, 9, 5, 5, 7, 10, 1, 2, 10],
  [10, 1, 2, 10, 9, 5, 5, 0, 5, 3, 5, 7],
  [2, 8, 2, 5, 5, 8, 5, 7, 10, 5, 2, 10],
  [2, 10, 10, 5, 5, 2, 5, 3, 5, 7],
  [7, 9, 9, 5, 5, 7, 11, 2],
  [9, 5, 5, 7, 7, 9, 7, 2, 2, 9, 11, 2],
  [11, 2, 1, 8, 1, 7, 1, 5, 5, 7],
  [11, 2, 1, 11, 1, 7, 1, 5, 5, 7],
  [9, 5, 5, 8, 5, 7, 10, 1, 3, 10, 11, 10],
  [5, 7, 7, 0, 0, 5, 9, 5, 11, 0, 0, 10, 10, 1, 11, 10],
  [11, 10, 10, 0, 0, 11, 10, 5, 5, 0, 0, 7, 5, 7],
  [11, 10, 10, 5, 5, 11, 5, 7],
  [10, 6, 6, 5, 5, 10],
  [5, 10, 10, 6, 6, 5],
  [1, 9, 5, 10, 10, 6, 6, 5],
  [1, 8, 1, 9, 5, 10, 10, 6, 6, 5],
  [1, 6, 6, 5, 5, 1, 2, 6],
  [1, 6, 6, 5, 5, 1, 2, 6],
  [9, 6, 6, 5, 5, 9, 0, 6, 2, 6],
  [5, 9, 8, 5, 8, 2, 2, 5, 2, 6, 6, 5],
  [11, 2, 10, 6, 6, 5, 5, 10],
  [11, 0, 11, 2, 10, 6, 6, 5, 5, 10],
  [1, 9, 11, 2, 5, 10, 10, 6, 6, 5],
  [5, 10, 10, 6, 6, 5, 1, 9, 9, 2, 9, 11, 11, 2],
  [6, 3, 11, 6, 6, 5, 5, 3, 5, 1],
  [11, 0, 11, 5, 5, 0, 5, 1, 11, 6, 6, 5],
  [11, 6, 6, 3, 6, 0, 6, 5, 5, 0, 5, 9],
  [6, 5, 5, 9, 9, 6, 9, 11, 11, 6],
  [5, 10, 10, 6, 6, 5, 4, 7],
  [4, 3, 4, 7, 6, 5, 5, 10, 10, 6],
  [1, 9, 5, 10, 10, 6, 6, 5, 4, 7],
  [10, 6, 6, 5, 5, 10, 1, 9, 9, 7, 7, 1, 4, 7],
  [6, 1, 2, 6, 6, 5, 5, 1, 4, 7],
  [2, 5, 5, 1, 2, 6, 6, 5, 4, 3, 4, 7],
  [4, 7, 0, 5, 5, 9, 0, 6, 6, 5, 2, 6],
  [3, 9, 9, 7, 4, 7, 2, 9, 5, 9, 9, 6, 6, 5, 2, 6],
  [11, 2, 4, 7, 10, 6, 6, 5, 5, 10],
  [5, 10, 10, 6, 6, 5, 4, 7, 7, 2, 2, 4, 11, 2],
  [1, 9, 4, 7, 11, 2, 5, 10, 10, 6, 6, 5],
  [9, 2, 1, 9, 9, 11, 11, 2, 4, 11, 4, 7, 5, 10, 10, 6, 6, 5],
  [4, 7, 11, 5, 5, 3, 5, 1, 11, 6, 6, 5],
  [5, 1, 1, 11, 11, 5, 11, 6, 6, 5, 0, 11, 11, 4, 4, 7],
  [0, 5, 5, 9, 0, 6, 6, 5, 3, 6, 11, 6, 4, 7],
  [6, 5, 5, 9, 9, 6, 9, 11, 11, 6, 4, 7, 7, 9],
  [10, 4, 9, 10, 6, 4, 10, 6],
  [4, 10, 10, 6, 6, 4, 9, 10],
  [10, 0, 1, 10, 10, 6, 6, 0, 6, 4],
  [1, 8, 1, 6, 6, 8, 6, 4, 1, 10, 10, 6],
  [1, 4, 9, 1, 2, 4, 2, 6, 6, 4],
  [2, 9, 9, 1, 2, 4, 2, 6, 6, 4],
  [2, 4, 2, 6, 6, 4],
  [2, 8, 2, 4, 2, 6, 6, 4],
  [10, 4, 9, 10, 10, 6, 6, 4, 11, 2],
  [8, 2, 11, 2, 9, 10, 10, 4, 10, 6, 6, 4],
  [11, 2, 1, 6, 6, 0, 6, 4, 1, 10, 10, 6],
  [6, 4, 4, 1, 1, 6, 1, 10, 10, 6, 8, 1, 1, 11, 11, 2],
  [9, 6, 6, 4, 9, 3, 3, 6, 9, 1, 11, 6],
  [11, 1, 1, 8, 11, 6, 6, 1, 9, 1, 1, 4, 6, 4],
  [11, 6, 6, 3, 6, 0, 6, 4],
  [6, 4, 8, 6, 11, 6],
  [7, 10, 10, 6, 6, 7, 8, 10, 9, 10],
  [0, 7, 0, 10, 10, 7, 9, 10, 6, 7, 10, 6],
  [10, 6, 6, 7, 7, 10, 1, 10, 7, 1, 8, 1],
  [10, 6, 6, 7, 7, 10, 7, 1, 1, 10],
  [2, 6, 6, 1, 6, 8, 8, 1, 9, 1, 6, 7],
  [2, 6, 6, 9, 9, 2, 9, 1, 6, 7, 7, 9, 9, 3],
  [0, 7, 0, 6, 6, 7, 2, 6],
  [2, 7, 6, 7, 2, 6],
  [11, 2, 10, 6, 6, 8, 8, 10, 9, 10, 6, 7],
  [0, 7, 7, 2, 11, 2, 9, 7, 6, 7, 7, 10, 10, 6, 9, 10],
  [1, 8, 1, 7, 1, 10, 10, 7, 6, 7, 10, 6, 11, 2],
  [11, 2, 1, 11, 1, 7, 10, 6, 6, 1, 1, 10, 6, 7],
  [9, 6, 6, 8, 6, 7, 9, 1, 1, 6, 11, 6, 6, 3],
  [9, 1, 11, 6, 6, 7],
  [0, 7, 0, 6, 6, 7, 11, 0, 11, 6],
  [11, 6, 6, 7],
  [7, 6, 6, 11],
  [7, 6, 6, 11],
  [1, 9, 7, 6, 6, 11],
  [8, 1, 1, 9, 7, 6, 6, 11],
  [10, 1, 2, 10, 6, 11, 7, 6],
  [2, 10, 10, 1, 6, 11, 7, 6],
  [2, 9, 2, 10, 10, 9, 6, 11, 7, 6],
  [6, 11, 7, 6, 2, 10, 10, 3, 10, 8, 10, 9],
  [7, 2, 6, 2, 7, 6],
  [7, 0, 7, 6, 6, 0, 6, 2],
  [2, 7, 7, 6, 6, 2, 1, 9],
  [1, 6, 6, 2, 1, 8, 8, 6, 1, 9, 7, 6],
  [10, 7, 7, 6, 6, 10, 10, 1, 1, 7],
  [10, 7, 7, 6, 6, 10, 1, 7, 10, 1, 1, 8],
  [7, 0, 7, 10, 10, 0, 10, 9, 6, 10, 7, 6],
  [7, 6, 6, 10, 10, 7, 10, 8, 10, 9],
  [6, 8, 4, 6, 6, 11],
  [3, 6, 6, 11, 0, 6, 4, 6],
  [8, 6, 6, 11, 4, 6, 1, 9],
  [4, 6, 6, 9, 6, 3, 3, 9, 1, 9, 6, 11],
  [6, 8, 4, 6, 6, 11, 2, 10, 10, 1],
  [2, 10, 10, 1, 0, 11, 0, 6, 6, 11, 4, 6],
  [4, 11, 4, 6, 6, 11, 2, 9, 2, 10, 10, 9],
  [10, 9, 9, 3, 3, 10, 2, 10, 4, 3, 3, 6, 6, 11, 4, 6],
  [8, 2, 4, 2, 4, 6, 6, 2],
  [4, 2, 4, 6, 6, 2],
  [1, 9, 3, 4, 4, 2, 4, 6, 6, 2],
  [1, 9, 4, 1, 4, 2, 4, 6, 6, 2],
  [8, 1, 8, 6, 6, 1, 4, 6, 6, 10, 10, 1],
  [10, 1, 0, 10, 0, 6, 6, 10, 4, 6],
  [4, 6, 6, 3, 3, 4, 6, 10, 10, 3, 3, 9, 10, 9],
  [10, 9, 4, 10, 6, 10, 4, 6],
  [9, 5, 5, 4, 7, 6, 6, 11],
  [9, 5, 5, 4, 7, 6, 6, 11],
  [5, 0, 1, 5, 5, 4, 7, 6, 6, 11],
  [7, 6, 6, 11, 3, 4, 3, 5, 5, 4, 1, 5],
  [9, 5, 5, 4, 10, 1, 2, 10, 7, 6, 6, 11],
  [6, 11, 7, 6, 2, 10, 10, 1, 9, 5, 5, 4],
  [7, 6, 6, 11, 5, 4, 4, 10, 10, 5, 4, 2, 2, 10],
  [3, 4, 3, 5, 5, 4, 2, 5, 10, 5, 2, 10, 7, 6, 6, 11],
  [7, 2, 7, 6, 6, 2, 5, 4, 9, 5],
  [9, 5, 5, 4, 8, 6, 6, 0, 6, 2, 7, 6],
  [3, 6, 6, 2, 7, 6, 1, 5, 5, 0, 5, 4],
  [6, 2, 2, 8, 8, 6, 7, 6, 1, 8, 8, 5, 5, 4, 1, 5],
  [9, 5, 5, 4, 10, 1, 1, 6, 6, 10, 1, 7, 7, 6],
  [1, 6, 6, 10, 10, 1, 1, 7, 7, 6, 0, 7, 9, 5, 5, 4],
  [0, 10, 10, 4, 10, 5, 5, 4, 3, 10, 6, 10, 10, 7, 7, 6],
  [7, 6, 6, 10, 10, 7, 10, 8, 5, 4, 4, 10, 10, 5],
  [6, 9, 9, 5, 5, 6, 6, 11, 11, 9],
  [3, 6, 6, 11, 0, 6, 0, 5, 5, 6, 9, 5],
  [0, 11, 0, 5, 5, 11, 1, 5, 5, 6, 6, 11],
  [6, 11, 3, 6, 3, 5, 5, 6, 1, 5],
  [2, 10, 10, 1, 9, 5, 5, 11, 11, 9, 5, 6, 6, 11],
  [0, 11, 0, 6, 6, 11, 9, 6, 5, 6, 9, 5, 2, 10, 10, 1],
  [8, 5, 5, 11, 5, 6, 6, 11, 0, 5, 10, 5, 5, 2, 2, 10],
  [6, 11, 3, 6, 3, 5, 5, 6, 2, 10, 10, 3, 10, 5],
  [5, 8, 9, 5, 5, 2, 2, 8, 5, 6, 6, 2],
  [9, 5, 5, 6, 6, 9, 6, 0, 6, 2],
  [1, 5, 5, 8, 8, 1, 5, 6, 6, 8, 8, 2, 6, 2],
  [1, 5, 5, 6, 6, 1, 6, 2],
  [3, 6, 6, 1, 6, 10, 10, 1, 8, 6, 5, 6, 6, 9, 9, 5],
  [10, 1, 0, 10, 0, 6, 6, 10, 9, 5, 5, 0, 5, 6],
  [5, 6, 6, 10, 10, 5],
  [10, 5, 5, 6, 6, 10],
  [11, 5, 5, 10, 10, 11, 7, 5],
  [11, 5, 5, 10, 10, 11, 7, 5],
  [5, 11, 7, 5, 5, 10, 10, 11, 1, 9],
  [10, 7, 7, 5, 5, 10, 10, 11, 8, 1, 1, 9],
  [11, 1, 2, 11, 7, 1, 7, 5, 5, 1],
  [2, 7, 7, 1, 7, 5, 5, 1, 2, 11],
  [9, 7, 7, 5, 5, 9, 9, 2, 2, 7, 2, 11],
  [7, 5, 5, 2, 2, 7, 2, 11, 5, 9, 9, 2, 2, 8],
  [2, 5, 5, 10, 10, 2, 3, 5, 7, 5],
  [8, 2, 8, 5, 5, 2, 7, 5, 10, 2, 5, 10],
  [1, 9, 5, 10, 10, 3, 3, 5, 7, 5, 10, 2],
  [8, 2, 2, 9, 1, 9, 7, 2, 10, 2, 2, 5, 5, 10, 7, 5],
  [3, 5, 5, 1, 7, 5],
  [7, 0, 7, 1, 7, 5, 5, 1],
  [3, 9, 3, 5, 5, 9, 7, 5],
  [7, 9, 5, 9, 7, 5],
  [5, 8, 4, 5, 5, 10, 10, 8, 10, 11],
  [5, 0, 4, 5, 5, 11, 11, 0, 5, 10, 10, 11],
  [1, 9, 4, 10, 10, 8, 10, 11, 4, 5, 5, 10],
  [10, 11, 11, 4, 4, 10, 4, 5, 5, 10, 3, 4, 4, 1, 1, 9],
  [2, 5, 5, 1, 2, 8, 8, 5, 2, 11, 4, 5],
  [4, 11, 11, 0, 4, 5, 5, 11, 2, 11, 11, 1, 5, 1],
  [2, 5, 5, 0, 5, 9, 2, 11, 11, 5, 4, 5, 5, 8],
  [4, 5, 5, 9, 2, 11],
  [2, 5, 5, 10, 10, 2, 3, 5, 3, 4, 4, 5],
  [5, 10, 10, 2, 2, 5, 2, 4, 4, 5],
  [3, 10, 10, 2, 3, 5, 5, 10, 8, 5, 4, 5, 1, 9],
  [5, 10, 10, 2, 2, 5, 2, 4, 4, 5, 1, 9, 9, 2],
  [4, 5, 5, 8, 5, 3, 5, 1],
  [4, 5, 5, 0, 5, 1],
  [4, 5, 5, 8, 5, 3, 0, 5, 5, 9],
  [4, 5, 5, 9],
  [4, 11, 7, 4, 9, 11, 9, 10, 10, 11],
  [9, 7, 7, 4, 9, 11, 9, 10, 10, 11],
  [1, 10, 10, 11, 11, 1, 11, 4, 4, 1, 7, 4],
  [1, 4, 4, 3, 1, 10, 10, 4, 7, 4, 4, 11, 10, 11],
  [4, 11, 7, 4, 9, 11, 9, 2, 2, 11, 9, 1],
  [9, 7, 7, 4, 9, 11, 9, 1, 1, 11, 2, 11],
  [7, 4, 4, 11, 4, 2, 2, 11],
  [7, 4, 4, 11, 4, 2, 2, 11, 3, 4],
  [2, 9, 9, 10, 10, 2, 2, 7, 7, 9, 7, 4],
  [9, 10, 10, 7, 7, 9, 7, 4, 10, 2, 2, 7, 7, 0],
  [7, 10, 10, 3, 10, 2, 7, 4, 4, 10, 1, 10, 10, 0],
  [1, 10, 10, 2, 7, 4],
  [9, 1, 1, 4, 1, 7, 7, 4],
  [9, 1, 1, 4, 1, 7, 7, 4, 8, 1],
  [3, 4, 7, 4],
  [7, 4],
  [9, 10, 10, 8, 10, 11],
  [9, 3, 9, 11, 9, 10, 10, 11],
  [1, 10, 10, 0, 10, 8, 10, 11],
  [1, 10, 10, 3, 10, 11],
  [2, 11, 11, 1, 11, 9, 9, 1],
  [9, 3, 9, 11, 2, 9, 9, 1, 2, 11],
  [2, 11, 11, 0],
  [2, 11],
  [8, 2, 8, 10, 10, 2, 9, 10],
  [9, 10, 10, 2, 2, 9],
  [8, 2, 8, 10, 10, 2, 1, 8, 1, 10],
  [1, 10, 10, 2],
  [8, 1, 9, 1],
  [9, 1],
  [],
  []];

var segTable2 = [
  [],
  [],
  [1, 9],
  [1, 9],
  [2, 10, 10, 1],
  [2, 10, 10, 1],
  [2, 10, 10, 9],
  [2, 10, 10, 9],
  [11, 2],
  [11, 2],
  [1, 9, 11, 2],
  [11, 2, 1, 9],
  [10, 1, 11, 10],
  [10, 1, 11, 10],
  [11, 10, 10, 9],
  [10, 9, 11, 10],
  [4, 7],
  [4, 7],
  [1, 9, 4, 7],
  [1, 9, 4, 7],
  [2, 10, 10, 1, 4, 7],
  [4, 7, 2, 10, 10, 1],
  [2, 10, 10, 9, 4, 7],
  [2, 10, 10, 9, 4, 7],
  [4, 7, 11, 2],
  [4, 7, 11, 2],
  [1, 9, 4, 7, 11, 2],
  [4, 7, 11, 2, 1, 9],
  [10, 1, 11, 10, 4, 7],
  [11, 10, 10, 1, 4, 7],
  [4, 7, 11, 10, 10, 9],
  [4, 7, 11, 10, 10, 9],
  [9, 5, 5, 4],
  [9, 5, 5, 4],
  [5, 4, 1, 5],
  [5, 4, 1, 5],
  [2, 10, 10, 1, 9, 5, 5, 4],
  [2, 10, 10, 1, 9, 5, 5, 4],
  [2, 10, 10, 5, 5, 4],
  [2, 10, 10, 5, 5, 4],
  [9, 5, 5, 4, 11, 2],
  [11, 2, 9, 5, 5, 4],
  [5, 4, 1, 5, 11, 2],
  [1, 5, 11, 2, 5, 4],
  [11, 10, 10, 1, 9, 5, 5, 4],
  [9, 5, 5, 4, 10, 1, 11, 10],
  [5, 4, 11, 10, 10, 5],
  [5, 4, 10, 5, 11, 10],
  [5, 7, 9, 5],
  [9, 5, 5, 7],
  [1, 5, 5, 7],
  [1, 5, 5, 7],
  [9, 5, 5, 7, 10, 1, 2, 10],
  [10, 1, 2, 10, 9, 5, 5, 7],
  [5, 7, 10, 5, 2, 10],
  [2, 10, 10, 5, 5, 7],
  [9, 5, 5, 7, 11, 2],
  [9, 5, 5, 7, 11, 2],
  [11, 2, 1, 5, 5, 7],
  [11, 2, 1, 5, 5, 7],
  [9, 5, 5, 7, 10, 1, 11, 10],
  [5, 7, 9, 5, 10, 1, 11, 10],
  [11, 10, 10, 5, 5, 7],
  [11, 10, 10, 5, 5, 7],
  [10, 6, 6, 5, 5, 10],
  [5, 10, 10, 6, 6, 5],
  [1, 9, 5, 10, 10, 6, 6, 5],
  [1, 9, 5, 10, 10, 6, 6, 5],
  [6, 5, 5, 1, 2, 6],
  [6, 5, 5, 1, 2, 6],
  [6, 5, 5, 9, 2, 6],
  [5, 9, 2, 6, 6, 5],
  [11, 2, 10, 6, 6, 5, 5, 10],
  [11, 2, 10, 6, 6, 5, 5, 10],
  [1, 9, 11, 2, 5, 10, 10, 6, 6, 5],
  [5, 10, 10, 6, 6, 5, 1, 9, 11, 2],
  [11, 6, 6, 5, 5, 1],
  [5, 1, 11, 6, 6, 5],
  [11, 6, 6, 5, 5, 9],
  [6, 5, 5, 9, 11, 6],
  [5, 10, 10, 6, 6, 5, 4, 7],
  [4, 7, 6, 5, 5, 10, 10, 6],
  [1, 9, 5, 10, 10, 6, 6, 5, 4, 7],
  [10, 6, 6, 5, 5, 10, 1, 9, 4, 7],
  [2, 6, 6, 5, 5, 1, 4, 7],
  [5, 1, 2, 6, 6, 5, 4, 7],
  [4, 7, 5, 9, 6, 5, 2, 6],
  [4, 7, 5, 9, 6, 5, 2, 6],
  [11, 2, 4, 7, 10, 6, 6, 5, 5, 10],
  [5, 10, 10, 6, 6, 5, 4, 7, 11, 2],
  [1, 9, 4, 7, 11, 2, 5, 10, 10, 6, 6, 5],
  [1, 9, 11, 2, 4, 7, 5, 10, 10, 6, 6, 5],
  [4, 7, 5, 1, 11, 6, 6, 5],
  [5, 1, 11, 6, 6, 5, 4, 7],
  [5, 9, 6, 5, 11, 6, 4, 7],
  [6, 5, 5, 9, 11, 6, 4, 7],
  [9, 10, 6, 4, 10, 6],
  [10, 6, 6, 4, 9, 10],
  [1, 10, 10, 6, 6, 4],
  [6, 4, 1, 10, 10, 6],
  [9, 1, 2, 6, 6, 4],
  [9, 1, 2, 6, 6, 4],
  [2, 6, 6, 4],
  [2, 6, 6, 4],
  [9, 10, 10, 6, 6, 4, 11, 2],
  [11, 2, 9, 10, 10, 6, 6, 4],
  [11, 2, 6, 4, 1, 10, 10, 6],
  [6, 4, 1, 10, 10, 6, 11, 2],
  [6, 4, 9, 1, 11, 6],
  [11, 6, 9, 1, 6, 4],
  [11, 6, 6, 4],
  [6, 4, 11, 6],
  [10, 6, 6, 7, 9, 10],
  [9, 10, 6, 7, 10, 6],
  [10, 6, 6, 7, 1, 10],
  [10, 6, 6, 7, 1, 10],
  [2, 6, 9, 1, 6, 7],
  [2, 6, 9, 1, 6, 7],
  [6, 7, 2, 6],
  [6, 7, 2, 6],
  [11, 2, 10, 6, 9, 10, 6, 7],
  [11, 2, 6, 7, 10, 6, 9, 10],
  [1, 10, 6, 7, 10, 6, 11, 2],
  [11, 2, 10, 6, 1, 10, 6, 7],
  [6, 7, 9, 1, 11, 6],
  [9, 1, 11, 6, 6, 7],
  [6, 7, 11, 6],
  [11, 6, 6, 7],
  [7, 6, 6, 11],
  [7, 6, 6, 11],
  [1, 9, 7, 6, 6, 11],
  [1, 9, 7, 6, 6, 11],
  [10, 1, 2, 10, 6, 11, 7, 6],
  [2, 10, 10, 1, 6, 11, 7, 6],
  [2, 10, 10, 9, 6, 11, 7, 6],
  [6, 11, 7, 6, 2, 10, 10, 9],
  [6, 2, 7, 6],
  [7, 6, 6, 2],
  [7, 6, 6, 2, 1, 9],
  [6, 2, 1, 9, 7, 6],
  [7, 6, 6, 10, 10, 1],
  [7, 6, 6, 10, 10, 1],
  [10, 9, 6, 10, 7, 6],
  [7, 6, 6, 10, 10, 9],
  [4, 6, 6, 11],
  [6, 11, 4, 6],
  [6, 11, 4, 6, 1, 9],
  [4, 6, 1, 9, 6, 11],
  [4, 6, 6, 11, 2, 10, 10, 1],
  [2, 10, 10, 1, 6, 11, 4, 6],
  [4, 6, 6, 11, 2, 10, 10, 9],
  [10, 9, 2, 10, 6, 11, 4, 6],
  [4, 6, 6, 2],
  [4, 6, 6, 2],
  [1, 9, 4, 6, 6, 2],
  [1, 9, 4, 6, 6, 2],
  [4, 6, 6, 10, 10, 1],
  [10, 1, 6, 10, 4, 6],
  [4, 6, 6, 10, 10, 9],
  [10, 9, 6, 10, 4, 6],
  [9, 5, 5, 4, 7, 6, 6, 11],
  [9, 5, 5, 4, 7, 6, 6, 11],
  [1, 5, 5, 4, 7, 6, 6, 11],
  [7, 6, 6, 11, 5, 4, 1, 5],
  [9, 5, 5, 4, 10, 1, 2, 10, 7, 6, 6, 11],
  [6, 11, 7, 6, 2, 10, 10, 1, 9, 5, 5, 4],
  [7, 6, 6, 11, 5, 4, 10, 5, 2, 10],
  [5, 4, 10, 5, 2, 10, 7, 6, 6, 11],
  [7, 6, 6, 2, 5, 4, 9, 5],
  [9, 5, 5, 4, 6, 2, 7, 6],
  [6, 2, 7, 6, 1, 5, 5, 4],
  [6, 2, 7, 6, 5, 4, 1, 5],
  [9, 5, 5, 4, 10, 1, 6, 10, 7, 6],
  [6, 10, 10, 1, 7, 6, 9, 5, 5, 4],
  [10, 5, 5, 4, 6, 10, 7, 6],
  [7, 6, 6, 10, 5, 4, 10, 5],
  [9, 5, 5, 6, 6, 11],
  [6, 11, 5, 6, 9, 5],
  [1, 5, 5, 6, 6, 11],
  [6, 11, 5, 6, 1, 5],
  [2, 10, 10, 1, 9, 5, 5, 6, 6, 11],
  [6, 11, 5, 6, 9, 5, 2, 10, 10, 1],
  [5, 6, 6, 11, 10, 5, 2, 10],
  [6, 11, 5, 6, 2, 10, 10, 5],
  [9, 5, 5, 6, 6, 2],
  [9, 5, 5, 6, 6, 2],
  [1, 5, 5, 6, 6, 2],
  [1, 5, 5, 6, 6, 2],
  [6, 10, 10, 1, 5, 6, 9, 5],
  [10, 1, 6, 10, 9, 5, 5, 6],
  [5, 6, 6, 10, 10, 5],
  [10, 5, 5, 6, 6, 10],
  [5, 10, 10, 11, 7, 5],
  [5, 10, 10, 11, 7, 5],
  [7, 5, 5, 10, 10, 11, 1, 9],
  [7, 5, 5, 10, 10, 11, 1, 9],
  [2, 11, 7, 5, 5, 1],
  [7, 5, 5, 1, 2, 11],
  [7, 5, 5, 9, 2, 11],
  [7, 5, 2, 11, 5, 9],
  [5, 10, 10, 2, 7, 5],
  [7, 5, 10, 2, 5, 10],
  [1, 9, 5, 10, 7, 5, 10, 2],
  [1, 9, 10, 2, 5, 10, 7, 5],
  [5, 1, 7, 5],
  [7, 5, 5, 1],
  [5, 9, 7, 5],
  [5, 9, 7, 5],
  [4, 5, 5, 10, 10, 11],
  [4, 5, 5, 10, 10, 11],
  [1, 9, 10, 11, 4, 5, 5, 10],
  [10, 11, 4, 5, 5, 10, 1, 9],
  [5, 1, 2, 11, 4, 5],
  [4, 5, 2, 11, 5, 1],
  [5, 9, 2, 11, 4, 5],
  [4, 5, 5, 9, 2, 11],
  [5, 10, 10, 2, 4, 5],
  [5, 10, 10, 2, 4, 5],
  [10, 2, 5, 10, 4, 5, 1, 9],
  [5, 10, 10, 2, 4, 5, 1, 9],
  [4, 5, 5, 1],
  [4, 5, 5, 1],
  [4, 5, 5, 9],
  [4, 5, 5, 9],
  [7, 4, 9, 10, 10, 11],
  [7, 4, 9, 10, 10, 11],
  [1, 10, 10, 11, 7, 4],
  [1, 10, 7, 4, 10, 11],
  [7, 4, 2, 11, 9, 1],
  [7, 4, 9, 1, 2, 11],
  [7, 4, 2, 11],
  [7, 4, 2, 11],
  [9, 10, 10, 2, 7, 4],
  [9, 10, 7, 4, 10, 2],
  [10, 2, 7, 4, 1, 10],
  [1, 10, 10, 2, 7, 4],
  [9, 1, 7, 4],
  [9, 1, 7, 4],
  [7, 4],
  [7, 4],
  [9, 10, 10, 11],
  [9, 10, 10, 11],
  [1, 10, 10, 11],
  [1, 10, 10, 11],
  [2, 11, 9, 1],
  [9, 1, 2, 11],
  [2, 11],
  [2, 11],
  [10, 2, 9, 10],
  [9, 10, 10, 2],
  [10, 2, 1, 10],
  [1, 10, 10, 2],
  [9, 1],
  [9, 1],
  [],
  []];

var cubeVerts = [[0,0,0], [1,0,0], [1,1,0], [0,1,0],
                   [0,0,1], [1,0,1], [1,1,1], [0,1,1]];
var edgeIndex = [[0,1], [1,2], [2,3], [3,0], [4,5], [5,6],
                   [6,7], [7,4], [0,4], [1,5], [2,6], [3,7]];
// edge directions: [x, y, -x, -y, x, y, -x, -y, z, z, z, z]

// return offsets relative to vertex [0,0,0]
function calculateVertOffsets(dims) {
  var vert_offsets = [];
  for (var i = 0; i < 8; ++i) {
    var v = cubeVerts[i];
    vert_offsets.push(v[0] + dims[2] * (v[1] + dims[1] * v[2]));
  }
  return vert_offsets;
}


function marchingCubes(dims, values, points, isolevel, method) {
  var snap = (method === 'snapped MC');
  var seg_table = (method === 'squarish' ? segTable2 : segTable);
  var vlist = new Array(12);
  var vert_offsets = calculateVertOffsets(dims);
  var vertex_values = new Float32Array(8);
  var p0 /*:Num3*/ = [0, 0, 0]; // unused initial value - to make Flow happy
  var vertex_points = [p0, p0, p0, p0, p0, p0, p0, p0];
  var size_x = dims[0];
  var size_y = dims[1];
  var size_z = dims[2];
  if (values == null || points == null) { return; }
  var vertices = [];
  var segments = [];
  var vertex_count = 0;
  for (var x = 0; x < size_x - 1; x++) {
    for (var y = 0; y < size_y - 1; y++) {
      for (var z = 0; z < size_z - 1; z++) {
        var offset0 = z + size_z * (y + size_y * x);
        var cubeindex = 0;
        var i = (void 0);
        var j = (void 0);
        for (i = 0; i < 8; ++i) {
          j = offset0 + vert_offsets[i];
          cubeindex |= (values[j] < isolevel) ? 1 << i : 0;
        }
        if (cubeindex === 0 || cubeindex === 255) { continue; }
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
              if (mu > 0.85) { mu = 1; }
              else if (mu < 0.15) { mu = 0; }
            }
            var p1 = vertex_points[e[0]];
            var p2 = vertex_points[e[1]];
            // The number of added vertices could be roughly halved
            // if we avoided duplicates between neighbouring cells.
            // Using a map for lookups is too slow, perhaps a big
            // array would do?
            vertices.push(p1[0] + (p2[0] - p1[0]) * mu,
                          p1[1] + (p2[1] - p1[1]) * mu,
                          p1[2] + (p2[2] - p1[2]) * mu);
            vlist[i] = vertex_count++;
          }
        }
        var t = seg_table[cubeindex];
        for (i = 0; i < t.length; i++) {
          segments.push(vlist[t[i]]);
        }
      }
    }
  }
  return { vertices: vertices, segments: segments };
}

// @flow

function modulo(a, b) {
  var reminder = a % b;
  return reminder >= 0 ? reminder : reminder + b;
}

var GridArray = function GridArray(dim /*:number[]*/) {
  this.dim = dim; // dimensions of the grid for the entire unit cell
  this.values = new Float32Array(dim[0] * dim[1] * dim[2]);
};

GridArray.prototype.grid2index = function grid2index (i/*:number*/, j/*:number*/, k/*:number*/) {
  i = modulo(i, this.dim[0]);
  j = modulo(j, this.dim[1]);
  k = modulo(k, this.dim[2]);
  return this.dim[2] * (this.dim[1] * i + j) + k;
};

GridArray.prototype.grid2index_unchecked = function grid2index_unchecked (i/*:number*/, j/*:number*/, k/*:number*/) {
  return this.dim[2] * (this.dim[1] * i + j) + k;
};

GridArray.prototype.grid2frac = function grid2frac (i/*:number*/, j/*:number*/, k/*:number*/) {
  return [i / this.dim[0], j / this.dim[1], k / this.dim[2]];
};

// return grid coordinates (rounded down) for the given fractional coordinates
GridArray.prototype.frac2grid = function frac2grid (xyz/*:number[]*/) {
  // at one point "| 0" here made extract_block() 40% faster on V8 3.14,
  // but I don't see any effect now
  return [Math.floor(xyz[0] * this.dim[0]) | 0,
          Math.floor(xyz[1] * this.dim[1]) | 0,
          Math.floor(xyz[2] * this.dim[2]) | 0];
};

GridArray.prototype.set_grid_value = function set_grid_value (i/*:number*/, j/*:number*/, k/*:number*/, value/*:number*/) {
  var idx = this.grid2index(i, j, k);
  this.values[idx] = value;
};

GridArray.prototype.get_grid_value = function get_grid_value (i/*:number*/, j/*:number*/, k/*:number*/) {
  var idx = this.grid2index(i, j, k);
  return this.values[idx];
};

function calculate_stddev(a, offset) {
  var sum = 0;
  var sq_sum = 0;
  var alen = a.length;
  for (var i = offset; i < alen; i++) {
    sum += a[i];
    sq_sum += a[i] * a[i];
  }
  var mean = sum / (alen - offset);
  var variance = sq_sum / (alen - offset) - mean * mean;
  return {mean: mean, rms: Math.sqrt(variance)};
}

var ElMap = function ElMap() {
  this.unit_cell = null;
  this.grid = null;
  this.stats = { mean: 0.0, rms: 1.0 };
  this.block = new Block();
};

ElMap.prototype.abs_level = function abs_level (sigma /*:number*/) {
  return sigma * this.stats.rms + this.stats.mean;
};

// http://www.ccp4.ac.uk/html/maplib.html#description
// eslint-disable-next-line complexity
ElMap.prototype.from_ccp4 = function from_ccp4 (buf /*:ArrayBuffer*/, expand_symmetry /*:?boolean*/) {
  if (expand_symmetry === undefined) { expand_symmetry = true; }
  if (buf.byteLength < 1024) { throw Error('File shorter than 1024 bytes.'); }
  //console.log('buf type: ' + Object.prototype.toString.call(buf));
  // for now we assume both file and host are little endian
  var iview = new Int32Array(buf, 0, 256);
  // word 53 - character string 'MAP ' to identify file type
  if (iview[52] !== 0x2050414d) { throw Error('not a CCP4 map'); }
  // map has 3 dimensions referred to as columns (fastest changing), rows
  // and sections (c-r-s)
  var n_crs = [iview[0], iview[1], iview[2]];
  var mode = iview[3];
  var nb;
  if (mode === 2) { nb = 4; }
  else if (mode === 0) { nb = 1; }
  else { throw Error('Only Mode 2 and Mode 0 of CCP4 map is supported.'); }
  var start = [iview[4], iview[5], iview[6]];
  var n_grid = [iview[7], iview[8], iview[9]];
  var nsymbt = iview[23]; // size of extended header in bytes
  if (1024 + nsymbt + nb*n_crs[0]*n_crs[1]*n_crs[2] !== buf.byteLength) {
    throw Error('ccp4 file too short or too long');
  }
  var fview = new Float32Array(buf, 0, buf.byteLength / 4);
  this.unit_cell = new UnitCell(fview[10], fview[11], fview[12],
                                fview[13], fview[14], fview[15]);
  // MAPC, MAPR, MAPS - axis corresp to cols, rows, sections (1,2,3 for X,Y,Z)
  var map_crs = [iview[16], iview[17], iview[18]];
  var ax = map_crs.indexOf(1);
  var ay = map_crs.indexOf(2);
  var az = map_crs.indexOf(3);

  var min = fview[19];
  var max = fview[20];
  //const sg_number = iview[22];
  //const lskflg = iview[24];
  var grid = new GridArray(n_grid);
  if (nsymbt % 4 !== 0) {
    throw Error('CCP4 map with NSYMBT not divisible by 4 is not supported.');
  }
  var data_view;
  if (mode === 2) { data_view = fview; }
  else /* mode === 0 */ { data_view = new Int8Array(buf); }
  var idx = (1024 + nsymbt) / nb | 0;

  // We assume that if DMEAN and RMS from the header are not clearly wrong
  // they are what the user wants. Because the map can cover a small part
  // of the asu and its rmsd may be different than the total rmsd.
  this.stats.mean = fview[21];
  this.stats.rms = fview[54];
  if (this.stats.mean < min || this.stats.mean > max || this.stats.rms <= 0) {
    this.stats = calculate_stddev(data_view, idx);
  }
  var b1 = 1;
  var b0 = 0;
  // if the file was converted by mapmode2to0 - scale the data
  if (mode === 0 && iview[39] === -128 && iview[40] === 127) {
    // scaling f(x)=b1*x+b0 such that f(-128)=min and f(127)=max
    b1 = (max - min) / 255.0;
    b0 = 0.5 * (min + max + b1);
  }

  var end = [start[0] + n_crs[0], start[1] + n_crs[1], start[2] + n_crs[2]];
  var it = [0, 0, 0];
  for (it[2] = start[2]; it[2] < end[2]; it[2]++) { // sections
    for (it[1] = start[1]; it[1] < end[1]; it[1]++) { // rows
      for (it[0] = start[0]; it[0] < end[0]; it[0]++) { // cols
        grid.set_grid_value(it[ax], it[ay], it[az], b1 * data_view[idx] + b0);
        idx++;
      }
    }
  }
  if (expand_symmetry && nsymbt > 0) {
    var u8view = new Uint8Array(buf);
    for (var i = 0; i+80 <= nsymbt; i += 80) {
      var j = (void 0);
      var symop = '';
      for (j = 0; j < 80; ++j) {
        symop += String.fromCharCode(u8view[1024 + i + j]);
      }
      if (/^\s*x\s*,\s*y\s*,\s*z\s*$/i.test(symop)) { continue; }// skip x,y,z
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
            grid.set_grid_value(xyz[0], xyz[1], xyz[2],
                                b1 * data_view[idx] + b0);
            idx++;
          }
        }
      }
    }
  }
  this.grid = grid;
};

// DSN6 MAP FORMAT
// http://www.uoxray.uoregon.edu/tnt/manual/node104.html
// Density values are stored as bytes.
ElMap.prototype.from_dsn6 = function from_dsn6 (buf /*: ArrayBuffer*/) {
  //console.log('buf type: ' + Object.prototype.toString.call(buf));
  var u8data = new Uint8Array(buf);
  var iview = new Int16Array(u8data.buffer);
  if (iview[18] !== 100) {
    var len = iview.length;// or only header, 256?
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
  var grid = new GridArray(n_grid);
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
                grid.set_grid_value(origin[0] + x,
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
  this.stats = calculate_stddev(grid.values, 0);
  this.grid = grid;
  //this.show_debug_info();
};

ElMap.prototype.show_debug_info = function show_debug_info () {
  console.log('unit cell:', this.unit_cell && this.unit_cell.parameters);
  console.log('grid:', this.grid && this.grid.dim);
};

// Extract a block of density for calculating an isosurface using the
// separate marching cubes implementation.
ElMap.prototype.extract_block = function extract_block (radius/*:number*/, center /*:[number,number,number]*/) {
  var grid = this.grid;
  var unit_cell = this.unit_cell;
  if (grid == null || unit_cell == null) { return; }
  var fc = unit_cell.fractionalize(center);
  var r = [radius / unit_cell.parameters[0],
             radius / unit_cell.parameters[1],
             radius / unit_cell.parameters[2]];
  var grid_min = grid.frac2grid([fc[0] - r[0], fc[1] - r[1], fc[2] - r[2]]);
  var grid_max = grid.frac2grid([fc[0] + r[0], fc[1] + r[1], fc[2] + r[2]]);
  var size = [grid_max[0] - grid_min[0] + 1,
                grid_max[1] - grid_min[1] + 1,
                grid_max[2] - grid_min[2] + 1];
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
  this.block.set(points, values, size);
};

ElMap.prototype.isomesh_in_block = function isomesh_in_block (sigma/*:number*/, method/*:string*/) {
  var abs_level = this.abs_level(sigma);
  return this.block.isosurface(abs_level, method);
};

ElMap.prototype.unit = 'e/\u212B\u00B3';

// symop -> matrix ([x,y,z] = matrix * [x,y,z,1])
function parse_symop(symop) {
  var ops = symop.toLowerCase().replace(/\s+/g, '').split(',');
  if (ops.length !== 3) { throw Error('Unexpected symop: ' + symop); }
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
        if (!m) { throw Error('What is ' + terms[j] + ' in ' + symop); }
        row[3] = sign * Number(m[1]) / Number(m[2]);
      }
    }
    mat.push(row);
  }
  return mat;
}

// a small subset of THREE.js v83 that is used by UglyMol
// modified with eslint --fix and manually,
// LICENSE: threejs.org/license (MIT)

/* eslint-disable max-len, one-var, guard-for-in */
/* eslint-disable prefer-rest-params, no-invalid-this, no-useless-escape */
/* eslint-disable new-cap, no-extend-native */

// Polyfills

if ( Function.prototype.name === undefined ) {
  // Missing in IE9-11.
  // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Function/name
  Object.defineProperty( Function.prototype, 'name', {
    get: function () {
      return this.toString().match( /^\s*function\s*([^\(\s]*)/ )[1];
    },
  } );
}

if ( Object.assign === undefined ) {
  // Missing in IE.
  // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Object/assign
  ( function () {
    Object.assign = function ( target ) {
      if ( target === undefined || target === null ) {
        throw new TypeError( 'Cannot convert undefined or null to object' );
      }
      var output = Object( target );
      for ( var index = 1; index < arguments.length; index ++ ) {
        var source = arguments[index];
        if ( source !== undefined && source !== null ) {
          for ( var nextKey in source ) {
            if ( Object.prototype.hasOwnProperty.call( source, nextKey ) ) {
              output[nextKey] = source[nextKey];
            }
          }
        }
      }
      return output;
    };
  } )();
}

var NoBlending = 0;
var NormalBlending = 1;
var TrianglesDrawMode = 0;
var TriangleStripDrawMode = 1;
var TriangleFanDrawMode = 2;

/**
* @author alteredq / http://alteredqualia.com/
* @author mrdoob / http://mrdoob.com/
*/

var _Math = {

  generateUUID: function () {
    // http://www.broofa.com/Tools/Math.uuid.htm

    var chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'.split( '' );
    var uuid = new Array( 36 );
    var rnd = 0, r;

    return function generateUUID() {
      for ( var i = 0; i < 36; i ++ ) {
        if ( i === 8 || i === 13 || i === 18 || i === 23 ) {
          uuid[i] = '-';
        } else if ( i === 14 ) {
          uuid[i] = '4';
        } else {
          if ( rnd <= 0x02 ) { rnd = 0x2000000 + ( Math.random() * 0x1000000 ) | 0; }
          r = rnd & 0xf;
          rnd = rnd >> 4;
          uuid[i] = chars[( i === 19 ) ? ( r & 0x3 ) | 0x8 : r];
        }
      }

      return uuid.join( '' );
    };
  }(),

  clamp: function ( value, min, max ) {
    return Math.max( min, Math.min( max, value ) );
  },

  // compute euclidian modulo of m % n
  // https://en.wikipedia.org/wiki/Modulo_operation

  euclideanModulo: function ( n, m ) {
    return ( ( n % m ) + m ) % m;
  },
};

/**
* @author mikael emtinger / http://gomo.se/
* @author alteredq / http://alteredqualia.com/
* @author WestLangley / http://github.com/WestLangley
* @author bhouston / http://clara.io
*/

function Quaternion( x, y, z, w ) {
  this._x = x || 0;
  this._y = y || 0;
  this._z = z || 0;
  this._w = ( w !== undefined ) ? w : 1;
}

Quaternion.prototype = {

  constructor: Quaternion,

  get x() {
    return this._x;
  },

  get y() {
    return this._y;
  },

  get z() {
    return this._z;
  },

  get w() {
    return this._w;
  },

  setFromAxisAngle: function ( axis, angle ) {
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm

    // assumes axis is normalized
    var halfAngle = angle / 2, s = Math.sin( halfAngle );

    this._x = axis.x * s;
    this._y = axis.y * s;
    this._z = axis.z * s;
    this._w = Math.cos( halfAngle );

    return this;
  },

  setFromRotationMatrix: function ( m ) {
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm

    // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

    var te = m.elements,

      m11 = te[0], m12 = te[4], m13 = te[8],
      m21 = te[1], m22 = te[5], m23 = te[9],
      m31 = te[2], m32 = te[6], m33 = te[10],

      trace = m11 + m22 + m33,
      s;

    if ( trace > 0 ) {
      s = 0.5 / Math.sqrt( trace + 1.0 );

      this._w = 0.25 / s;
      this._x = ( m32 - m23 ) * s;
      this._y = ( m13 - m31 ) * s;
      this._z = ( m21 - m12 ) * s;
    } else if ( m11 > m22 && m11 > m33 ) {
      s = 2.0 * Math.sqrt( 1.0 + m11 - m22 - m33 );

      this._w = ( m32 - m23 ) / s;
      this._x = 0.25 * s;
      this._y = ( m12 + m21 ) / s;
      this._z = ( m13 + m31 ) / s;
    } else if ( m22 > m33 ) {
      s = 2.0 * Math.sqrt( 1.0 + m22 - m11 - m33 );

      this._w = ( m13 - m31 ) / s;
      this._x = ( m12 + m21 ) / s;
      this._y = 0.25 * s;
      this._z = ( m23 + m32 ) / s;
    } else {
      s = 2.0 * Math.sqrt( 1.0 + m33 - m11 - m22 );

      this._w = ( m21 - m12 ) / s;
      this._x = ( m13 + m31 ) / s;
      this._y = ( m23 + m32 ) / s;
      this._z = 0.25 * s;
    }

    return this;
  },

  setFromUnitVectors: function () {
    // http://lolengine.net/blog/2014/02/24/quaternion-from-two-vectors-final

    // assumes direction vectors vFrom and vTo are normalized

    var v1, r;

    var EPS = 0.000001;

    return function setFromUnitVectors( vFrom, vTo ) {
      if ( v1 === undefined ) { v1 = new Vector3(); }

      r = vFrom.dot( vTo ) + 1;

      if ( r < EPS ) {
        r = 0;

        if ( Math.abs( vFrom.x ) > Math.abs( vFrom.z ) ) {
          v1.set( - vFrom.y, vFrom.x, 0 );
        } else {
          v1.set( 0, - vFrom.z, vFrom.y );
        }
      } else {
        v1.crossVectors( vFrom, vTo );
      }

      this._x = v1.x;
      this._y = v1.y;
      this._z = v1.z;
      this._w = r;

      return this.normalize();
    };
  }(),

  length: function () {
    return Math.sqrt( this._x * this._x + this._y * this._y + this._z * this._z + this._w * this._w );
  },

  normalize: function () {
    var l = this.length();

    if ( l === 0 ) {
      this._x = 0;
      this._y = 0;
      this._z = 0;
      this._w = 1;
    } else {
      l = 1 / l;

      this._x = this._x * l;
      this._y = this._y * l;
      this._z = this._z * l;
      this._w = this._w * l;
    }

    return this;
  },

};

/**
* @author mrdoob / http://mrdoob.com/
* @author *kile / http://kile.stravaganza.org/
* @author philogb / http://blog.thejit.org/
* @author mikael emtinger / http://gomo.se/
* @author egraether / http://egraether.com/
* @author WestLangley / http://github.com/WestLangley
*/

function Vector3( x, y, z ) {
  this.x = x || 0;
  this.y = y || 0;
  this.z = z || 0;
}

Vector3.prototype = {

  constructor: Vector3,

  isVector3: true,

  set: function ( x, y, z ) {
    this.x = x;
    this.y = y;
    this.z = z;

    return this;
  },

  clone: function () {
    return new this.constructor( this.x, this.y, this.z );
  },

  copy: function ( v ) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;

    return this;
  },

  add: function ( v ) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;

    return this;
  },

  addVectors: function ( a, b ) {
    this.x = a.x + b.x;
    this.y = a.y + b.y;
    this.z = a.z + b.z;

    return this;
  },

  addScaledVector: function ( v, s ) {
    this.x += v.x * s;
    this.y += v.y * s;
    this.z += v.z * s;

    return this;
  },

  sub: function ( v ) {
    this.x -= v.x;
    this.y -= v.y;
    this.z -= v.z;

    return this;
  },

  subVectors: function ( a, b ) {
    this.x = a.x - b.x;
    this.y = a.y - b.y;
    this.z = a.z - b.z;

    return this;
  },

  multiplyScalar: function ( scalar ) {
    if ( isFinite( scalar ) ) {
      this.x *= scalar;
      this.y *= scalar;
      this.z *= scalar;
    } else {
      this.x = 0;
      this.y = 0;
      this.z = 0;
    }

    return this;
  },

  applyMatrix4: function ( m ) {
    // input: Matrix4 affine matrix

    var x = this.x, y = this.y, z = this.z;
    var e = m.elements;

    this.x = e[0] * x + e[4] * y + e[8] * z + e[12];
    this.y = e[1] * x + e[5] * y + e[9] * z + e[13];
    this.z = e[2] * x + e[6] * y + e[10] * z + e[14];

    return this;
  },

  applyProjection: function ( m ) {
    // input: Matrix4 projection matrix

    var x = this.x, y = this.y, z = this.z;
    var e = m.elements;
    var d = 1 / ( e[3] * x + e[7] * y + e[11] * z + e[15] ); // perspective divide

    this.x = ( e[0] * x + e[4] * y + e[8] * z + e[12] ) * d;
    this.y = ( e[1] * x + e[5] * y + e[9] * z + e[13] ) * d;
    this.z = ( e[2] * x + e[6] * y + e[10] * z + e[14] ) * d;

    return this;
  },

  applyQuaternion: function ( q ) {
    var x = this.x, y = this.y, z = this.z;
    var qx = q.x, qy = q.y, qz = q.z, qw = q.w;

    // calculate quat * vector

    var ix = qw * x + qy * z - qz * y;
    var iy = qw * y + qz * x - qx * z;
    var iz = qw * z + qx * y - qy * x;
    var iw = - qx * x - qy * y - qz * z;

    // calculate result * inverse quat

    this.x = ix * qw + iw * - qx + iy * - qz - iz * - qy;
    this.y = iy * qw + iw * - qy + iz * - qx - ix * - qz;
    this.z = iz * qw + iw * - qz + ix * - qy - iy * - qx;

    return this;
  },

  unproject: function () {
    var matrix;

    return function unproject( camera ) {
      if ( matrix === undefined ) { matrix = new Matrix4(); }

      matrix.multiplyMatrices( camera.matrixWorld, matrix.getInverse( camera.projectionMatrix ) );
      return this.applyProjection( matrix );
    };
  }(),

  transformDirection: function ( m ) {
    // input: Matrix4 affine matrix
    // vector interpreted as a direction

    var x = this.x, y = this.y, z = this.z;
    var e = m.elements;

    this.x = e[0] * x + e[4] * y + e[8] * z;
    this.y = e[1] * x + e[5] * y + e[9] * z;
    this.z = e[2] * x + e[6] * y + e[10] * z;

    return this.normalize();
  },

  divideScalar: function ( scalar ) {
    return this.multiplyScalar( 1 / scalar );
  },

  dot: function ( v ) {
    return this.x * v.x + this.y * v.y + this.z * v.z;
  },

  lengthSq: function () {
    return this.x * this.x + this.y * this.y + this.z * this.z;
  },

  length: function () {
    return Math.sqrt( this.x * this.x + this.y * this.y + this.z * this.z );
  },

  normalize: function () {
    return this.divideScalar( this.length() );
  },

  setLength: function ( length ) {
    return this.multiplyScalar( length / this.length() );
  },

  lerp: function ( v, alpha ) {
    this.x += ( v.x - this.x ) * alpha;
    this.y += ( v.y - this.y ) * alpha;
    this.z += ( v.z - this.z ) * alpha;

    return this;
  },

  cross: function ( v ) {
    var x = this.x, y = this.y, z = this.z;

    this.x = y * v.z - z * v.y;
    this.y = z * v.x - x * v.z;
    this.z = x * v.y - y * v.x;

    return this;
  },

  crossVectors: function ( a, b ) {
    var ax = a.x, ay = a.y, az = a.z;
    var bx = b.x, by = b.y, bz = b.z;

    this.x = ay * bz - az * by;
    this.y = az * bx - ax * bz;
    this.z = ax * by - ay * bx;

    return this;
  },

  projectOnVector: function ( vector ) {
    var scalar = vector.dot( this ) / vector.lengthSq();
    return this.copy( vector ).multiplyScalar( scalar );
  },

  projectOnPlane: function () {
    var v1;
    return function projectOnPlane( planeNormal ) {
      if ( v1 === undefined ) { v1 = new Vector3(); }
      v1.copy( this ).projectOnVector( planeNormal );
      return this.sub( v1 );
    };
  }(),

  distanceTo: function ( v ) {
    return Math.sqrt( this.distanceToSquared( v ) );
  },

  distanceToSquared: function ( v ) {
    var dx = this.x - v.x, dy = this.y - v.y, dz = this.z - v.z;
    return dx * dx + dy * dy + dz * dz;
  },

  setFromMatrixPosition: function ( m ) {
    return this.setFromMatrixColumn( m, 3 );
  },

  setFromMatrixColumn: function ( m, index ) {
    return this.fromArray( m.elements, index * 4 );
  },

  equals: function ( v ) {
    return ( ( v.x === this.x ) && ( v.y === this.y ) && ( v.z === this.z ) );
  },

  fromArray: function ( array, offset ) {
    if ( offset === undefined ) { offset = 0; }

    this.x = array[offset];
    this.y = array[offset + 1];
    this.z = array[offset + 2];

    return this;
  },
};

/**
* @author mrdoob / http://mrdoob.com/
* @author supereggbert / http://www.paulbrunt.co.uk/
* @author philogb / http://blog.thejit.org/
* @author jordi_ros / http://plattsoft.com
* @author D1plo1d / http://github.com/D1plo1d
* @author alteredq / http://alteredqualia.com/
* @author mikael emtinger / http://gomo.se/
* @author timknip / http://www.floorplanner.com/
* @author bhouston / http://clara.io
* @author WestLangley / http://github.com/WestLangley
*/

function Matrix4() {
  this.elements = new Float32Array( [
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1 ] );
}

Matrix4.prototype = {

  constructor: Matrix4,

  isMatrix4: true,


  copy: function ( m ) {
    this.elements.set( m.elements );

    return this;
  },

  makeRotationFromQuaternion: function ( q ) {
    var te = this.elements;

    var x = q.x, y = q.y, z = q.z, w = q.w;
    var x2 = x + x, y2 = y + y, z2 = z + z;
    var xx = x * x2, xy = x * y2, xz = x * z2;
    var yy = y * y2, yz = y * z2, zz = z * z2;
    var wx = w * x2, wy = w * y2, wz = w * z2;

    te[0] = 1 - ( yy + zz );
    te[4] = xy - wz;
    te[8] = xz + wy;

    te[1] = xy + wz;
    te[5] = 1 - ( xx + zz );
    te[9] = yz - wx;

    te[2] = xz - wy;
    te[6] = yz + wx;
    te[10] = 1 - ( xx + yy );

    // last column
    te[3] = 0;
    te[7] = 0;
    te[11] = 0;

    // bottom row
    te[12] = 0;
    te[13] = 0;
    te[14] = 0;
    te[15] = 1;

    return this;
  },

  lookAt: function () {
    var x, y, z;

    return function lookAt( eye, target, up ) {
      if ( x === undefined ) {
        x = new Vector3();
        y = new Vector3();
        z = new Vector3();
      }

      var te = this.elements;

      z.subVectors( eye, target ).normalize();

      if ( z.lengthSq() === 0 ) {
        z.z = 1;
      }

      x.crossVectors( up, z ).normalize();

      if ( x.lengthSq() === 0 ) {
        z.z += 0.0001;
        x.crossVectors( up, z ).normalize();
      }

      y.crossVectors( z, x );


      te[0] = x.x; te[4] = y.x; te[8] = z.x;
      te[1] = x.y; te[5] = y.y; te[9] = z.y;
      te[2] = x.z; te[6] = y.z; te[10] = z.z;

      return this;
    };
  }(),

  multiplyMatrices: function ( a, b ) {
    var ae = a.elements;
    var be = b.elements;
    var te = this.elements;

    var a11 = ae[0], a12 = ae[4], a13 = ae[8], a14 = ae[12];
    var a21 = ae[1], a22 = ae[5], a23 = ae[9], a24 = ae[13];
    var a31 = ae[2], a32 = ae[6], a33 = ae[10], a34 = ae[14];
    var a41 = ae[3], a42 = ae[7], a43 = ae[11], a44 = ae[15];

    var b11 = be[0], b12 = be[4], b13 = be[8], b14 = be[12];
    var b21 = be[1], b22 = be[5], b23 = be[9], b24 = be[13];
    var b31 = be[2], b32 = be[6], b33 = be[10], b34 = be[14];
    var b41 = be[3], b42 = be[7], b43 = be[11], b44 = be[15];

    te[0] = a11 * b11 + a12 * b21 + a13 * b31 + a14 * b41;
    te[4] = a11 * b12 + a12 * b22 + a13 * b32 + a14 * b42;
    te[8] = a11 * b13 + a12 * b23 + a13 * b33 + a14 * b43;
    te[12] = a11 * b14 + a12 * b24 + a13 * b34 + a14 * b44;

    te[1] = a21 * b11 + a22 * b21 + a23 * b31 + a24 * b41;
    te[5] = a21 * b12 + a22 * b22 + a23 * b32 + a24 * b42;
    te[9] = a21 * b13 + a22 * b23 + a23 * b33 + a24 * b43;
    te[13] = a21 * b14 + a22 * b24 + a23 * b34 + a24 * b44;

    te[2] = a31 * b11 + a32 * b21 + a33 * b31 + a34 * b41;
    te[6] = a31 * b12 + a32 * b22 + a33 * b32 + a34 * b42;
    te[10] = a31 * b13 + a32 * b23 + a33 * b33 + a34 * b43;
    te[14] = a31 * b14 + a32 * b24 + a33 * b34 + a34 * b44;

    te[3] = a41 * b11 + a42 * b21 + a43 * b31 + a44 * b41;
    te[7] = a41 * b12 + a42 * b22 + a43 * b32 + a44 * b42;
    te[11] = a41 * b13 + a42 * b23 + a43 * b33 + a44 * b43;
    te[15] = a41 * b14 + a42 * b24 + a43 * b34 + a44 * b44;

    return this;
  },

  setPosition: function ( v ) {
    var te = this.elements;

    te[12] = v.x;
    te[13] = v.y;
    te[14] = v.z;

    return this;
  },

  getInverse: function ( m, throwOnDegenerate ) {
    // based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
    var te = this.elements,
      me = m.elements,

      n11 = me[0], n21 = me[1], n31 = me[2], n41 = me[3],
      n12 = me[4], n22 = me[5], n32 = me[6], n42 = me[7],
      n13 = me[8], n23 = me[9], n33 = me[10], n43 = me[11],
      n14 = me[12], n24 = me[13], n34 = me[14], n44 = me[15],

      t11 = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44,
      t12 = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44,
      t13 = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44,
      t14 = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;

    var det = n11 * t11 + n21 * t12 + n31 * t13 + n41 * t14;

    if ( det === 0 ) {
      var msg = 'Matrix4.getInverse(): can\'t invert matrix, determinant is 0';

      if ( throwOnDegenerate === true ) {
        throw new Error( msg );
      } else {
        console.warn( msg );
      }

      return this.identity();
    }

    var detInv = 1 / det;

    te[0] = t11 * detInv;
    te[1] = ( n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 + n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44 ) * detInv;
    te[2] = ( n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 - n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44 ) * detInv;
    te[3] = ( n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 + n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43 ) * detInv;

    te[4] = t12 * detInv;
    te[5] = ( n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 - n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44 ) * detInv;
    te[6] = ( n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 + n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44 ) * detInv;
    te[7] = ( n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 - n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43 ) * detInv;

    te[8] = t13 * detInv;
    te[9] = ( n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 + n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44 ) * detInv;
    te[10] = ( n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 - n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44 ) * detInv;
    te[11] = ( n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 + n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43 ) * detInv;

    te[12] = t14 * detInv;
    te[13] = ( n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 - n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34 ) * detInv;
    te[14] = ( n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 + n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34 ) * detInv;
    te[15] = ( n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33 ) * detInv;

    return this;
  },

  scale: function ( v ) {
    var te = this.elements;
    var x = v.x, y = v.y, z = v.z;

    te[0] *= x; te[4] *= y; te[8] *= z;
    te[1] *= x; te[5] *= y; te[9] *= z;
    te[2] *= x; te[6] *= y; te[10] *= z;
    te[3] *= x; te[7] *= y; te[11] *= z;

    return this;
  },

  getMaxScaleOnAxis: function () {
    var te = this.elements;

    var scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
    var scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
    var scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];

    return Math.sqrt( Math.max( scaleXSq, scaleYSq, scaleZSq ) );
  },

  compose: function ( position, quaternion, scale ) {
    this.makeRotationFromQuaternion( quaternion );
    this.scale( scale );
    this.setPosition( position );

    return this;
  },

  makeOrthographic: function ( left, right, top, bottom, near, far ) {
    var te = this.elements;
    var w = 1.0 / ( right - left );
    var h = 1.0 / ( top - bottom );
    var p = 1.0 / ( far - near );

    var x = ( right + left ) * w;
    var y = ( top + bottom ) * h;
    var z = ( far + near ) * p;

    te[0] = 2 * w; te[4] = 0; te[8] = 0; te[12] = - x;
    te[1] = 0; te[5] = 2 * h; te[9] = 0; te[13] = - y;
    te[2] = 0; te[6] = 0; te[10] = - 2 * p; te[14] = - z;
    te[3] = 0; te[7] = 0; te[11] = 0; te[15] = 1;

    return this;
  },
};

/**
* https://github.com/mrdoob/eventdispatcher.js/
*/

function EventDispatcher() {}

Object.assign( EventDispatcher.prototype, {

  addEventListener: function ( type, listener ) {
    if ( this._listeners === undefined ) { this._listeners = {}; }

    var listeners = this._listeners;

    if ( listeners[type] === undefined ) {
      listeners[type] = [];
    }

    if ( listeners[type].indexOf( listener ) === - 1 ) {
      listeners[type].push( listener );
    }
  },

  removeEventListener: function ( type, listener ) {
    if ( this._listeners === undefined ) { return; }

    var listeners = this._listeners;
    var listenerArray = listeners[type];

    if ( listenerArray !== undefined ) {
      var index = listenerArray.indexOf( listener );

      if ( index !== - 1 ) {
        listenerArray.splice( index, 1 );
      }
    }
  },

  dispatchEvent: function ( event ) {
    if ( this._listeners === undefined ) { return; }

    var listeners = this._listeners;
    var listenerArray = listeners[event.type];

    if ( listenerArray !== undefined ) {
      event.target = this;

      var array = [], i = 0;
      var length = listenerArray.length;

      for ( i = 0; i < length; i ++ ) {
        array[i] = listenerArray[i];
      }

      for ( i = 0; i < length; i ++ ) {
        array[i].call( this, event );
      }
    }
  },

} );

/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
* @author szimek / https://github.com/szimek/
*/

var textureId = 0;

function Texture( image ) {
  Object.defineProperty( this, 'id', { value: textureId ++ } );

  this.uuid = _Math.generateUUID();

  this.name = '';

  this.image = image;

  this.version = 0;
}

Texture.prototype = {

  constructor: Texture,

  isTexture: true,

  set needsUpdate( value ) {
    if ( value === true ) { this.version ++; }
  },

  dispose: function () {
    this.dispatchEvent( { type: 'dispose' } );
  },
};

Object.assign( Texture.prototype, EventDispatcher.prototype );


/**
* @author tschw
*
* Uniforms of a program.
* Those form a tree structure with a special top-level container for the root,
* which you get by calling 'new WebGLUniforms( gl, program, renderer )'.
*
*
* Properties of inner nodes including the top-level container:
*
* .seq - array of nested uniforms
* .map - nested uniforms by name
*
*
* Methods of all nodes except the top-level container:
*
* .setValue( gl, value, [renderer] )
*
*     uploads a uniform value(s)
*   the 'renderer' parameter is needed for sampler uniforms
*
*
* Static methods of the top-level container (renderer factorizations):
*
* .upload( gl, seq, values, renderer )
*
*     sets uniforms in 'seq' to 'values[id].value'
*
* .seqWithValue( seq, values ) : filteredSeq
*
*     filters 'seq' entries with corresponding entry in values
*
*
* Methods of the top-level container (renderer factorizations):
*
* .setValue( gl, name, value )
*
*     sets uniform with  name 'name' to 'value'
*
* .set( gl, obj, prop )
*
*     sets uniform from object and property with same name than uniform
*
* .setOptional( gl, obj, prop )
*
*     like .set for an optional property of the object
*
*/

var emptyTexture = new Texture();

// --- Base for inner nodes (including the root) ---

function UniformContainer() {
  this.seq = [];
  this.map = {};
}

// --- Setters ---

// Note: Defining these methods externally, because they come in a bunch
// and this way their names minify.

// Single scalar

function setValue1f( gl, v ) { gl.uniform1f( this.addr, v ); }
function setValue1i( gl, v ) { gl.uniform1i( this.addr, v ); }

// Single float vector (from flat array or VectorN)

function setValue2fv( gl, v ) {
  if ( v.x === undefined ) { gl.uniform2fv( this.addr, v ); }
  else { gl.uniform2f( this.addr, v.x, v.y ); }
}

function setValue3fv( gl, v ) {
  if ( v.x !== undefined ) { gl.uniform3f( this.addr, v.x, v.y, v.z ); } else if ( v.r !== undefined ) { gl.uniform3f( this.addr, v.r, v.g, v.b ); } else { gl.uniform3fv( this.addr, v ); }
}

function setValue4fv( gl, v ) {
  if ( v.x === undefined ) { gl.uniform4fv( this.addr, v ); }
  else { gl.uniform4f( this.addr, v.x, v.y, v.z, v.w ); }
}

// Single matrix (from flat array or MatrixN)

function setValue2fm( gl, v ) {
  gl.uniformMatrix2fv( this.addr, false, v.elements || v );
}

function setValue3fm( gl, v ) {
  gl.uniformMatrix3fv( this.addr, false, v.elements || v );
}

function setValue4fm( gl, v ) {
  gl.uniformMatrix4fv( this.addr, false, v.elements || v );
}

// Single texture (2D / Cube)

function setValueT1( gl, v, renderer ) {
  var unit = renderer.allocTextureUnit();
  gl.uniform1i( this.addr, unit );
  renderer.setTexture2D( v || emptyTexture, unit );
}

// Integer / Boolean vectors or arrays thereof (always flat arrays)

function setValue2iv( gl, v ) { gl.uniform2iv( this.addr, v ); }
function setValue3iv( gl, v ) { gl.uniform3iv( this.addr, v ); }
function setValue4iv( gl, v ) { gl.uniform4iv( this.addr, v ); }

// Helper to pick the right setter for the singular case

function getSingularSetter( type ) {
  switch ( type ) {
    case 0x1406: return setValue1f; // FLOAT
    case 0x8b50: return setValue2fv; // _VEC2
    case 0x8b51: return setValue3fv; // _VEC3
    case 0x8b52: return setValue4fv; // _VEC4

    case 0x8b5a: return setValue2fm; // _MAT2
    case 0x8b5b: return setValue3fm; // _MAT3
    case 0x8b5c: return setValue4fm; // _MAT4

    case 0x8b5e: return setValueT1; // SAMPLER_2D

    case 0x1404: case 0x8b56: return setValue1i; // INT, BOOL
    case 0x8b53: case 0x8b57: return setValue2iv; // _VEC2
    case 0x8b54: case 0x8b58: return setValue3iv; // _VEC3
    case 0x8b55: case 0x8b59: return setValue4iv; // _VEC4
  }
}


// --- Uniform Classes ---

function SingleUniform( id, activeInfo, addr ) {
  this.id = id;
  this.addr = addr;
  this.setValue = getSingularSetter( activeInfo.type );

  // this.path = activeInfo.name; // DEBUG
}

// --- Top-level ---

// Parser - builds up the property tree from the path strings

var RePathPart = /([\w\d_]+)(\])?(\[|\.)?/g;

// extracts
//  - the identifier (member name or array index)
//  - followed by an optional right bracket (found when array index)
//  - followed by an optional left bracket or dot (type of subscript)
//
// Note: These portions can be read in a non-overlapping fashion and
// allow straightforward parsing of the hierarchy that WebGL encodes
// in the uniform names.

function addUniform( container, uniformObject ) {
  container.seq.push( uniformObject );
  container.map[uniformObject.id] = uniformObject;
}

function parseUniform( activeInfo, addr, container ) {
  var path = activeInfo.name;

  // reset RegExp object, because of the early exit of a previous run
  RePathPart.lastIndex = 0;

  for (; ;) {
    var match = RePathPart.exec( path ),
      id = match[1],
      idIsIndex = match[2] === ']',
      subscript = match[3];

    if ( idIsIndex ) { id = id | 0; } // convert to integer

    if ( subscript === undefined ) {
      addUniform( container, new SingleUniform( id, activeInfo, addr ) );
      break;
    }
  }
}

// Root Container

function WebGLUniforms( gl, program, renderer ) {
  UniformContainer.call( this );

  this.renderer = renderer;

  var n = gl.getProgramParameter( program, gl.ACTIVE_UNIFORMS );

  for ( var i = 0; i !== n; ++ i ) {
    var info = gl.getActiveUniform( program, i ),
      path = info.name,
      addr = gl.getUniformLocation( program, path );

    parseUniform( info, addr, this );
  }
}

WebGLUniforms.prototype.setValue = function ( gl, name, value ) {
  var u = this.map[name];
  if ( u !== undefined ) { u.setValue( gl, value, this.renderer ); }
};

WebGLUniforms.prototype.set = function ( gl, object, name ) {
  var u = this.map[name];
  if ( u !== undefined ) { u.setValue( gl, object[name], this.renderer ); }
};

// Static interface

WebGLUniforms.upload = function ( gl, seq, values, renderer ) {
  for ( var i = 0, n = seq.length; i !== n; ++ i ) {
    var u = seq[i],
      v = values[u.id];

    if ( v.needsUpdate !== false ) {
      // note: always updating when .needsUpdate is undefined

      u.setValue( gl, v.value, renderer );
    }
  }
};

WebGLUniforms.seqWithValue = function ( seq, values ) {
  var r = [];

  for ( var i = 0, n = seq.length; i !== n; ++ i ) {
    var u = seq[i];
    if ( u.id in values ) { r.push( u ); }
  }

  return r;
};

/**
* @author supereggbert / http://www.paulbrunt.co.uk/
* @author philogb / http://blog.thejit.org/
* @author mikael emtinger / http://gomo.se/
* @author egraether / http://egraether.com/
* @author WestLangley / http://github.com/WestLangley
*/

function Vector4( x, y, z, w ) {
  this.x = x || 0;
  this.y = y || 0;
  this.z = z || 0;
  this.w = ( w !== undefined ) ? w : 1;
}

Vector4.prototype = {

  constructor: Vector4,

  isVector4: true,

  set: function ( x, y, z, w ) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.w = w;

    return this;
  },

  copy: function ( v ) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;
    this.w = ( v.w !== undefined ) ? v.w : 1;

    return this;
  },

  multiplyScalar: function ( scalar ) {
    if ( isFinite( scalar ) ) {
      this.x *= scalar;
      this.y *= scalar;
      this.z *= scalar;
      this.w *= scalar;
    } else {
      this.x = 0;
      this.y = 0;
      this.z = 0;
      this.w = 0;
    }

    return this;
  },

  equals: function ( v ) {
    return ( ( v.x === this.x ) && ( v.y === this.y ) && ( v.z === this.z ) && ( v.w === this.w ) );
  },
};

/**
* @author mrdoob / http://mrdoob.com/
*/

function Color( r, g, b ) {
  if ( g === undefined && b === undefined ) {
    // r is Color, hex or string
    return this.set( r );
  }

  return this.setRGB( r, g, b );
}

Color.prototype = {

  constructor: Color,

  isColor: true,

  r: 1, g: 1, b: 1,

  set: function ( value ) {
    if ( value && value.isColor ) {
      this.copy( value );
    } else if ( typeof value === 'number' ) {
      this.setHex( value );
    }

    return this;
  },

  setHex: function ( hex ) {
    hex = Math.floor( hex );

    this.r = ( hex >> 16 & 255 ) / 255;
    this.g = ( hex >> 8 & 255 ) / 255;
    this.b = ( hex & 255 ) / 255;

    return this;
  },

  setRGB: function ( r, g, b ) {
    this.r = r;
    this.g = g;
    this.b = b;

    return this;
  },

  setHSL: function () {
    function hue2rgb( p, q, t ) {
      if ( t < 0 ) { t += 1; }
      if ( t > 1 ) { t -= 1; }
      if ( t < 1 / 6 ) { return p + ( q - p ) * 6 * t; }
      if ( t < 1 / 2 ) { return q; }
      if ( t < 2 / 3 ) { return p + ( q - p ) * 6 * ( 2 / 3 - t ); }
      return p;
    }

    return function setHSL( h, s, l ) {
      // h,s,l ranges are in 0.0 - 1.0
      h = _Math.euclideanModulo( h, 1 );
      s = _Math.clamp( s, 0, 1 );
      l = _Math.clamp( l, 0, 1 );

      if ( s === 0 ) {
        this.r = this.g = this.b = l;
      } else {
        var p = l <= 0.5 ? l * ( 1 + s ) : l + s - ( l * s );
        var q = ( 2 * l ) - p;

        this.r = hue2rgb( q, p, h + 1 / 3 );
        this.g = hue2rgb( q, p, h );
        this.b = hue2rgb( q, p, h - 1 / 3 );
      }

      return this;
    };
  }(),

  getHSL: function () {
    // h,s,l ranges are in 0.0 - 1.0
    var hsl = { h: 0, s: 0, l: 0 };
    var r = this.r, g = this.g, b = this.b;
    var max = Math.max( r, g, b );
    var min = Math.min( r, g, b );
    var hue, saturation;
    var lightness = ( min + max ) / 2.0;
    if ( min === max ) {
      hue = 0;
      saturation = 0;
    } else {
      var delta = max - min;
      saturation = lightness <= 0.5 ? delta / ( max + min ) : delta / ( 2 - max - min );
      switch ( max ) {
        case r: hue = ( g - b ) / delta + ( g < b ? 6 : 0 ); break;
        case g: hue = ( b - r ) / delta + 2; break;
        case b: hue = ( r - g ) / delta + 4; break;
      }
      hue /= 6;
    }
    hsl.h = hue;
    hsl.s = saturation;
    hsl.l = lightness;
    return hsl;
  },

  clone: function () {
    return new this.constructor( this.r, this.g, this.b );
  },

  copy: function ( color ) {
    this.r = color.r;
    this.g = color.g;
    this.b = color.b;

    return this;
  },

  getHex: function () {
    return ( this.r * 255 ) << 16 ^ ( this.g * 255 ) << 8 ^ ( this.b * 255 ) << 0;
  },

  getHexString: function () {
    return ( '000000' + this.getHex().toString( 16 ) ).slice( - 6 );
  },
};

/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
*/

var materialId = 0;

function Material() {
  Object.defineProperty( this, 'id', { value: materialId ++ } );

  this.uuid = _Math.generateUUID();

  this.name = '';
  this.type = 'Material';

  this.fog = true;

  this.opacity = 1;
  this.transparent = false;

  this.depthTest = true;
  this.depthWrite = true;

  this.precision = null; // override the renderer's default precision for this material

  this.premultipliedAlpha = false;

  this.visible = true;

  this._needsUpdate = true;
}

Material.prototype = {

  constructor: Material,

  isMaterial: true,

  get needsUpdate() {
    return this._needsUpdate;
  },

  set needsUpdate( value ) {
    if ( value === true ) { this.update(); }
    this._needsUpdate = value;
  },

  setValues: function ( values ) {
    for ( var key in values ) {
      var newValue = values[key];
      this[key] = newValue;
    }
  },

  update: function () {
    this.dispatchEvent( { type: 'update' } );
  },

  dispose: function () {
    this.dispatchEvent( { type: 'dispose' } );
  },

};

Object.assign( Material.prototype, EventDispatcher.prototype );

/**
* @author alteredq / http://alteredqualia.com/
*
* parameters = {
*  uniforms: { "parameter1": { value: 1.0 }, "parameter2": { value2: 2 } },
*
*  fragmentShader: <string>,
*  vertexShader: <string>,
* }
*/

function ShaderMaterial( parameters ) {
  Material.call( this );

  this.type = 'ShaderMaterial';

  this.uniforms = {};

  this.vertexShader = '';
  this.fragmentShader = '';

  this.linewidth = 1;

  this.fog = false; // set to use scene fog

  this.extensions = {
    fragDepth: false, // set to use fragment depth values
  };

  this.setValues( parameters );
}

ShaderMaterial.prototype = Object.create( Material.prototype );
ShaderMaterial.prototype.constructor = ShaderMaterial;

ShaderMaterial.prototype.isShaderMaterial = true;


/**
* @author bhouston / http://clara.io
*/

function Ray( origin, direction ) {
  this.origin = ( origin !== undefined ) ? origin : new Vector3();
  this.direction = ( direction !== undefined ) ? direction : new Vector3();
}

Ray.prototype = {

  constructor: Ray,

  copy: function ( ray ) {
    this.origin.copy( ray.origin );
    this.direction.copy( ray.direction );

    return this;
  },

  distanceSqToPoint: function () {
    var v1 = new Vector3();

    return function distanceSqToPoint( point ) {
      var directionDistance = v1.subVectors( point, this.origin ).dot( this.direction );
      // point behind the ray
      if ( directionDistance < 0 ) {
        return this.origin.distanceToSquared( point );
      }
      v1.copy( this.direction ).multiplyScalar( directionDistance ).add( this.origin );
      return v1.distanceToSquared( point );
    };
  }(),

  distanceSqToSegment: function () {
    var segCenter = new Vector3();
    var segDir = new Vector3();
    var diff = new Vector3();

    return function distanceSqToSegment( v0, v1, optionalPointOnRay, optionalPointOnSegment ) {
      // from http://www.geometrictools.com/GTEngine/Include/Mathematics/GteDistRaySegment.h
      // It returns the min distance between the ray and the segment
      // defined by v0 and v1
      // It can also set two optional targets :
      // - The closest point on the ray
      // - The closest point on the segment

      segCenter.copy( v0 ).add( v1 ).multiplyScalar( 0.5 );
      segDir.copy( v1 ).sub( v0 ).normalize();
      diff.copy( this.origin ).sub( segCenter );

      var segExtent = v0.distanceTo( v1 ) * 0.5;
      var a01 = - this.direction.dot( segDir );
      var b0 = diff.dot( this.direction );
      var b1 = - diff.dot( segDir );
      var c = diff.lengthSq();
      var det = Math.abs( 1 - a01 * a01 );
      var s0, s1, sqrDist, extDet;

      if ( det > 0 ) {
        // The ray and segment are not parallel.

        s0 = a01 * b1 - b0;
        s1 = a01 * b0 - b1;
        extDet = segExtent * det;

        if ( s0 >= 0 ) {
          if ( s1 >= - extDet ) {
            if ( s1 <= extDet ) {
              // region 0
              // Minimum at interior points of ray and segment.

              var invDet = 1 / det;
              s0 *= invDet;
              s1 *= invDet;
              sqrDist = s0 * ( s0 + a01 * s1 + 2 * b0 ) + s1 * ( a01 * s0 + s1 + 2 * b1 ) + c;
            } else {
              // region 1

              s1 = segExtent;
              s0 = Math.max( 0, - ( a01 * s1 + b0 ) );
              sqrDist = - s0 * s0 + s1 * ( s1 + 2 * b1 ) + c;
            }
          } else {
            // region 5

            s1 = - segExtent;
            s0 = Math.max( 0, - ( a01 * s1 + b0 ) );
            sqrDist = - s0 * s0 + s1 * ( s1 + 2 * b1 ) + c;
          }
        } else {
          if ( s1 <= - extDet ) {
            // region 4

            s0 = Math.max( 0, - ( - a01 * segExtent + b0 ) );
            s1 = ( s0 > 0 ) ? - segExtent : Math.min( Math.max( - segExtent, - b1 ), segExtent );
            sqrDist = - s0 * s0 + s1 * ( s1 + 2 * b1 ) + c;
          } else if ( s1 <= extDet ) {
            // region 3

            s0 = 0;
            s1 = Math.min( Math.max( - segExtent, - b1 ), segExtent );
            sqrDist = s1 * ( s1 + 2 * b1 ) + c;
          } else {
            // region 2

            s0 = Math.max( 0, - ( a01 * segExtent + b0 ) );
            s1 = ( s0 > 0 ) ? segExtent : Math.min( Math.max( - segExtent, - b1 ), segExtent );
            sqrDist = - s0 * s0 + s1 * ( s1 + 2 * b1 ) + c;
          }
        }
      } else {
        // Ray and segment are parallel.

        s1 = ( a01 > 0 ) ? - segExtent : segExtent;
        s0 = Math.max( 0, - ( a01 * s1 + b0 ) );
        sqrDist = - s0 * s0 + s1 * ( s1 + 2 * b1 ) + c;
      }

      if ( optionalPointOnRay ) {
        optionalPointOnRay.copy( this.direction ).multiplyScalar( s0 ).add( this.origin );
      }

      if ( optionalPointOnSegment ) {
        optionalPointOnSegment.copy( segDir ).multiplyScalar( s1 ).add( segCenter );
      }

      return sqrDist;
    };
  }(),

  applyMatrix4: function ( matrix4 ) {
    this.direction.add( this.origin ).applyMatrix4( matrix4 );
    this.origin.applyMatrix4( matrix4 );
    this.direction.sub( this.origin );
    this.direction.normalize();

    return this;
  },

};

/**
* @author mrdoob / http://mrdoob.com/
* @author mikael emtinger / http://gomo.se/
* @author alteredq / http://alteredqualia.com/
* @author WestLangley / http://github.com/WestLangley
* @author elephantatwork / www.elephantatwork.ch
*/

var object3DId = 0;

function Object3D() {
  Object.defineProperty( this, 'id', { value: object3DId ++ } );

  this.uuid = _Math.generateUUID();

  this.name = '';
  this.type = 'Object3D';

  this.parent = null;
  this.children = [];

  this.up = Object3D.DefaultUp.clone();

  var position = new Vector3();
  var quaternion = new Quaternion();
  var scale = new Vector3( 1, 1, 1 );


  Object.defineProperties( this, {
    position: {
      enumerable: true,
      value: position,
    },
    quaternion: {
      enumerable: true,
      value: quaternion,
    },
    scale: {
      enumerable: true,
      value: scale,
    },
    modelViewMatrix: {
      value: new Matrix4(),
    },
  } );

  this.matrix = new Matrix4();
  this.matrixWorld = new Matrix4();

  this.matrixAutoUpdate = Object3D.DefaultMatrixAutoUpdate;
  this.matrixWorldNeedsUpdate = false;

  this.visible = true;

  this.frustumCulled = true;
  this.renderOrder = 0;

  this.userData = {};
}

Object3D.DefaultUp = new Vector3( 0, 1, 0 );
Object3D.DefaultMatrixAutoUpdate = true;

Object.assign( Object3D.prototype, EventDispatcher.prototype, {

  isObject3D: true,

  add: function ( object ) {
    if ( arguments.length > 1 ) {
      for ( var i = 0; i < arguments.length; i ++ ) {
        this.add( arguments[i] );
      }

      return this;
    }

    if ( ( object && object.isObject3D ) ) {
      if ( object.parent !== null ) {
        object.parent.remove( object );
      }

      object.parent = this;
      object.dispatchEvent( { type: 'added' } );

      this.children.push( object );
    }

    return this;
  },

  remove: function ( object ) {
    if ( arguments.length > 1 ) {
      for ( var i = 0; i < arguments.length; i ++ ) {
        this.remove( arguments[i] );
      }
    }

    var index = this.children.indexOf( object );

    if ( index !== - 1 ) {
      object.parent = null;

      object.dispatchEvent( { type: 'removed' } );

      this.children.splice( index, 1 );
    }
  },

  updateMatrix: function () {
    this.matrix.compose( this.position, this.quaternion, this.scale );

    this.matrixWorldNeedsUpdate = true;
  },

  updateMatrixWorld: function ( force ) {
    if ( this.matrixAutoUpdate === true ) { this.updateMatrix(); }

    if ( this.matrixWorldNeedsUpdate === true || force === true ) {
      if ( this.parent === null ) {
        this.matrixWorld.copy( this.matrix );
      } else {
        this.matrixWorld.multiplyMatrices( this.parent.matrixWorld, this.matrix );
      }

      this.matrixWorldNeedsUpdate = false;

      force = true;
    }

    // update children

    var children = this.children;

    for ( var i = 0, l = children.length; i < l; i ++ ) {
      children[i].updateMatrixWorld( force );
    }
  },
} );


/**
* @author mrdoob / http://mrdoob.com/
*/

function BufferAttribute( array, itemSize, normalized ) {
  if ( Array.isArray( array ) ) {
    throw new TypeError( 'BufferAttribute: array should be a Typed Array.' );
  }

  this.uuid = _Math.generateUUID();

  this.array = array;
  this.itemSize = itemSize;
  this.count = array !== undefined ? array.length / itemSize : 0;
  this.normalized = normalized === true;

  this.dynamic = false;
  this.updateRange = { offset: 0, count: - 1 };

  this.onUploadCallback = function () {};

  this.version = 0;
}

BufferAttribute.prototype = {
  constructor: BufferAttribute,
  isBufferAttribute: true,
};


var count = 0;
function GeometryIdCount() { return count++; }

/**
* @author alteredq / http://alteredqualia.com/
* @author mrdoob / http://mrdoob.com/
*/

function BufferGeometry() {
  Object.defineProperty( this, 'id', { value: GeometryIdCount() } );

  this.uuid = _Math.generateUUID();

  this.name = '';
  this.type = 'BufferGeometry';

  this.index = null;
  this.attributes = {};

  this.groups = [];

  this.boundingBox = null;
  this.boundingSphere = null;

  this.drawRange = { start: 0, count: Infinity };
}

Object.assign( BufferGeometry.prototype, EventDispatcher.prototype, {

  isBufferGeometry: true,

  setIndex: function ( index ) {
    this.index = index;
  },

  addAttribute: function ( name, attribute ) {
    this.attributes[name] = attribute;
    return this;
  },

  dispose: function () {
    this.dispatchEvent( { type: 'dispose' } );
  },

} );

BufferGeometry.MaxIndex = 65535;

/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
* @author mikael emtinger / http://gomo.se/
* @author jonobr1 / http://jonobr1.com/
*/

function Mesh( geometry, material ) {
  Object3D.call( this );

  this.type = 'Mesh';

  this.geometry = geometry !== undefined ? geometry : new BufferGeometry();
  this.material = material;

  this.drawMode = TrianglesDrawMode;
}

Mesh.prototype = Object.assign( Object.create( Object3D.prototype ), {

  constructor: Mesh,

  isMesh: true,
} );


/**
* @author mrdoob / http://mrdoob.com/
* @author mikael emtinger / http://gomo.se/
* @author WestLangley / http://github.com/WestLangley
*/

function Camera() {
  Object3D.call( this );

  this.type = 'Camera';

  this.matrixWorldInverse = new Matrix4();
  this.projectionMatrix = new Matrix4();
}

Camera.prototype = Object.create( Object3D.prototype );
Camera.prototype.constructor = Camera;

Camera.prototype.isCamera = true;

Camera.prototype.lookAt = function () {
// This routine does not support cameras with rotated and/or translated parent(s)

var m1 = new Matrix4();

return function lookAt( vector ) {
  m1.lookAt( this.position, vector, this.up );

  this.quaternion.setFromRotationMatrix( m1 );
};
}();

/**
* @author alteredq / http://alteredqualia.com/
* @author arose / http://github.com/arose
*/

function OrthographicCamera( left, right, top, bottom, near, far ) {
  Camera.call( this );

  this.type = 'OrthographicCamera';

  this.zoom = 1;

  this.left = left;
  this.right = right;
  this.top = top;
  this.bottom = bottom;

  this.near = ( near !== undefined ) ? near : 0.1;
  this.far = ( far !== undefined ) ? far : 2000;

  this.updateProjectionMatrix();
}

OrthographicCamera.prototype = Object.assign( Object.create( Camera.prototype ), {

  constructor: OrthographicCamera,

  updateProjectionMatrix: function () {
    var dx = ( this.right - this.left ) / ( 2 * this.zoom );
    var dy = ( this.top - this.bottom ) / ( 2 * this.zoom );
    var cx = ( this.right + this.left ) / 2;
    var cy = ( this.top + this.bottom ) / 2;

    var left = cx - dx;
    var right = cx + dx;
    var top = cy + dy;
    var bottom = cy - dy;

    this.projectionMatrix.makeOrthographic( left, right, top, bottom, this.near, this.far );
  },
} );

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLIndexedBufferRenderer( gl, extensions, infoRender ) {
  var mode;

  function setMode( value ) {
    mode = value;
  }

  var type, size;

  function setIndex( index ) {
    if ( index.array instanceof Uint32Array && extensions.get( 'OES_element_index_uint' ) ) {
      type = gl.UNSIGNED_INT;
      size = 4;
    } else if ( index.array instanceof Uint16Array ) {
      type = gl.UNSIGNED_SHORT;
      size = 2;
    } else {
      type = gl.UNSIGNED_BYTE;
      size = 1;
    }
  }

  function render( start, count ) {
    gl.drawElements( mode, count, type, start * size );

    infoRender.calls ++;
    infoRender.vertices += count;

    if ( mode === gl.TRIANGLES ) { infoRender.faces += count / 3; }
  }

  return {
    setMode: setMode,
    setIndex: setIndex,
    render: render,
  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLBufferRenderer( gl, extensions, infoRender ) {
  var mode;

  function setMode( value ) {
    mode = value;
  }

  function render( start, count ) {
    gl.drawArrays( mode, start, count );

    infoRender.calls ++;
    infoRender.vertices += count;

    if ( mode === gl.TRIANGLES ) { infoRender.faces += count / 3; }
  }

  return {
    setMode: setMode,
    render: render,
  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/


function WebGLShader( gl, type, string ) {
  var shader = gl.createShader( type );

  gl.shaderSource( shader, string );
  gl.compileShader( shader );

  if ( gl.getShaderParameter( shader, gl.COMPILE_STATUS ) === false ) {
    console.error( 'WebGLShader: Shader couldn\'t compile.' );
  }

  if ( gl.getShaderInfoLog( shader ) !== '' ) {
    var info = gl.getShaderInfoLog( shader );
    // workaround for https://github.com/mrdoob/three.js/issues/9716
    if (info.indexOf('GL_ARB_gpu_shader5') === -1) {
      console.warn( 'WebGLShader: gl.getShaderInfoLog()', type === gl.VERTEX_SHADER ? 'vertex' : 'fragment', info, string );
    }
  }

  return shader;
}

/**
* @author mrdoob / http://mrdoob.com/
*/

var programIdCount = 0;

function generateExtensions( extensions, parameters, rendererExtensions ) {
  extensions = extensions || {};

  var chunks = [
  ( extensions.fragDepth ) && rendererExtensions.get( 'EXT_frag_depth' ) ? '#extension GL_EXT_frag_depth : enable' : '' ];

  return chunks.join( '\n' );
}

function fetchAttributeLocations( gl, program, identifiers ) {
  var attributes = {};

  var n = gl.getProgramParameter( program, gl.ACTIVE_ATTRIBUTES );

  for ( var i = 0; i < n; i ++ ) {
    var info = gl.getActiveAttrib( program, i );
    var name = info.name;

    // console.log("WebGLProgram: ACTIVE VERTEX ATTRIBUTE:", name, i );

    attributes[name] = gl.getAttribLocation( program, name );
  }

  return attributes;
}

function WebGLProgram( renderer, code, material, parameters ) {
  var gl = renderer.context;

  var extensions = material.extensions;

  var vertexShader = material.__webglShader.vertexShader;
  var fragmentShader = material.__webglShader.fragmentShader;

  // console.log( 'building new program ' );

  //

  var customExtensions = generateExtensions( extensions, parameters, renderer.extensions );

  //

  var program = gl.createProgram();

  var prefixVertex, prefixFragment;

  prefixVertex = [
    'precision ' + parameters.precision + ' float;',
    'precision ' + parameters.precision + ' int;',
    '#define SHADER_NAME ' + material.__webglShader.name,
    'uniform mat4 modelMatrix;',
    'uniform mat4 modelViewMatrix;',
    'uniform mat4 projectionMatrix;',
    'uniform mat4 viewMatrix;',
    'attribute vec3 position;',
    '' ].join( '\n' );

  prefixFragment = [
    customExtensions,
    'precision ' + parameters.precision + ' float;',
    'precision ' + parameters.precision + ' int;',
    '#define SHADER_NAME ' + material.__webglShader.name,
    ( parameters.useFog && parameters.fog ) ? '#define USE_FOG' : '',
    '' ].join( '\n' );

  var vertexGlsl = prefixVertex + vertexShader;
  var fragmentGlsl = prefixFragment + fragmentShader;

  // console.log( '*VERTEX*', vertexGlsl );
  // console.log( '*FRAGMENT*', fragmentGlsl );

  var glVertexShader = WebGLShader( gl, gl.VERTEX_SHADER, vertexGlsl );
  var glFragmentShader = WebGLShader( gl, gl.FRAGMENT_SHADER, fragmentGlsl );

  gl.attachShader( program, glVertexShader );
  gl.attachShader( program, glFragmentShader );

  gl.linkProgram( program );

  var programLog = gl.getProgramInfoLog( program );
  var vertexLog = gl.getShaderInfoLog( glVertexShader );
  var fragmentLog = gl.getShaderInfoLog( glFragmentShader );

  // console.log( '**VERTEX**', gl.getExtension( 'WEBGL_debug_shaders' ).getTranslatedShaderSource( glVertexShader ) );
  // console.log( '**FRAGMENT**', gl.getExtension( 'WEBGL_debug_shaders' ).getTranslatedShaderSource( glFragmentShader ) );

  if ( gl.getProgramParameter( program, gl.LINK_STATUS ) === false ) {
    console.error( 'WebGLProgram: shader error: ', gl.getError(), 'gl.VALIDATE_STATUS', gl.getProgramParameter( program, gl.VALIDATE_STATUS ), 'gl.getProgramInfoLog', programLog, vertexLog, fragmentLog );
  } else if ( programLog !== '' ) {
    console.warn( 'WebGLProgram: gl.getProgramInfoLog()', programLog );
  }

  // clean up

  gl.deleteShader( glVertexShader );
  gl.deleteShader( glFragmentShader );

  // set up caching for uniform locations

  var cachedUniforms;

  this.getUniforms = function () {
    if ( cachedUniforms === undefined ) {
      cachedUniforms =
      new WebGLUniforms( gl, program, renderer );
    }

    return cachedUniforms;
  };

  // set up caching for attribute locations

  var cachedAttributes;

  this.getAttributes = function () {
    if ( cachedAttributes === undefined ) {
      cachedAttributes = fetchAttributeLocations( gl, program );
    }

    return cachedAttributes;
  };

  // free resource

  this.destroy = function () {
    gl.deleteProgram( program );
    this.program = undefined;
  };

  //

  this.id = programIdCount ++;
  this.code = code;
  this.usedTimes = 1;
  this.program = program;
  this.vertexShader = glVertexShader;
  this.fragmentShader = glFragmentShader;

  return this;
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLPrograms( renderer, capabilities ) {
  var programs = [];

  var parameterNames = [
    'precision',
    'fog', 'useFog',
    'premultipliedAlpha' ];

  this.getParameters = function ( material, fog, object ) {
    var precision = renderer.getPrecision();

    if ( material.precision !== null ) {
      precision = capabilities.getMaxPrecision( material.precision );

      if ( precision !== material.precision ) {
        console.warn( 'WebGLProgram.getParameters:', material.precision, 'not supported, using', precision, 'instead.' );
      }
    }

    var parameters = {
      precision: precision,
      fog: !! fog,
      useFog: material.fog,
      premultipliedAlpha: material.premultipliedAlpha,
    };

    return parameters;
  };

  this.getProgramCode = function ( material, parameters ) {
    var array = [];

    array.push( material.fragmentShader );
    array.push( material.vertexShader );

    for ( var i = 0; i < parameterNames.length; i ++ ) {
      array.push( parameters[parameterNames[i]] );
    }

    return array.join();
  };

  this.acquireProgram = function ( material, parameters, code ) {
    var program;

    // Check if code has been already compiled
    for ( var p = 0, pl = programs.length; p < pl; p ++ ) {
      var programInfo = programs[p];

      if ( programInfo.code === code ) {
        program = programInfo;
        ++ program.usedTimes;

        break;
      }
    }

    if ( program === undefined ) {
      program = new WebGLProgram( renderer, code, material, parameters );
      programs.push( program );
    }

    return program;
  };

  this.releaseProgram = function ( program ) {
    if ( -- program.usedTimes === 0 ) {
      // Remove from unordered set
      var i = programs.indexOf( program );
      programs[i] = programs[programs.length - 1];
      programs.pop();

      // Free WebGL resources
      program.destroy();
    }
  };

  // Exposed for resource monitoring & error feedback via renderer.info:
  this.programs = programs;
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLGeometries( gl, properties, info ) {
  var geometries = {};

  function onGeometryDispose( event ) {
    var geometry = event.target;
    var buffergeometry = geometries[geometry.id];

    if ( buffergeometry.index !== null ) {
      deleteAttribute( buffergeometry.index );
    }

    deleteAttributes( buffergeometry.attributes );

    geometry.removeEventListener( 'dispose', onGeometryDispose );

    delete geometries[geometry.id];

    properties.delete( geometry );

    properties.delete( buffergeometry );
  }

  function getAttributeBuffer( attribute ) {
    return properties.get( attribute ).__webglBuffer;
  }

  function deleteAttribute( attribute ) {
    var buffer = getAttributeBuffer( attribute );

    if ( buffer !== undefined ) {
      gl.deleteBuffer( buffer );
      removeAttributeBuffer( attribute );
    }
  }

  function deleteAttributes( attributes ) {
    for ( var name in attributes ) {
      deleteAttribute( attributes[name] );
    }
  }

  function removeAttributeBuffer( attribute ) {
    properties.delete( attribute );
  }

  return {

    get: function ( object ) {
      var geometry = object.geometry;

      if ( geometries[geometry.id] !== undefined ) {
        return geometries[geometry.id];
      }

      geometry.addEventListener( 'dispose', onGeometryDispose );

      var buffergeometry;

      if ( geometry.isBufferGeometry ) {
        buffergeometry = geometry;
      }

      geometries[geometry.id] = buffergeometry;

      return buffergeometry;
    },

  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLObjects( gl, properties, info ) {
  var geometries = new WebGLGeometries( gl, properties, info );

  //

  function update( object ) {
    var geometry = geometries.get( object );

    var index = geometry.index;
    var attributes = geometry.attributes;

    if ( index !== null ) {
      updateAttribute( index, gl.ELEMENT_ARRAY_BUFFER );
    }

    for ( var name in attributes ) {
      updateAttribute( attributes[name], gl.ARRAY_BUFFER );
    }

    return geometry;
  }

  function updateAttribute( attribute, bufferType ) {
    var data = attribute;

    var attributeProperties = properties.get( data );

    if ( attributeProperties.__webglBuffer === undefined ) {
      createBuffer( attributeProperties, data, bufferType );
    } else if ( attributeProperties.version !== data.version ) {
      updateBuffer( attributeProperties, data, bufferType );
    }
  }

  function createBuffer( attributeProperties, data, bufferType ) {
    attributeProperties.__webglBuffer = gl.createBuffer();
    gl.bindBuffer( bufferType, attributeProperties.__webglBuffer );

    var usage = data.dynamic ? gl.DYNAMIC_DRAW : gl.STATIC_DRAW;

    gl.bufferData( bufferType, data.array, usage );

    var type = gl.FLOAT;
    var array = data.array;

    if ( array instanceof Float32Array ) {
      type = gl.FLOAT;
    } else if ( array instanceof Float64Array ) {
      console.warn( 'Unsupported data buffer format: Float64Array' );
    } else if ( array instanceof Uint16Array ) {
      type = gl.UNSIGNED_SHORT;
    } else if ( array instanceof Int16Array ) {
      type = gl.SHORT;
    } else if ( array instanceof Uint32Array ) {
      type = gl.UNSIGNED_INT;
    } else if ( array instanceof Int32Array ) {
      type = gl.INT;
    } else if ( array instanceof Int8Array ) {
      type = gl.BYTE;
    } else if ( array instanceof Uint8Array ) {
      type = gl.UNSIGNED_BYTE;
    }

    attributeProperties.bytesPerElement = array.BYTES_PER_ELEMENT;
    attributeProperties.type = type;
    attributeProperties.version = data.version;

    data.onUploadCallback();
  }

  function updateBuffer( attributeProperties, data, bufferType ) {
    gl.bindBuffer( bufferType, attributeProperties.__webglBuffer );

    if ( data.dynamic === false ) {
      gl.bufferData( bufferType, data.array, gl.STATIC_DRAW );
    } else if ( data.updateRange.count === - 1 ) {
      // Not using update ranges

      gl.bufferSubData( bufferType, 0, data.array );
    } else if ( data.updateRange.count === 0 ) {
      console.error( 'WebGLObjects.updateBuffer: updateRange.count is 0.' );
    } else {
      gl.bufferSubData( bufferType, data.updateRange.offset * data.array.BYTES_PER_ELEMENT,
                        data.array.subarray( data.updateRange.offset, data.updateRange.offset + data.updateRange.count ) );

      data.updateRange.count = 0; // reset range
    }

    attributeProperties.version = data.version;
  }

  function getAttributeBuffer( attribute ) {
    return properties.get( attribute ).__webglBuffer;
  }

  function getAttributeProperties( attribute ) {
    return properties.get( attribute );
  }


  return {

    getAttributeBuffer: getAttributeBuffer,
    getAttributeProperties: getAttributeProperties,

    update: update,

  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLTextures( _gl, extensions, state, properties, capabilities, info ) {
  function onTextureDispose( event ) {
    var texture = event.target;

    texture.removeEventListener( 'dispose', onTextureDispose );

    deallocateTexture( texture );
  }

  //

  function deallocateTexture( texture ) {
    var textureProperties = properties.get( texture );

    // 2D texture

    if ( textureProperties.__webglInit === undefined ) { return; }

    _gl.deleteTexture( textureProperties.__webglTexture );

    // remove all webgl properties
    properties.delete( texture );
  }


  function setTexture2D( texture, slot ) {
    var textureProperties = properties.get( texture );

    if ( texture.version > 0 && textureProperties.__version !== texture.version ) {
      var image = texture.image;

      if ( image === undefined ) {
        console.warn( 'WebGLRenderer: Texture marked for update but image is undefined', texture );
      } else if ( image.complete === false ) {
        console.warn( 'WebGLRenderer: Texture marked for update but image is incomplete', texture );
      } else {
        uploadTexture( textureProperties, texture, slot );
        return;
      }
    }

    state.activeTexture( _gl.TEXTURE0 + slot );
    state.bindTexture( _gl.TEXTURE_2D, textureProperties.__webglTexture );
  }

  function setTextureParameters( textureType ) {
    _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_S, _gl.CLAMP_TO_EDGE );
    _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_T, _gl.CLAMP_TO_EDGE );
    _gl.texParameteri( textureType, _gl.TEXTURE_MAG_FILTER, _gl.LINEAR );
    _gl.texParameteri( textureType, _gl.TEXTURE_MIN_FILTER, _gl.LINEAR_MIPMAP_LINEAR );
  }

  function uploadTexture( textureProperties, texture, slot ) {
    if ( textureProperties.__webglInit === undefined ) {
      textureProperties.__webglInit = true;

      texture.addEventListener( 'dispose', onTextureDispose );

      textureProperties.__webglTexture = _gl.createTexture();
    }

    state.activeTexture( _gl.TEXTURE0 + slot );
    state.bindTexture( _gl.TEXTURE_2D, textureProperties.__webglTexture );

    _gl.pixelStorei( _gl.UNPACK_FLIP_Y_WEBGL, true );
    _gl.pixelStorei( _gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, false );
    _gl.pixelStorei( _gl.UNPACK_ALIGNMENT, 4 );

    var image = texture.image;

    var glFormat = _gl.RGBA;
    var glType = _gl.UNSIGNED_BYTE;

    setTextureParameters( _gl.TEXTURE_2D );

    state.texImage2D( _gl.TEXTURE_2D, 0, glFormat, glFormat, glType, image );

    _gl.generateMipmap( _gl.TEXTURE_2D );

    textureProperties.__version = texture.version;
  }

  this.setTexture2D = setTexture2D;
}

/**
* @author fordacious / fordacious.github.io
*/

function WebGLProperties() {
  var properties = {};

  return {

    get: function ( object ) {
      var uuid = object.uuid;
      var map = properties[uuid];

      if ( map === undefined ) {
        map = {};
        properties[uuid] = map;
      }

      return map;
    },

    delete: function ( object ) {
      delete properties[object.uuid];
    },

  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLState( gl, extensions ) {
  function ColorBuffer() {
    var color = new Vector4();
    var currentColorClear = new Vector4();

    return {

      setClear: function ( r, g, b, a, premultipliedAlpha ) {
        if ( premultipliedAlpha === true ) {
          r *= a; g *= a; b *= a;
        }

        color.set( r, g, b, a );

        if ( currentColorClear.equals( color ) === false ) {
          gl.clearColor( r, g, b, a );
          currentColorClear.copy( color );
        }
      },

    };
  }

  function DepthBuffer() {
    var currentDepthMask = null;
    var currentDepthClear = null;

    return {

      setTest: function ( depthTest ) {
        if ( depthTest ) {
          enable( gl.DEPTH_TEST );
        } else {
          disable( gl.DEPTH_TEST );
        }
      },

      setMask: function ( depthMask ) {
        if ( currentDepthMask !== depthMask ) {
          gl.depthMask( depthMask );
          currentDepthMask = depthMask;
        }
      },

      setClear: function ( depth ) {
        if ( currentDepthClear !== depth ) {
          gl.clearDepth( depth );
          currentDepthClear = depth;
        }
      },

      reset: function () {
        currentDepthMask = null;
        currentDepthClear = null;
      },

    };
  }


  //

  var colorBuffer = new ColorBuffer();
  var depthBuffer = new DepthBuffer();

  var maxVertexAttributes = gl.getParameter( gl.MAX_VERTEX_ATTRIBS );
  var newAttributes = new Uint8Array( maxVertexAttributes );
  var enabledAttributes = new Uint8Array( maxVertexAttributes );

  var capabilities = {};

  var currentBlending = null;
  var currentPremultipledAlpha = false;

  var currentLineWidth = null;

  var version = parseFloat( /^WebGL\ ([0-9])/.exec( gl.getParameter( gl.VERSION ) )[1] );
  var lineWidthAvailable = parseFloat( version ) >= 1.0;

  var currentTextureSlot = null;
  var currentBoundTextures = {};

  var currentViewport = new Vector4();

  function createTexture( type, target, count ) {
    var data = new Uint8Array( 4 ); // 4 is required to match default unpack alignment of 4.
    var texture = gl.createTexture();

    gl.bindTexture( type, texture );
    gl.texParameteri( type, gl.TEXTURE_MIN_FILTER, gl.NEAREST );
    gl.texParameteri( type, gl.TEXTURE_MAG_FILTER, gl.NEAREST );

    for ( var i = 0; i < count; i ++ ) {
      gl.texImage2D( target + i, 0, gl.RGBA, 1, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, data );
    }

    return texture;
  }

  var emptyTextures = {};
  emptyTextures[gl.TEXTURE_2D] = createTexture( gl.TEXTURE_2D, gl.TEXTURE_2D, 1 );

  //

  function init() {
    colorBuffer.setClear( 0, 0, 0, 1 );
    depthBuffer.setClear( 1 );

    enable( gl.DEPTH_TEST );
    gl.depthFunc( gl.LEQUAL );

    enable( gl.BLEND );
    setBlending( NormalBlending );
  }

  function initAttributes() {
    for ( var i = 0, l = newAttributes.length; i < l; i ++ ) {
      newAttributes[i] = 0;
    }
  }

  function enableAttribute( attribute ) {
    newAttributes[attribute] = 1;

    if ( enabledAttributes[attribute] === 0 ) {
      gl.enableVertexAttribArray( attribute );
      enabledAttributes[attribute] = 1;
    }
  }

  function disableUnusedAttributes() {
    for ( var i = 0, l = enabledAttributes.length; i !== l; ++ i ) {
      if ( enabledAttributes[i] !== newAttributes[i] ) {
        gl.disableVertexAttribArray( i );
        enabledAttributes[i] = 0;
      }
    }
  }

  function enable( id ) {
    if ( capabilities[id] !== true ) {
      gl.enable( id );
      capabilities[id] = true;
    }
  }

  function disable( id ) {
    if ( capabilities[id] !== false ) {
      gl.disable( id );
      capabilities[id] = false;
    }
  }

  function setBlending( blending, premultipliedAlpha ) {
    if ( blending !== NoBlending ) {
      enable( gl.BLEND );
    } else {
      disable( gl.BLEND );
    }

    if ( blending !== currentBlending || premultipliedAlpha !== currentPremultipledAlpha ) {
      if ( premultipliedAlpha ) {
        gl.blendEquationSeparate( gl.FUNC_ADD, gl.FUNC_ADD );
        gl.blendFuncSeparate( gl.ONE, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA );
      } else {
        gl.blendEquationSeparate( gl.FUNC_ADD, gl.FUNC_ADD );
        gl.blendFuncSeparate( gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA );
      }

      currentBlending = blending;
      currentPremultipledAlpha = premultipliedAlpha;
    }
  }

  function setDepthTest( depthTest ) {
    depthBuffer.setTest( depthTest );
  }

  function setDepthWrite( depthWrite ) {
    depthBuffer.setMask( depthWrite );
  }

  //

  function setLineWidth( width ) {
    if ( width !== currentLineWidth ) {
      if ( lineWidthAvailable ) { gl.lineWidth( width ); }

      currentLineWidth = width;
    }
  }

  // texture

  function activeTexture( webglSlot ) {
    if ( currentTextureSlot !== webglSlot ) {
      gl.activeTexture( webglSlot );
      currentTextureSlot = webglSlot;
    }
  }

  function bindTexture( webglType, webglTexture ) {
    var boundTexture = currentBoundTextures[currentTextureSlot];
    if ( boundTexture === undefined ) {
      boundTexture = { type: undefined, texture: undefined };
      currentBoundTextures[currentTextureSlot] = boundTexture;
    }

    if ( boundTexture.type !== webglType || boundTexture.texture !== webglTexture ) {
      gl.bindTexture( webglType, webglTexture || emptyTextures[webglType] );

      boundTexture.type = webglType;
      boundTexture.texture = webglTexture;
    }
  }

  function texImage2D() {
    try {
      gl.texImage2D.apply( gl, arguments );
    } catch ( error ) {
      console.error( error );
    }
  }

  //

  function viewport( viewport ) {
    if ( currentViewport.equals( viewport ) === false ) {
      gl.viewport( viewport.x, viewport.y, viewport.z, viewport.w );
      currentViewport.copy( viewport );
    }
  }

  //

  return {
    buffers: {
      color: colorBuffer,
      depth: depthBuffer,
    },

    init: init,
    initAttributes: initAttributes,
    enableAttribute: enableAttribute,
    disableUnusedAttributes: disableUnusedAttributes,
    enable: enable,
    disable: disable,

    setBlending: setBlending,

    setDepthTest: setDepthTest,
    setDepthWrite: setDepthWrite,

    setLineWidth: setLineWidth,

    activeTexture: activeTexture,
    bindTexture: bindTexture,
    texImage2D: texImage2D,

    viewport: viewport,
  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLCapabilities( gl, extensions, parameters ) {
  function getMaxPrecision( precision ) {
    if ( precision === 'highp' ) {
      if ( gl.getShaderPrecisionFormat( gl.VERTEX_SHADER, gl.HIGH_FLOAT ).precision > 0 &&
         gl.getShaderPrecisionFormat( gl.FRAGMENT_SHADER, gl.HIGH_FLOAT ).precision > 0 ) {
        return 'highp';
      }

      precision = 'mediump';
    }

    if ( precision === 'mediump' ) {
      if ( gl.getShaderPrecisionFormat( gl.VERTEX_SHADER, gl.MEDIUM_FLOAT ).precision > 0 &&
         gl.getShaderPrecisionFormat( gl.FRAGMENT_SHADER, gl.MEDIUM_FLOAT ).precision > 0 ) {
        return 'mediump';
      }
    }

    return 'lowp';
  }

  var precision = parameters.precision !== undefined ? parameters.precision : 'highp';
  var maxPrecision = getMaxPrecision( precision );

  if ( maxPrecision !== precision ) {
    console.warn( 'WebGLRenderer:', precision, 'not supported, using', maxPrecision, 'instead.' );
    precision = maxPrecision;
  }

  var maxTextures = gl.getParameter( gl.MAX_TEXTURE_IMAGE_UNITS );

  return {
    getMaxPrecision: getMaxPrecision,
    precision: precision,
    maxTextures: maxTextures,
  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLExtensions( gl ) {
  var extensions = {};

  return {

    get: function ( name ) {
      if ( extensions[name] !== undefined ) {
        return extensions[name];
      }

      var extension;

      switch ( name ) {
        case 'WEBGL_depth_texture':
          extension = gl.getExtension( 'WEBGL_depth_texture' ) || gl.getExtension( 'MOZ_WEBGL_depth_texture' ) || gl.getExtension( 'WEBKIT_WEBGL_depth_texture' );
          break;
        default:
          extension = gl.getExtension( name );
      }

      if ( extension === null ) {
        console.warn( 'WebGLRenderer: ' + name + ' extension not supported.' );
      }

      extensions[name] = extension;

      return extension;
    },

  };
}

/**
* @author supereggbert / http://www.paulbrunt.co.uk/
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
* @author szimek / https://github.com/szimek/
* @author tschw
*/

function WebGLRenderer( parameters ) {
  parameters = parameters || {};

  var _canvas = parameters.canvas !== undefined ? parameters.canvas : document.createElementNS( 'http://www.w3.org/1999/xhtml', 'canvas' ),
    _context = parameters.context !== undefined ? parameters.context : null,
    _antialias = parameters.antialias !== undefined ? parameters.antialias : false,
    _premultipliedAlpha = parameters.premultipliedAlpha !== undefined ? parameters.premultipliedAlpha : true;

  var opaqueObjects = [];
  var opaqueObjectsLastIndex = - 1;
  var transparentObjects = [];
  var transparentObjectsLastIndex = - 1;

  // public properties

  this.domElement = _canvas;
  this.context = null;

  // clearing

  this.autoClear = true;
  this.autoClearColor = true;
  this.autoClearDepth = true;

  // scene graph

  this.sortObjects = true;

  // internal properties

  var _this = this,

    // internal state cache

    _currentProgram = null,
    _currentRenderTarget = null,
    _currentFramebuffer = null,
    _currentMaterialId = - 1,
    _currentGeometryProgram = '',
    _currentCamera = null,

    _currentViewport = new Vector4(),

    //

    _usedTextureUnits = 0,

    //

    _clearColor = new Color( 0x000000 ),
    _clearAlpha = 0,

    _width = _canvas.width,
    _height = _canvas.height,

    _pixelRatio = 1,

    _viewport = new Vector4( 0, 0, _width, _height ),

    // camera matrices cache

    _projScreenMatrix = new Matrix4(),

    _vector3 = new Vector3(),

    // info

    _infoRender = {
      calls: 0,
      vertices: 0,
      faces: 0,
      points: 0,
    };

  this.info = {
    render: _infoRender,
    programs: null,
  };


  // initialize

  var _gl;

  try {
    var attributes = {
      alpha: false,
      depth: true,
      antialias: _antialias,
      premultipliedAlpha: _premultipliedAlpha,
      preserveDrawingBuffer: false,
    };

    _gl = _context || _canvas.getContext( 'webgl', attributes ) || _canvas.getContext( 'experimental-webgl', attributes );

    if ( _gl === null ) {
      if ( _canvas.getContext( 'webgl' ) !== null ) {
        throw new Error('Error creating WebGL context with your selected attributes.');
      } else {
        throw new Error('Error creating WebGL context.');
      }
    }

    // Some experimental-webgl implementations do not have getShaderPrecisionFormat

    if ( _gl.getShaderPrecisionFormat === undefined ) {
      _gl.getShaderPrecisionFormat = function () {
        return { 'rangeMin': 1, 'rangeMax': 1, 'precision': 1 };
      };
    }

    _canvas.addEventListener( 'webglcontextlost', onContextLost, false );
  } catch ( error ) {
    console.error( 'WebGLRenderer: ' + error );
  }

  var extensions = new WebGLExtensions( _gl );
  extensions.get( 'WEBGL_depth_texture' );
  if ( extensions.get( 'OES_element_index_uint' ) ) {
    BufferGeometry.MaxIndex = 4294967296;
  }

  var capabilities = new WebGLCapabilities( _gl, extensions, parameters );

  var state = new WebGLState( _gl, extensions );
  var properties = new WebGLProperties();
  var textures = new WebGLTextures( _gl, extensions, state, properties, capabilities, this.info );
  var objects = new WebGLObjects( _gl, properties, this.info );
  var programCache = new WebGLPrograms( this, capabilities );

  this.info.programs = programCache.programs;

  var bufferRenderer = new WebGLBufferRenderer( _gl, extensions, _infoRender );
  var indexedBufferRenderer = new WebGLIndexedBufferRenderer( _gl, extensions, _infoRender );

  //

  function getTargetPixelRatio() {
    return _currentRenderTarget === null ? _pixelRatio : 1;
  }

  function setDefaultGLState() {
    state.init();

    state.viewport( _currentViewport.copy( _viewport ).multiplyScalar( _pixelRatio ) );

    state.buffers.color.setClear( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha, _premultipliedAlpha );
  }

  function resetGLState() {
    _currentProgram = null;
    _currentCamera = null;

    _currentGeometryProgram = '';
    _currentMaterialId = - 1;

    state.reset();
  }

  setDefaultGLState();

  this.context = _gl;
  this.capabilities = capabilities;
  this.extensions = extensions;
  this.properties = properties;
  this.state = state;

  // API

  this.getContext = function () {
    return _gl;
  };

  this.getPrecision = function () {
    return capabilities.precision;
  };

  this.getPixelRatio = function () {
    return _pixelRatio;
  };

  this.setPixelRatio = function ( value ) {
    if ( value === undefined ) { return; }

    _pixelRatio = value;

    this.setSize( _viewport.z, _viewport.w, false );
  };

  this.setSize = function ( width, height, updateStyle ) {
    _width = width;
    _height = height;

    _canvas.width = width * _pixelRatio;
    _canvas.height = height * _pixelRatio;

    if ( updateStyle !== false ) {
      _canvas.style.width = width + 'px';
      _canvas.style.height = height + 'px';
    }

    this.setViewport( 0, 0, width, height );
  };

  this.setViewport = function ( x, y, width, height ) {
    state.viewport( _viewport.set( x, y, width, height ) );
  };

  // Clearing

  this.getClearColor = function () {
    return _clearColor;
  };

  this.setClearColor = function ( color, alpha ) {
    _clearColor.set( color );

    _clearAlpha = alpha !== undefined ? alpha : 1;

    state.buffers.color.setClear( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha, _premultipliedAlpha );
  };

  this.clear = function ( color, depth ) {
    var bits = 0;

    if ( color === undefined || color ) { bits |= _gl.COLOR_BUFFER_BIT; }
    if ( depth === undefined || depth ) { bits |= _gl.DEPTH_BUFFER_BIT; }

    _gl.clear( bits );
  };

  this.clearColor = function () {
    this.clear( true, false );
  };

  this.clearDepth = function () {
    this.clear( false, true );
  };

  // Reset

  this.resetGLState = resetGLState;

  this.dispose = function () {
    transparentObjects = [];
    transparentObjectsLastIndex = -1;
    opaqueObjects = [];
    opaqueObjectsLastIndex = -1;

    _canvas.removeEventListener( 'webglcontextlost', onContextLost, false );
  };

  // Events

  function onContextLost( event ) {
    event.preventDefault();

    resetGLState();
    setDefaultGLState();

    properties.clear();
  }

  function onMaterialDispose( event ) {
    var material = event.target;

    material.removeEventListener( 'dispose', onMaterialDispose );

    deallocateMaterial( material );
  }

  // Buffer deallocation

  function deallocateMaterial( material ) {
    releaseMaterialProgramReference( material );

    properties.delete( material );
  }


  function releaseMaterialProgramReference( material ) {
    var programInfo = properties.get( material ).program;

    material.program = undefined;

    if ( programInfo !== undefined ) {
      programCache.releaseProgram( programInfo );
    }
  }

  // Buffer rendering


  this.renderBufferDirect = function ( camera, fog, geometry, material, object, group ) {
    setMaterial( material );

    var program = setProgram( camera, fog, material, object );

    var updateBuffers = false;
    var geometryProgram = geometry.id + '_' + program.id;

    if ( geometryProgram !== _currentGeometryProgram ) {
      _currentGeometryProgram = geometryProgram;
      updateBuffers = true;
    }

    //

    var index = geometry.index;
    var position = geometry.attributes.position;
    var rangeFactor = 1;

    var renderer;

    if ( index !== null ) {
      renderer = indexedBufferRenderer;
      renderer.setIndex( index );
    } else {
      renderer = bufferRenderer;
    }

    if ( updateBuffers ) {
      setupVertexAttributes( material, program, geometry );

      if ( index !== null ) {
        _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, objects.getAttributeBuffer( index ) );
      }
    }

    //

    var dataCount = 0;

    if ( index !== null ) {
      dataCount = index.count;
    } else if ( position !== undefined ) {
      dataCount = position.count;
    }

    var rangeStart = geometry.drawRange.start * rangeFactor;
    var rangeCount = geometry.drawRange.count * rangeFactor;

    var groupStart = group !== null ? group.start * rangeFactor : 0;
    var groupCount = group !== null ? group.count * rangeFactor : Infinity;

    var drawStart = Math.max( rangeStart, groupStart );
    var drawEnd = Math.min( dataCount, rangeStart + rangeCount, groupStart + groupCount ) - 1;

    var drawCount = Math.max( 0, drawEnd - drawStart + 1 );

    if ( drawCount === 0 ) { return; }

    //

    if ( object.isMesh ) {
      switch ( object.drawMode ) {
        case TrianglesDrawMode:
          renderer.setMode( _gl.TRIANGLES );
          break;

        case TriangleStripDrawMode:
          renderer.setMode( _gl.TRIANGLE_STRIP );
          break;

        case TriangleFanDrawMode:
          renderer.setMode( _gl.TRIANGLE_FAN );
          break;
      }
    } else if ( object.isLine ) {
      var lineWidth = material.linewidth;
      if ( lineWidth === undefined ) { lineWidth = 1; } // Not using Line*Material
      state.setLineWidth( lineWidth * getTargetPixelRatio() );

      if ( object.isLineSegments ) {
        renderer.setMode( _gl.LINES );
      } else {
        renderer.setMode( _gl.LINE_STRIP );
      }
    } else if ( object.isPoints ) {
      renderer.setMode( _gl.POINTS );
    }

    renderer.render( drawStart, drawCount );
  };

  function setupVertexAttributes( material, program, geometry, startIndex ) {
    if ( startIndex === undefined ) { startIndex = 0; }

    state.initAttributes();

    var geometryAttributes = geometry.attributes;

    var programAttributes = program.getAttributes();

    for ( var name in programAttributes ) {
      var programAttribute = programAttributes[name];

      if ( programAttribute >= 0 ) {
        var geometryAttribute = geometryAttributes[name];

        if ( geometryAttribute !== undefined ) {
          var normalized = geometryAttribute.normalized;
          var size = geometryAttribute.itemSize;

          var attributeProperties = objects.getAttributeProperties( geometryAttribute );

          var buffer = attributeProperties.__webglBuffer;
          var type = attributeProperties.type;
          var bytesPerElement = attributeProperties.bytesPerElement;

          state.enableAttribute( programAttribute );

          _gl.bindBuffer( _gl.ARRAY_BUFFER, buffer );
          _gl.vertexAttribPointer( programAttribute, size, type, normalized, 0, startIndex * size * bytesPerElement );
        } else {
          console.error( 'undefined geometryAttribute' );
        }
      }
    }
    state.disableUnusedAttributes();
  }

  // Sorting

  function painterSortStable( a, b ) {
    if ( a.object.renderOrder !== b.object.renderOrder ) {
      return a.object.renderOrder - b.object.renderOrder;
    } else if ( a.material.program && b.material.program && a.material.program !== b.material.program ) {
      return a.material.program.id - b.material.program.id;
    } else if ( a.material.id !== b.material.id ) {
      return a.material.id - b.material.id;
    } else if ( a.z !== b.z ) {
      return a.z - b.z;
    } else {
      return a.id - b.id;
    }
  }

  function reversePainterSortStable( a, b ) {
    if ( a.object.renderOrder !== b.object.renderOrder ) {
      return a.object.renderOrder - b.object.renderOrder;
    } if ( a.z !== b.z ) {
      return b.z - a.z;
    } else {
      return a.id - b.id;
    }
  }

  // Rendering

  this.render = function ( scene, camera, renderTarget, forceClear ) {
    if ( camera !== undefined && camera.isCamera !== true ) {
      console.error( 'camera is not an instance of Camera.' );
      return;
    }

    // reset caching for this frame

    _currentGeometryProgram = '';
    _currentMaterialId = - 1;
    _currentCamera = null;

    // update scene graph

    if ( scene.autoUpdate === true ) { scene.updateMatrixWorld(); }

    // update camera matrices and frustum

    if ( camera.parent === null ) { camera.updateMatrixWorld(); }

    camera.matrixWorldInverse.getInverse( camera.matrixWorld );

    _projScreenMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse );

    opaqueObjectsLastIndex = - 1;
    transparentObjectsLastIndex = - 1;

    projectObject( scene );

    opaqueObjects.length = opaqueObjectsLastIndex + 1;
    transparentObjects.length = transparentObjectsLastIndex + 1;

    if ( _this.sortObjects === true ) {
      opaqueObjects.sort( painterSortStable );
      transparentObjects.sort( reversePainterSortStable );
    }

    //

    _infoRender.calls = 0;
    _infoRender.vertices = 0;
    _infoRender.faces = 0;
    _infoRender.points = 0;

    if ( renderTarget === undefined ) {
      renderTarget = null;
    }

    this.setRenderTarget( renderTarget );

    state.buffers.color.setClear( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha, _premultipliedAlpha );

    if ( this.autoClear || forceClear ) {
      this.clear( this.autoClearColor, this.autoClearDepth );
    }

    // opaque pass (front-to-back order)

    state.setBlending( NoBlending );
    renderObjects( opaqueObjects, scene, camera );

    // transparent pass (back-to-front order)

    renderObjects( transparentObjects, scene, camera );

    // Generate mipmap if we're using any kind of mipmap filtering

    if ( renderTarget ) {
      textures.updateRenderTargetMipmap( renderTarget );
    }

    // Ensure depth buffer writing is enabled so it can be cleared on next render

    state.setDepthTest( true );
    state.setDepthWrite( true );

    // _gl.finish();
  };

  function pushRenderItem( object, geometry, material, z, group ) {
    var array, index;

    // allocate the next position in the appropriate array

    if ( material.transparent ) {
      array = transparentObjects;
      index = ++ transparentObjectsLastIndex;
    } else {
      array = opaqueObjects;
      index = ++ opaqueObjectsLastIndex;
    }

    // recycle existing render item or grow the array

    var renderItem = array[index];

    if ( renderItem !== undefined ) {
      renderItem.id = object.id;
      renderItem.object = object;
      renderItem.geometry = geometry;
      renderItem.material = material;
      renderItem.z = _vector3.z;
      renderItem.group = group;
    } else {
      renderItem = {
        id: object.id,
        object: object,
        geometry: geometry,
        material: material,
        z: _vector3.z,
        group: group,
      };

      // assert( index === array.length );
      array.push( renderItem );
    }
  }

  function projectObject( object ) {
    if ( object.visible === false ) { return; }

    if ( object.isMesh || object.isLine || object.isPoints ) {
      var material = object.material;

      if ( material.visible === true ) {
        if ( _this.sortObjects === true ) {
          _vector3.setFromMatrixPosition( object.matrixWorld );
          _vector3.applyProjection( _projScreenMatrix );
        }

        var geometry = objects.update( object );

        pushRenderItem( object, geometry, material, _vector3.z, null );
      }
    }

    var children = object.children;

    for ( var i = 0, l = children.length; i < l; i ++ ) {
      projectObject( children[i] );
    }
  }

  function renderObjects( renderList, scene, camera, overrideMaterial ) {
    for ( var i = 0, l = renderList.length; i < l; i ++ ) {
      var renderItem = renderList[i];
      var object = renderItem.object;
      var geometry = renderItem.geometry;
      var material = overrideMaterial === undefined ? renderItem.material : overrideMaterial;
      var group = renderItem.group;

      object.modelViewMatrix.multiplyMatrices( camera.matrixWorldInverse, object.matrixWorld );
      _this.renderBufferDirect( camera, scene.fog, geometry, material, object, group );
    }
  }

  function initMaterial( material, fog, object ) {
    var materialProperties = properties.get( material );

    var parameters = programCache.getParameters(
      material, fog, object );

    var code = programCache.getProgramCode( material, parameters );

    var program = materialProperties.program;
    var programChange = true;

    if ( program === undefined ) {
      // new material
      material.addEventListener( 'dispose', onMaterialDispose );
    } else if ( program.code !== code ) {
      // changed glsl or parameters
      releaseMaterialProgramReference( material );
    } else {
      // only rebuild uniform list
      programChange = false;
    }

    if ( programChange ) {
      materialProperties.__webglShader = {
        name: material.type,
        uniforms: material.uniforms,
        vertexShader: material.vertexShader,
        fragmentShader: material.fragmentShader,
      };

      material.__webglShader = materialProperties.__webglShader;

      program = programCache.acquireProgram( material, parameters, code );

      materialProperties.program = program;
      material.program = program;
    }

    var uniforms = materialProperties.__webglShader.uniforms;

    materialProperties.fog = fog;

    var progUniforms = materialProperties.program.getUniforms(),
      uniformsList =
      WebGLUniforms.seqWithValue( progUniforms.seq, uniforms );

    materialProperties.uniformsList = uniformsList;
  }

  function setMaterial( material ) {
    material.transparent === true ?
      state.setBlending( NormalBlending, material.premultipliedAlpha )
      : state.setBlending( NoBlending );

    state.setDepthTest( material.depthTest );
    state.setDepthWrite( material.depthWrite );
  }

  function setProgram( camera, fog, material, object ) {
    _usedTextureUnits = 0;

    var materialProperties = properties.get( material );

    if ( material.needsUpdate === false ) {
      if ( materialProperties.program === undefined ) {
        material.needsUpdate = true;
      } else if ( material.fog && materialProperties.fog !== fog ) {
        material.needsUpdate = true;
      }
    }

    if ( material.needsUpdate ) {
      initMaterial( material, fog, object );
      material.needsUpdate = false;
    }

    var refreshProgram = false;
    var refreshMaterial = false;

    var program = materialProperties.program,
      p_uniforms = program.getUniforms(),
      m_uniforms = materialProperties.__webglShader.uniforms;

    if ( program.id !== _currentProgram ) {
      _gl.useProgram( program.program );
      _currentProgram = program.id;

      refreshProgram = true;
      refreshMaterial = true;
    }

    if ( material.id !== _currentMaterialId ) {
      _currentMaterialId = material.id;

      refreshMaterial = true;
    }

    if ( refreshProgram || camera !== _currentCamera ) {
      p_uniforms.set( _gl, camera, 'projectionMatrix' );

      if ( camera !== _currentCamera ) {
        _currentCamera = camera;

        // lighting uniforms depend on the camera so enforce an update
        // now, in case this material supports lights - or later, when
        // the next material that does gets activated:

        refreshMaterial = true;   // set to true on material change
      }

      // load material specific uniforms
      // (shader material also gets them for the sake of genericity)

      if ( material.isShaderMaterial ) {
        p_uniforms.setValue( _gl, 'viewMatrix', camera.matrixWorldInverse );
      }
    }

    if ( refreshMaterial ) {
      // refresh uniforms common to several materials

      if ( fog && material.fog ) {
        refreshUniformsFog( m_uniforms, fog );
      }

      // refresh single material specific uniforms

      WebGLUniforms.upload(
        _gl, materialProperties.uniformsList, m_uniforms, _this );
    }


    // common matrices

    p_uniforms.set( _gl, object, 'modelViewMatrix' );
    p_uniforms.setValue( _gl, 'modelMatrix', object.matrixWorld );

    return program;
  }

  // Uniforms (refresh uniforms objects)

  function refreshUniformsFog( uniforms, fog ) {
    uniforms.fogColor.value = fog.color;

    if ( fog.isFog ) {
      uniforms.fogNear.value = fog.near;
      uniforms.fogFar.value = fog.far;
    }
  }


  // Textures

  function allocTextureUnit() {
    var textureUnit = _usedTextureUnits;
    _usedTextureUnits += 1;
    return textureUnit;
  }

  this.allocTextureUnit = allocTextureUnit;

  // this.setTexture2D = setTexture2D;
  this.setTexture2D = ( function () {
    var warned = false;

    // backwards compatibility: peel texture.texture
    return function setTexture2D( texture, slot ) {
      if ( texture && texture.isWebGLRenderTarget ) {
        if ( ! warned ) {
          console.warn( 'WebGLRenderer.setTexture2D: don\'t use render targets as textures. Use their .texture property instead.' );
          warned = true;
        }

        texture = texture.texture;
      }

      textures.setTexture2D( texture, slot );
    };
  }() );

  this.setRenderTarget = function ( renderTarget ) {
    _currentRenderTarget = renderTarget;

    if ( renderTarget && properties.get( renderTarget ).__webglFramebuffer === undefined ) {
      textures.setupRenderTarget( renderTarget );
    }

    var framebuffer;

    if ( renderTarget ) {
      var renderTargetProperties = properties.get( renderTarget );

      framebuffer = renderTargetProperties.__webglFramebuffer;

      _currentViewport.copy( renderTarget.viewport );
    } else {
      framebuffer = null;

      _currentViewport.copy( _viewport ).multiplyScalar( _pixelRatio );
    }

    if ( _currentFramebuffer !== framebuffer ) {
      _gl.bindFramebuffer( _gl.FRAMEBUFFER, framebuffer );
      _currentFramebuffer = framebuffer;
    }

    state.viewport( _currentViewport );
  };
}

/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
*/

function Fog( color, near, far ) {
  this.name = '';

  this.color = new Color( color );

  this.near = ( near !== undefined ) ? near : 1;
  this.far = ( far !== undefined ) ? far : 1000;
}

Fog.prototype.isFog = true;


/**
* @author mrdoob / http://mrdoob.com/
*/

function Scene() {
  Object3D.call( this );

  this.type = 'Scene';

  this.fog = null;

  this.autoUpdate = true; // checked by the renderer
}

Scene.prototype = Object.create( Object3D.prototype );

Scene.prototype.constructor = Scene;


/**
* @author mrdoob / http://mrdoob.com/
*/

function Line( geometry, material, mode ) {
  Object3D.call( this );

  this.type = 'Line';

  this.geometry = geometry !== undefined ? geometry : new BufferGeometry();
  this.material = material;
}

Line.prototype = Object.assign( Object.create( Object3D.prototype ), {

  constructor: Line,

  isLine: true,

} );

/**
* @author mrdoob / http://mrdoob.com/
*/

function LineSegments( geometry, material ) {
  Line.call( this, geometry, material );

  this.type = 'LineSegments';
}

LineSegments.prototype = Object.assign( Object.create( Line.prototype ), {

  constructor: LineSegments,

  isLineSegments: true,

} );


/**
* @author alteredq / http://alteredqualia.com/
*/

function Points( geometry, material ) {
  Object3D.call( this );

  this.type = 'Points';

  this.geometry = geometry !== undefined ? geometry : new BufferGeometry();
  this.material = material;
}

Points.prototype = Object.assign( Object.create( Object3D.prototype ), {

  constructor: Points,

  isPoints: true,

} );

// kept for compatibility with THREE
function AmbientLight( color ) {}


/**
* @author zz85 / http://www.lab4games.net/zz85/blog
* Extensible curve object
*
**/

/**************************************************************
* Abstract Curve base class
**************************************************************/

function Curve() {}

Curve.prototype = {

  constructor: Curve,

  // Get sequence of points using getPoint( t )

  getPoints: function ( divisions ) {
    var points = [];
    for ( var d = 0; d <= divisions; d ++ ) {
      points.push( this.getPoint( d / divisions ) );
    }
    return points;
  },

};

// TODO: Transformation for Curves?

/**************************************************************
* 3D Curves
**************************************************************/

// A Factory method for creating new curve subclasses

Curve.create = function ( constructor, getPointFunc ) {
  constructor.prototype = Object.create( Curve.prototype );
  constructor.prototype.constructor = constructor;
  constructor.prototype.getPoint = getPointFunc;

  return constructor;
};

/**
* @author zz85 https://github.com/zz85
*
* Centripetal CatmullRom Curve - which is useful for avoiding
* cusps and self-intersections in non-uniform catmull rom curves.
* http://www.cemyuksel.com/research/catmullrom_param/catmullrom.pdf
*
* curve.type accepts centripetal(default), chordal and catmullrom
* curve.tension is used for catmullrom which defaults to 0.5
*/

var CatmullRomCurve3 = ( function () {
var
  tmp = new Vector3(),
  px = new CubicPoly(),
  py = new CubicPoly(),
  pz = new CubicPoly();

/*
Based on an optimized c++ solution in
 - http://stackoverflow.com/questions/9489736/catmull-rom-curve-with-no-cusps-and-no-self-intersections/
 - http://ideone.com/NoEbVM

This CubicPoly class could be used for reusing some variables and calculations,
but for three.js curve use, it could be possible inlined and flatten into a single function call
which can be placed in CurveUtils.
*/

function CubicPoly() {}

/*
 * Compute coefficients for a cubic polynomial
 *   p(s) = c0 + c1*s + c2*s^2 + c3*s^3
 * such that
 *   p(0) = x0, p(1) = x1
 *  and
 *   p'(0) = t0, p'(1) = t1.
 */
CubicPoly.prototype.init = function ( x0, x1, t0, t1 ) {
  this.c0 = x0;
  this.c1 = t0;
  this.c2 = - 3 * x0 + 3 * x1 - 2 * t0 - t1;
  this.c3 = 2 * x0 - 2 * x1 + t0 + t1;
};

CubicPoly.prototype.initNonuniformCatmullRom = function ( x0, x1, x2, x3, dt0, dt1, dt2 ) {
  // compute tangents when parameterized in [t1,t2]
  var t1 = ( x1 - x0 ) / dt0 - ( x2 - x0 ) / ( dt0 + dt1 ) + ( x2 - x1 ) / dt1;
  var t2 = ( x2 - x1 ) / dt1 - ( x3 - x1 ) / ( dt1 + dt2 ) + ( x3 - x2 ) / dt2;

  // rescale tangents for parametrization in [0,1]
  t1 *= dt1;
  t2 *= dt1;

  // initCubicPoly
  this.init( x1, x2, t1, t2 );
};

CubicPoly.prototype.calc = function ( t ) {
  var t2 = t * t;
  var t3 = t2 * t;
  return this.c0 + this.c1 * t + this.c2 * t2 + this.c3 * t3;
};

// Subclass Three.js curve
return Curve.create(

  function ( p /* array of Vector3 */ ) {
    this.points = p || [];
    this.closed = false;
  },

  function ( t ) {
    var points = this.points,
      point, intPoint, weight, l;

    l = points.length;

    if ( l < 2 ) { console.log( 'duh, you need at least 2 points' ); }

    point = ( l - ( this.closed ? 0 : 1 ) ) * t;
    intPoint = Math.floor( point );
    weight = point - intPoint;

    if ( this.closed ) {
      intPoint += intPoint > 0 ? 0 : ( Math.floor( Math.abs( intPoint ) / points.length ) + 1 ) * points.length;
    } else if ( weight === 0 && intPoint === l - 1 ) {
      intPoint = l - 2;
      weight = 1;
    }

    var p0, p1, p2, p3; // 4 points

    if ( this.closed || intPoint > 0 ) {
      p0 = points[( intPoint - 1 ) % l];
    } else {
      // extrapolate first point
      tmp.subVectors( points[0], points[1] ).add( points[0] );
      p0 = tmp;
    }

    p1 = points[intPoint % l];
    p2 = points[( intPoint + 1 ) % l];

    if ( this.closed || intPoint + 2 < l ) {
      p3 = points[( intPoint + 2 ) % l];
    } else {
      // extrapolate last point
      tmp.subVectors( points[l - 1], points[l - 2] ).add( points[l - 1] );
      p3 = tmp;
    }

    // init Centripetal Catmull-Rom
    var pow = 0.25;
    var dt0 = Math.pow( p0.distanceToSquared( p1 ), pow );
    var dt1 = Math.pow( p1.distanceToSquared( p2 ), pow );
    var dt2 = Math.pow( p2.distanceToSquared( p3 ), pow );

    // safety check for repeated points
    if ( dt1 < 1e-4 ) { dt1 = 1.0; }
    if ( dt0 < 1e-4 ) { dt0 = dt1; }
    if ( dt2 < 1e-4 ) { dt2 = dt1; }

    px.initNonuniformCatmullRom( p0.x, p1.x, p2.x, p3.x, dt0, dt1, dt2 );
    py.initNonuniformCatmullRom( p0.y, p1.y, p2.y, p3.y, dt0, dt1, dt2 );
    pz.initNonuniformCatmullRom( p0.z, p1.z, p2.z, p3.z, dt0, dt1, dt2 );

    var v = new Vector3(
      px.calc( weight ),
      py.calc( weight ),
      pz.calc( weight )
    );

    return v;
  }

);
} )();

// @flow

/*:: type Num3 = [number, number, number] */
/*:: import type {AtomT} from './model.js' */

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

function makeColorAttribute(colors /*:Color[]*/) {
  var col = new Float32Array(colors.length * 3);
  for (var i = 0; i < colors.length; i++) {
    col[3*i+0] = colors[i].r;
    col[3*i+1] = colors[i].g;
    col[3*i+2] = colors[i].b;
  }
  return new BufferAttribute(col, 3);
}

var light_dir = new Vector3(-0.2, 0.3, 1.0); // length affects brightness

var fog_pars_fragment =
"#ifdef USE_FOG\nuniform vec3 fogColor;\nuniform float fogNear;\nuniform float fogFar;\n#endif";

var fog_end_fragment =
"#ifdef USE_FOG\n  float depth = gl_FragCoord.z / gl_FragCoord.w;\n  float fogFactor = smoothstep(fogNear, fogFar, depth);\n  gl_FragColor.rgb = mix(gl_FragColor.rgb, fogColor, fogFactor);\n#endif";


var varcolor_vert = "\nattribute vec3 color;\nvarying vec3 vcolor;\nvoid main() {\n  vcolor = color;\n  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);\n}\n";

var unicolor_vert = "\nvoid main() {\n  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);\n}\n";

var unicolor_frag = "\n" + fog_pars_fragment + "\nuniform vec3 vcolor;\nvoid main() {\n  gl_FragColor = vec4(vcolor, 1.0);\n" + fog_end_fragment + "\n}";

var varcolor_frag = "\n" + fog_pars_fragment + "\nvarying vec3 vcolor;\nvoid main() {\n  gl_FragColor = vec4(vcolor, 1.0);\n" + fog_end_fragment + "\n}";

function makeLines(pos /*:Float32Array*/, color /*:Color*/,
                          linewidth /*:number*/) {
  var material = new ShaderMaterial({
    uniforms: makeUniforms({vcolor: color}),
    vertexShader: unicolor_vert,
    fragmentShader: unicolor_frag,
    fog: true,
    linewidth: linewidth,
    type: 'um_lines',
  });
  var geometry = new BufferGeometry();
  geometry.addAttribute('position', new BufferAttribute(pos, 3));
  return new LineSegments(geometry, material);
}

function makeCube(size /*:number*/,
                         ctr /*:Vector3*/,
                         options /*:{[key:string]: any}*/) {
  var pos = new Float32Array(CUBE_EDGES.length * 3);
  for (var i = 0; i < CUBE_EDGES.length; i++) {
    var coor = CUBE_EDGES[i];
    pos[3*i+0] = ctr.x + size * (coor[0] - 0.5);
    pos[3*i+1] = ctr.y + size * (coor[1] - 0.5);
    pos[3*i+2] = ctr.z + size * (coor[2] - 0.5);
  }
  return makeLines(pos, options.color, options.linewidth);
}

function makeMultiColorLines(pos /*:Float32Array*/,
                                    colors /*:Color[]*/,
                                    linewidth /*:number*/) {
  var material = new ShaderMaterial({
    uniforms: makeUniforms({}),
    vertexShader: varcolor_vert,
    fragmentShader: varcolor_frag,
    fog: true,
    linewidth: linewidth,
    type: 'um_multicolor_lines',
  });
  var geometry = new BufferGeometry();
  geometry.addAttribute('position', new BufferAttribute(pos, 3));
  geometry.addAttribute('color', makeColorAttribute(colors));
  return new LineSegments(geometry, material);
}

// A cube with 3 edges (for x, y, z axes) colored in red, green and blue.
function makeRgbBox(transform_func /*:Num3 => Num3*/, color /*:Color*/) {
  var pos = new Float32Array(CUBE_EDGES.length * 3);
  for (var i = 0; i < CUBE_EDGES.length; i++) {
    var coor = transform_func(CUBE_EDGES[i]);
    pos[3*i+0] = coor[0];
    pos[3*i+1] = coor[1];
    pos[3*i+2] = coor[2];
  }
  var colors = [
    new Color(0xff0000), new Color(0xffaa00),
    new Color(0x00ff00), new Color(0xaaff00),
    new Color(0x0000ff), new Color(0x00aaff) ];
  for (var j = 6; j < CUBE_EDGES.length; j++) {
    colors.push(color);
  }
  return makeMultiColorLines(pos, colors, 1);
}

function double_pos(pos /*:Num3[]*/) {
  var double_pos = [];
  for (var i = 0; i < pos.length; i++) {
    var v = pos[i];
    double_pos.push(v[0], v[1], v[2]);
    double_pos.push(v[0], v[1], v[2]);
  }
  return double_pos;
}

function double_color(color_arr /*:Color[]*/) {
  var len = color_arr.length;
  var color = new Float32Array(6*len);
  for (var i = 0; i < len; i++) {
    var col = color_arr[i];
    color[6*i] = col.r;
    color[6*i+1] = col.g;
    color[6*i+2] = col.b;
    color[6*i+3] = col.r;
    color[6*i+4] = col.g;
    color[6*i+5] = col.b;
  }
  return color;
}

// draw quads as 2 triangles: 4 attributes / quad, 6 indices / quad
function make_quad_index_buffer(len) {
  var index = (4*len < 65536 ? new Uint16Array(6*len)
                             : new Uint32Array(6*len));
  var vert_order = [0, 1, 2, 0, 2, 3];
  for (var i = 0; i < len; i++) {
    for (var j = 0; j < 6; j++) {
      index[6*i+j] = 4*i + vert_order[j];
    }
  }
  return new BufferAttribute(index, 1);
}


var wide_line_vert = [
  'attribute vec3 color;',
  'attribute vec3 previous;',
  'attribute vec3 next;',
  'attribute float side;',
  'uniform vec2 win_size;',
  'uniform float linewidth;',
  'varying vec3 vcolor;',

  'void main() {',
  '  vcolor = color;',
  '  mat4 mat = projectionMatrix * modelViewMatrix;',
  '  vec2 dir1 = (mat * vec4(next - position, 0.0)).xy * win_size;',
  '  float len = length(dir1);',
  '  if (len > 0.0) dir1 /= len;',
  '  vec2 dir2 = (mat * vec4(position - previous, 0.0)).xy * win_size;',
  '  len = length(dir2);',
  '  dir2 = len > 0.0 ? dir2 / len : dir1;',
  '  vec2 tang = normalize(dir1 + dir2);',
  '  vec2 normal = vec2(-tang.y, tang.x);',
  // Now we have more or less a miter join. Bavel join could be more
  // appropriate, but it'd require one more triangle and more complex shader.
  // max() is a trade-off between too-long miters and too-thin lines.
  // The outer vertex should not go too far, the inner is not a problem.
  '  float outer = side * dot(dir2, normal);',
  '  float angle_factor = max(dot(tang, dir2), outer > 0.0 ? 0.5 : 0.1);',
  '  gl_Position = mat * vec4(position, 1.0);',
  '  gl_Position.xy += side * linewidth / angle_factor * normal / win_size;',
  '}'].join('\n');

var wide_segments_vert = "\nattribute vec3 color;\nattribute vec3 other;\nattribute float side;\nuniform vec2 win_size;\nuniform float linewidth;\nvarying vec3 vcolor;\n\nvoid main() {\n  vcolor = color;\n  mat4 mat = projectionMatrix * modelViewMatrix;\n  vec2 dir = normalize((mat * vec4(position - other, 0.0)).xy);\n  vec2 normal = vec2(-dir.y, dir.x);\n  gl_Position = mat * vec4(position, 1.0);\n  gl_Position.xy += side * linewidth * normal / win_size;\n}";

function interpolate_vertices(segment, smooth) /*:Vector3[]*/{
  var vertices = [];
  for (var i = 0; i < segment.length; i++) {
    var xyz = segment[i].xyz;
    vertices.push(new Vector3(xyz[0], xyz[1], xyz[2]));
  }
  if (!smooth || smooth < 2) { return vertices; }
  var curve = new CatmullRomCurve3(vertices);
  return curve.getPoints((segment.length - 1) * smooth);
}

function interpolate_colors(colors, smooth) {
  if (!smooth || smooth < 2) { return colors; }
  var ret = [];
  for (var i = 0; i < colors.length - 1; i++) {
    for (var j = 0; j < smooth; j++) {
      // currently we don't really interpolate colors
      ret.push(colors[i]);
    }
  }
  ret.push(colors[colors.length - 1]);
  return ret;
}

// a simplistic linear interpolation, no need to SLERP
function interpolate_directions(dirs, smooth) {
  smooth = smooth || 1;
  var ret = [];
  var i;
  for (i = 0; i < dirs.length - 1; i++) {
    var p = dirs[i];
    var n = dirs[i+1];
    for (var j = 0; j < smooth; j++) {
      var an = j / smooth;
      var ap = 1 - an;
      ret.push(ap*p[0] + an*n[0], ap*p[1] + an*n[1], ap*p[2] + an*n[2]);
    }
  }
  ret.push(dirs[i][0], dirs[i][1], dirs[i][2]);
  return ret;
}

function makeUniforms(params/*:{[id:string]:mixed}*/) {
  var uniforms = {
    fogNear: { value: null },  // will be updated in setProgram()
    fogFar: { value: null },
    fogColor: { value: null },
  };
  for (var p in params) {  // eslint-disable-line guard-for-in
    uniforms[p] = { value: params[p] };
  }
  return uniforms;
}

var ribbon_vert = "\nattribute vec3 color;\nattribute vec3 tan;\nuniform float shift;\nvarying vec3 vcolor;\nvoid main() {\n  vcolor = color;\n  vec3 pos = position + shift * normalize(tan);\n  gl_Position = projectionMatrix * modelViewMatrix * vec4(pos, 1.0);\n}";

// 9-line ribbon
function makeRibbon(vertices /*:AtomT[]*/,
                           colors /*:Color[]*/,
                           tangents /*:Num3[]*/,
                           smoothness /*:number*/) {
  var vertex_arr = interpolate_vertices(vertices, smoothness);
  var color_arr = interpolate_colors(colors, smoothness);
  var tang_arr = interpolate_directions(tangents, smoothness);
  var obj = new Object3D();
  var geometry = new BufferGeometry();
  var pos = new Float32Array(vertex_arr.length * 3);
  for (var i = 0; i < vertex_arr.length; i++) {
    var v = vertex_arr[i];
    pos[3*i+0] = v.x;
    pos[3*i+1] = v.y;
    pos[3*i+2] = v.z;
  }
  geometry.addAttribute('position', new BufferAttribute(pos, 3));
  geometry.addAttribute('color', makeColorAttribute(color_arr));
  var tan = new Float32Array(tang_arr);
  geometry.addAttribute('tan', new BufferAttribute(tan, 3));
  for (var n = -4; n < 5; n++) {
    var material = new ShaderMaterial({
      uniforms: makeUniforms({shift: 0.1 * n}),
      vertexShader: ribbon_vert,
      fragmentShader: varcolor_frag,
      fog: true,
      type: 'um_ribbon',
    });
    obj.add(new Line(geometry, material));
  }
  return obj;
}


function makeChickenWire(data /*:{vertices: number[], segments: number[]}*/,
                         options /*:{[key: string]: any}*/) {
  var geom = new BufferGeometry();
  var position = new Float32Array(data.vertices);
  geom.addAttribute('position', new BufferAttribute(position, 3));

  // Although almost all browsers support OES_element_index_uint nowadays,
  // use Uint32 indexes only when needed.
  var arr = (data.vertices.length < 3*65536 ? new Uint16Array(data.segments)
                                              : new Uint32Array(data.segments));
  //console.log('arr len:', data.vertices.length, data.segments.length);
  geom.setIndex(new BufferAttribute(arr, 1));
  var material = new ShaderMaterial({
    uniforms: makeUniforms({vcolor: options.color}),
    vertexShader: unicolor_vert,
    fragmentShader: unicolor_frag,
    fog: true,
    linewidth: options.linewidth,
    type: 'um_line_chickenwire',
  });
  return new LineSegments(geom, material);
}


var grid_vert = "\nuniform vec3 ucolor;\nuniform vec3 fogColor;\nvarying vec4 vcolor;\nvoid main() {\n  vec2 scale = vec2(projectionMatrix[0][0], projectionMatrix[1][1]);\n  float z = position.z;\n  float fogFactor = (z > 0.5 ? 0.2 : 0.7);\n  float alpha = 0.8 * smoothstep(z > 1.5 ? -10.0 : 0.01, 0.1, scale.y);\n  vcolor = vec4(mix(ucolor, fogColor, fogFactor), alpha);\n  gl_Position = vec4(position.xy * scale, -0.99, 1.0);\n}";

var grid_frag = "\nvarying vec4 vcolor;\nvoid main() {\n  gl_FragColor = vcolor;\n}";

function makeGrid() {
  var N = 50;
  var pos = [];
  for (var i = -N; i <= N; i++) {
    var z = 0; // z only marks major/minor axes
    if (i % 5 === 0) { z = i % 2 === 0 ? 2 : 1; }
    pos.push(-N, i, z, N, i, z);  // horizontal line
    pos.push(i, -N, z, i, N, z);  // vertical line
  }
  var geom = new BufferGeometry();
  var pos_arr = new Float32Array(pos);
  geom.addAttribute('position', new BufferAttribute(pos_arr, 3));
  var material = new ShaderMaterial({
    uniforms: makeUniforms({ucolor: new Color(0x888888)}),
    //linewidth: 3,
    vertexShader: grid_vert,
    fragmentShader: grid_frag,
    fog: true, // no really, but we use fogColor
    type: 'um_grid',
  });
  material.transparent = true;
  var obj = new LineSegments(geom, material);
  obj.frustumCulled = false;  // otherwise the renderer could skip it
  obj.color_value = material.uniforms.ucolor.value; // shortcut
  return obj;
}


function makeLineMaterial(options /*:{[key: string]: mixed}*/) {
  var uniforms = makeUniforms({
    linewidth: options.linewidth,
    win_size: options.win_size,
  });
  return new ShaderMaterial({
    uniforms: uniforms,
    vertexShader: options.segments ? wide_segments_vert : wide_line_vert,
    fragmentShader: varcolor_frag,
    fog: true,
    type: 'um_line',
  });
}

// vertex_arr and color_arr must be of the same length
function makeLine(material /*:ShaderMaterial*/,
                         vertex_arr /*:Num3[]*/,
                         color_arr /*:Color[]*/) {
  var len = vertex_arr.length;
  var pos = double_pos(vertex_arr);
  var position = new Float32Array(pos);
  // could we use three overlapping views of the same buffer?
  var previous = new Float32Array(6*len);
  var i;
  for (i = 0; i < 6; i++) { previous[i] = pos[i]; }
  for (; i < 6 * len; i++) { previous[i] = pos[i-6]; }
  var next = new Float32Array(6*len);
  for (i = 0; i < 6 * (len-1); i++) { next[i] = pos[i+6]; }
  for (; i < 6 * len; i++) { next[i] = pos[i]; }
  var side = new Float32Array(2*len);
  for (i = 0; i < len; i++) {
    side[2*i] = 1;
    side[2*i+1] = -1;
  }
  var color = double_color(color_arr);
  var geometry = new BufferGeometry();
  geometry.addAttribute('position', new BufferAttribute(position, 3));
  geometry.addAttribute('previous', new BufferAttribute(previous, 3));
  geometry.addAttribute('next', new BufferAttribute(next, 3));
  geometry.addAttribute('side', new BufferAttribute(side, 1));
  geometry.addAttribute('color', new BufferAttribute(color, 3));

  var mesh = new Mesh(geometry, material);
  mesh.drawMode = TriangleStripDrawMode;
  mesh.userData.bond_lines = true;
  return mesh;
}

// vertex_arr and color_arr must be of the same length
function makeLineSegments(material /*:ShaderMaterial*/,
                                 vertex_arr /*:Num3[]*/,
                                 color_arr /*:?Color[]*/) {
  // n input vertices => 2n output vertices, n triangles, 3n indexes
  var len = vertex_arr.length;
  var pos = double_pos(vertex_arr);
  var position = new Float32Array(pos);
  var other_vert = new Float32Array(6*len);
  for (var i = 0; i < 6 * len; i += 12) {
    var j = 0;
    for (; j < 6; j++) { other_vert[i+j] = pos[i+j+6]; }
    for (; j < 12; j++) { other_vert[i+j] = pos[i+j-6]; }
  }
  var side = new Float32Array(2*len);
  for (var k = 0; k < len; k++) {
    side[2*k] = -1;
    side[2*k+1] = 1;
  }
  var geometry = new BufferGeometry();
  geometry.addAttribute('position', new BufferAttribute(position, 3));
  geometry.addAttribute('other', new BufferAttribute(other_vert, 3));
  geometry.addAttribute('side', new BufferAttribute(side, 1));
  if (color_arr != null) {
    var color = double_color(color_arr);
    geometry.addAttribute('color', new BufferAttribute(color, 3));
  }
  geometry.setIndex(make_quad_index_buffer(len/2));

  var mesh = new Mesh(geometry, material);
  mesh.userData.bond_lines = true;
  return mesh;
}

var wheel_vert = "\nattribute vec3 color;\nuniform float size;\nvarying vec3 vcolor;\nvoid main() {\n  vcolor = color;\n  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);\n  gl_PointSize = size;\n}";

// not sure how portable it is
var wheel_frag = "\n" + fog_pars_fragment + "\nvarying vec3 vcolor;\nvoid main() {\n  vec2 diff = gl_PointCoord - vec2(0.5, 0.5);\n  if (dot(diff, diff) >= 0.25) discard;\n  gl_FragColor = vec4(vcolor, 1.0);\n" + fog_end_fragment + "\n}";

function makeWheels(atom_arr /*:AtomT[]*/,
                           color_arr /*:Color[]*/,
                           size /*:number*/) {
  var geometry = new BufferGeometry();
  var pos = new Float32Array(atom_arr.length * 3);
  for (var i = 0; i < atom_arr.length; i++) {
    var xyz = atom_arr[i].xyz;
    pos[3*i+0] = xyz[0];
    pos[3*i+1] = xyz[1];
    pos[3*i+2] = xyz[2];
  }
  geometry.addAttribute('position', new BufferAttribute(pos, 3));
  geometry.addAttribute('color', makeColorAttribute(color_arr));
  var material = new ShaderMaterial({
    uniforms: makeUniforms({size: size}),
    vertexShader: wheel_vert,
    fragmentShader: wheel_frag,
    fog: true,
    type: 'um_wheel',
  });
  var obj = new Points(geometry, material);
  return obj;
}

// For the ball-and-stick rendering we use so-called imposters.
// This technique was described in:
// http://doi.ieeecomputersociety.org/10.1109/TVCG.2006.115
// free copy here:
// http://vcg.isti.cnr.it/Publications/2006/TCM06/Tarini_FinalVersionElec.pdf
// and was nicely summarized in:
// http://www.sunsetlakesoftware.com/2011/05/08/enhancing-molecules-using-opengl-es-20

var sphere_vert = "\nattribute vec3 color;\nattribute vec2 corner;\nuniform float radius;\nvarying vec3 vcolor;\nvarying vec2 vcorner;\nvarying vec3 vpos;\n\nvoid main() {\n  vcolor = color;\n  vcorner = corner;\n  vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);\n  vpos = mvPosition.xyz;\n  mvPosition.xy += corner * radius;\n  gl_Position = projectionMatrix * mvPosition;\n}\n";

// based on 3Dmol imposter shaders
var sphere_frag = "\n" + fog_pars_fragment + "\nuniform mat4 projectionMatrix;\nuniform vec3 lightDir;\nuniform float radius;\nvarying vec3 vcolor;\nvarying vec2 vcorner;\nvarying vec3 vpos;\n\nvoid main() {\n  float sq = dot(vcorner, vcorner);\n  if (sq > 1.0) discard;\n  float z = sqrt(1.0-sq);\n  vec3 xyz = vec3(vcorner.x, vcorner.y, z);\n  vec4 projPos = projectionMatrix * vec4(vpos + radius * xyz, 1.0);\n  gl_FragDepthEXT = 0.5 * ((gl_DepthRange.diff * (projPos.z / projPos.w)) +\n                           gl_DepthRange.near + gl_DepthRange.far);\n  float weight = clamp(dot(xyz, lightDir), 0.0, 1.0) * 0.8 + 0.2;\n  gl_FragColor = vec4(weight * vcolor, 1.0);\n  " + fog_end_fragment + "\n}\n";

var stick_vert = "\nattribute vec3 color;\nattribute vec3 axis;\nattribute vec2 corner;\nuniform float radius;\nvarying vec3 vcolor;\nvarying vec2 vcorner;\nvarying vec3 vpos;\nvarying vec3 vaxis;\n\nvoid main() {\n  vcolor = color;\n  vcorner = corner;\n  vaxis = normalize((modelViewMatrix * vec4(axis, 0.0)).xyz);\n  vec2 normal = normalize(vec2(-vaxis.y, vaxis.x));\n  vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);\n  vpos = mvPosition.xyz;\n  mvPosition.xy += corner[1] * radius * normal;\n  gl_Position = projectionMatrix * mvPosition;\n}";

var stick_frag = "\n" + fog_pars_fragment + "\nuniform mat4 projectionMatrix;\nuniform vec3 lightDir;\nuniform float radius;\nvarying vec3 vcolor;\nvarying vec2 vcorner;\nvarying vec3 vpos;\nvarying vec3 vaxis;\nvoid main() {\n  float central = 1.0 - vcorner[1] * vcorner[1];\n  vec4 pos = vec4(vpos, 1.0);\n  pos.z += radius * vaxis.z * central;\n  vec4 projPos = projectionMatrix * pos;\n  gl_FragDepthEXT = 0.5 * ((gl_DepthRange.diff * (projPos.z / projPos.w)) +\n                           gl_DepthRange.near + gl_DepthRange.far);\n  float weight = length(cross(vaxis, lightDir)) * central * 0.8 + 0.2;\n  gl_FragColor = vec4(min(weight, 1.0) * vcolor, 1.0);\n" + fog_end_fragment + "\n}";

function makeSticks(vertex_arr /*:Num3[]*/,
                           color_arr /*:Color[]*/,
                           radius /*:number*/) {
  var uniforms = makeUniforms({
    radius: radius,
    lightDir: light_dir,
  });
  var material = new ShaderMaterial({
    uniforms: uniforms,
    vertexShader: stick_vert,
    fragmentShader: stick_frag,
    fog: true,
    type: 'um_stick',
  });
  material.extensions.fragDepth = true;

  var len = vertex_arr.length;
  var pos = double_pos(vertex_arr);
  var position = new Float32Array(pos);
  var axis = new Float32Array(6*len);
  for (var i = 0; i < 6 * len; i += 12) {
    for (var j = 0; j < 6; j++) { axis[i+j] = pos[i+j+6] - pos[i+j]; }
    for (var j$1 = 0; j$1 < 6; j$1++) { axis[i+j$1+6] = axis[i+j$1]; }
  }
  var geometry = new BufferGeometry();
  geometry.addAttribute('position', new BufferAttribute(position, 3));
  var corner = new Float32Array(4*len);
  for (var i$1 = 0; 2 * i$1 < len; i$1++) {
    corner[8*i$1 + 0] = -1;  // 0
    corner[8*i$1 + 1] = -1;  // 0
    corner[8*i$1 + 2] = -1;  // 1
    corner[8*i$1 + 3] = +1;  // 1
    corner[8*i$1 + 4] = +1;  // 2
    corner[8*i$1 + 5] = +1;  // 2
    corner[8*i$1 + 6] = +1;  // 3
    corner[8*i$1 + 7] = -1;  // 3
  }
  geometry.addAttribute('axis', new BufferAttribute(axis, 3));
  geometry.addAttribute('corner', new BufferAttribute(corner, 2));
  var color = double_color(color_arr);
  geometry.addAttribute('color', new BufferAttribute(color, 3));
  geometry.setIndex(make_quad_index_buffer(len/2));

  var mesh = new Mesh(geometry, material);
  mesh.userData.bond_lines = true;
  return mesh;
}

function makeBalls(atom_arr /*:AtomT[]*/,
                          color_arr /*:Color[]*/,
                          radius /*:number*/) {
  var N = atom_arr.length;
  var geometry = new BufferGeometry();

  var pos = new Float32Array(N * 4 * 3);
  for (var i = 0; i < N; i++) {
    var xyz = atom_arr[i].xyz;
    for (var j = 0; j < 4; j++) {
      for (var k = 0; k < 3; k++) {
        pos[3 * (4*i + j) + k] = xyz[k];
      }
    }
  }
  geometry.addAttribute('position', new BufferAttribute(pos, 3));

  var corner = new Float32Array(N * 4 * 2);
  for (var i$1 = 0; i$1 < N; i$1++) {
    corner[8*i$1 + 0] = -1;  // 0
    corner[8*i$1 + 1] = -1;  // 0
    corner[8*i$1 + 2] = -1;  // 1
    corner[8*i$1 + 3] = +1;  // 1
    corner[8*i$1 + 4] = +1;  // 2
    corner[8*i$1 + 5] = +1;  // 2
    corner[8*i$1 + 6] = +1;  // 3
    corner[8*i$1 + 7] = -1;  // 3
  }
  geometry.addAttribute('corner', new BufferAttribute(corner, 2));

  var colors = new Float32Array(N * 4 * 3);
  for (var i$2 = 0; i$2 < N; i$2++) {
    var col = color_arr[i$2];
    for (var j$1 = 0; j$1 < 4; j$1++) {
      colors[3 * (4*i$2 + j$1) + 0] = col.r;
      colors[3 * (4*i$2 + j$1) + 1] = col.g;
      colors[3 * (4*i$2 + j$1) + 2] = col.b;
    }
  }
  geometry.addAttribute('color', new BufferAttribute(colors, 3));

  geometry.setIndex(make_quad_index_buffer(N));

  var material = new ShaderMaterial({
    uniforms: makeUniforms({
      radius: radius,
      lightDir: light_dir,
    }),
    vertexShader: sphere_vert,
    fragmentShader: sphere_frag,
    fog: true,
    type: 'um_sphere',
  });
  material.extensions.fragDepth = true;
  var obj = new Mesh(geometry, material);
  return obj;
}

// based on Line.prototype.raycast(), but skipping duplicated points
var inverseMatrix = new Matrix4();
var ray = new Ray();
function line_raycast(mesh/*:Mesh*/, options/*:Object*/,
                             intersects/*:Object[]*/) {
  var precisionSq = options.precision * options.precision;
  inverseMatrix.getInverse(mesh.matrixWorld);
  ray.copy(options.ray).applyMatrix4(inverseMatrix);
  var vStart = new Vector3();
  var vEnd = new Vector3();
  var interSegment = new Vector3();
  var interRay = new Vector3();
  var step = mesh.drawMode === TriangleStripDrawMode ? 1 : 2;
  var positions = mesh.geometry.attributes.position.array;
  for (var i = 0, l = positions.length / 6 - 1; i < l; i += step) {
    vStart.fromArray(positions, 6 * i);
    vEnd.fromArray(positions, 6 * i + 6);
    var distSq = ray.distanceSqToSegment(vStart, vEnd, interRay, interSegment);
    if (distSq > precisionSq) { continue; }
    interRay.applyMatrix4(mesh.matrixWorld);
    var distance = options.ray.origin.distanceTo(interRay);
    if (distance < options.near || distance > options.far) { continue; }
    intersects.push({
      distance: distance,
      point: interSegment.clone().applyMatrix4(mesh.matrixWorld),
      index: i,
      object: mesh,
      line_dist: Math.sqrt(distSq), // extra property, not in Three.js
    });
  }
}

function makeCanvasWithText(text, options) {
  if (typeof document === 'undefined') { return; }  // for testing on node
  var canvas = document.createElement('canvas');
  // Canvas size should be 2^N.
  canvas.width = 256;  // arbitrary limit, to keep it simple
  canvas.height = 16;  // font size
  var context = canvas.getContext('2d');
  if (!context) { return null; }
  context.font = (options.font || 'bold 14px') + ' sans-serif';
  //context.fillStyle = 'green';
  //context.fillRect(0, 0, canvas.width, canvas.height);
  context.textBaseline = 'bottom';
  if (options.color) { context.fillStyle = options.color; }
  context.fillText(text, 0, canvas.height);
  return canvas;
}

var label_vert = "\nattribute vec2 uvs;\nuniform vec2 canvas_size;\nuniform vec2 win_size;\nuniform float z_shift;\nvarying vec2 vUv;\nvoid main() {\n  vUv = uvs;\n  vec2 rel_offset = vec2(0.02, -0.3);\n  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);\n  gl_Position.xy += (uvs + rel_offset) * 2.0 * canvas_size / win_size;\n  gl_Position.z += z_shift * projectionMatrix[2][2];\n}";

var label_frag = "\n" + fog_pars_fragment + "\nvarying vec2 vUv;\nuniform sampler2D map;\nvoid main() {\n  gl_FragColor = texture2D(map, vUv);\n" + fog_end_fragment + "\n}";


function makeLabel(text /*:string*/, options /*:{[key:string]: any}*/) {
  var canvas = makeCanvasWithText(text, options);
  if (!canvas) { return; }
  var texture = new Texture(canvas);
  texture.needsUpdate = true;

  // Rectangle geometry.
  var geometry = new BufferGeometry();
  var pos = options.pos;
  var position = new Float32Array([].concat(pos, pos, pos, pos));
  var uvs = new Float32Array([0, 1, 1, 1, 0, 0, 1, 0]);
  var indices = new Uint16Array([0, 2, 1, 2, 3, 1]);
  geometry.setIndex(new BufferAttribute(indices, 1));
  geometry.addAttribute('position', new BufferAttribute(position, 3));
  geometry.addAttribute('uvs', new BufferAttribute(uvs, 2));

  var material = new ShaderMaterial({
    uniforms: makeUniforms({map: texture,
                            canvas_size: [canvas.width, canvas.height],
                            win_size: options.win_size,
                            z_shift: options.z_shift}),
    vertexShader: label_vert,
    fragmentShader: label_frag,
    fog: true,
    type: 'um_label',
  });
  material.transparent = true;
  var mesh = new Mesh(geometry, material);
  mesh.remake = function (text, options) {
    texture.image = makeCanvasWithText(text, options);
    texture.needsUpdate = true;
  };
  return mesh;
}

// Add vertices of a 3d cross (representation of an unbonded atom)
function addXyzCross(vertices /*:Num3[]*/, xyz /*:Num3*/, r /*:number*/) {
  vertices.push([xyz[0]-r, xyz[1], xyz[2]], [xyz[0]+r, xyz[1], xyz[2]]);
  vertices.push([xyz[0], xyz[1]-r, xyz[2]], [xyz[0], xyz[1]+r, xyz[2]]);
  vertices.push([xyz[0], xyz[1], xyz[2]-r], [xyz[0], xyz[1], xyz[2]+r]);
}

// @flow

/*:: import type {OrthographicCamera} from './fromthree.js' */

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

var STATE = { NONE: -1, ROTATE: 0, PAN: 1, ZOOM: 2, PAN_ZOOM: 3,
                       SLAB: 4, ROLL: 5, AUTO_ROTATE: 6, GO: 7 };

var auto_speed = 1.0;

// based on three.js/examples/js/controls/OrthographicTrackballControls.js
var Controls = function Controls(camera /*:OrthographicCamera*/, target /*:Vector3*/) {
  this._camera = camera;
  this._target = target;
  this._state = STATE.NONE;
  this._rotate_start = new Vector3();
  this._rotate_end = new Vector3();
  this._zoom_start = [0, 0];
  this._zoom_end = [0, 0];
  this._pinch_start = 0;
  this._pinch_end = 0;
  this._pan_start = [0, 0];
  this._pan_end = [0, 0];
  this._panned = true;
  this._rotating = 0.0;
  this._auto_stamp = null;
  this._go_func = null;

  // the far plane is more distant from the target than the near plane (3:1)
  this.slab_width = [2.5, 7.5, null];
};

Controls.prototype._rotate_camera = function _rotate_camera (eye /*:Vector3*/) {
  var quat = new Quaternion();
  quat.setFromUnitVectors(this._rotate_end, this._rotate_start);
  eye.applyQuaternion(quat);
  this._camera.up.applyQuaternion(quat);
  this._rotate_end.applyQuaternion(quat);
  this._rotate_start.copy(this._rotate_end);
};

Controls.prototype._zoom_camera = function _zoom_camera (eye /*:Vector3*/) {
  var dx = this._zoom_end[0] - this._zoom_start[0];
  var dy = this._zoom_end[1] - this._zoom_start[1];
  if (this._state === STATE.ZOOM) {
    this._camera.zoom /= (1 - dx + dy);
  } else if (this._state === STATE.SLAB) {
    this._target.addScaledVector(eye, -5.0 / eye.length() * dy);
  } else if (this._state === STATE.ROLL) {
    var quat = new Quaternion();
    quat.setFromAxisAngle(eye, 0.05 * (dx - dy));
    this._camera.up.applyQuaternion(quat);
  }
  this._zoom_start[0] = this._zoom_end[0];
  this._zoom_start[1] = this._zoom_end[1];
  return this._state === STATE.SLAB ? 10*dx : null;
};

Controls.prototype._pan_camera = function _pan_camera (eye /*:Vector3*/) {
  var dx = this._pan_end[0] - this._pan_start[0];
  var dy = this._pan_end[1] - this._pan_start[1];
  dx *= 0.5 * (this._camera.right - this._camera.left) / this._camera.zoom;
  dy *= 0.5 * (this._camera.bottom - this._camera.top) / this._camera.zoom;
  var pan = eye.clone().cross(this._camera.up).setLength(dx);
  pan.addScaledVector(this._camera.up, dy / this._camera.up.length());
  this._camera.position.add(pan);
  this._target.add(pan);
  this._pan_start[0] = this._pan_end[0];
  this._pan_start[1] = this._pan_end[1];
};

Controls.prototype._auto_rotate = function _auto_rotate (eye /*:Vector3*/) {
  this._rotate_start.copy(eye).normalize();
  var now = Date.now();
  var elapsed = (this._auto_stamp !== null ? now - this._auto_stamp : 16.7);
  var speed = 1.8e-5 * elapsed * auto_speed;
  this._auto_stamp = now;
  if (this._rotating === true) {
    speed = -speed;
  } else if (this._rotating !== false) {
    this._rotating += 0.02;
    speed = 4e-5 * auto_speed * Math.cos(this._rotating);
  }
  this._rotate_end.crossVectors(this._camera.up, eye).multiplyScalar(speed)
    .add(this._rotate_start);
};

Controls.prototype.toggle_auto = function toggle_auto (param /*:number|boolean*/) {
  if (this._state === STATE.AUTO_ROTATE &&
      typeof param === typeof this._rotating) {
    this._state = STATE.NONE;
  } else {
    this._state = STATE.AUTO_ROTATE;
    this._auto_stamp = null;
    this._rotating = param;
  }
};

Controls.prototype.is_going = function is_going () { return this._state === STATE.GO; };

Controls.prototype.is_moving = function is_moving () { return this._state !== STATE.NONE; };

Controls.prototype.update = function update () {
  var changed = false;
  var eye = this._camera.position.clone().sub(this._target);
  if (this._state === STATE.AUTO_ROTATE) {
    this._auto_rotate(eye);
  }
  if (!this._rotate_start.equals(this._rotate_end)) {
    this._rotate_camera(eye);
    changed = true;
  }
  if (this._pinch_end !== this._pinch_start) {
    this._camera.zoom *= this._pinch_end / this._pinch_start;
    this._pinch_start = this._pinch_end;
    changed = true;
  }
  if (this._zoom_end[0] !== this._zoom_start[0] ||
      this._zoom_end[1] !== this._zoom_start[1]) {
    var dslab = this._zoom_camera(eye);
    if (dslab) {
      this.slab_width[0] = Math.max(this.slab_width[0] + dslab, 0.01);
      this.slab_width[1] = Math.max(this.slab_width[1] + dslab, 0.01);
    }
    changed = true;
  }
  if (this._pan_end[0] !== this._pan_start[0] ||
      this._pan_end[1] !== this._pan_start[1]) {
    this._pan_camera(eye);
    this._panned = true;
    changed = true;
  }
  this._camera.position.addVectors(this._target, eye);
  if (this._state === STATE.GO && this._go_func) {
    this._go_func();
    changed = true;
  }
  this._camera.lookAt(this._target);
  return changed;
};

Controls.prototype.start = function start (new_state /*:number*/, x /*:number*/, y /*:number*/, dist/*:?number*/) {
  if (this._state === STATE.NONE || this._state === STATE.AUTO_ROTATE) {
    this._state = new_state;
  }
  this.move(x, y, dist);
  switch (this._state) {
    case STATE.ROTATE:
      this._rotate_start.copy(this._rotate_end);
      break;
    case STATE.ZOOM:
    case STATE.SLAB:
    case STATE.ROLL:
      this._zoom_start[0] = this._zoom_end[0];
      this._zoom_start[1] = this._zoom_end[1];
      break;
    case STATE.PAN:
      this._pan_start[0] = this._pan_end[0];
      this._pan_start[1] = this._pan_end[1];
      this._panned = false;
      break;
    case STATE.PAN_ZOOM:
      this._pinch_start = this._pinch_end;
      this._pan_start[0] = this._pan_end[0];
      this._pan_start[1] = this._pan_end[1];
      break;
  }
};

Controls.prototype.move = function move (x /*:number*/, y /*:number*/, dist /*:?number*/) {
  switch (this._state) {
    case STATE.ROTATE: {
      var xyz = project_on_ball(x, y);
      //console.log(this._camera.projectionMatrix);
      //console.log(this._camera.matrixWorld);
      // TODO maybe use project()/unproject()/applyProjection()
      var eye = this._camera.position.clone().sub(this._target);
      var up = this._camera.up;
      this._rotate_end.crossVectors(up, eye).setLength(xyz[0]);
      this._rotate_end.addScaledVector(up, xyz[1] / up.length());
      this._rotate_end.addScaledVector(eye, xyz[2] / eye.length());
      break;
    }
    case STATE.ZOOM:
    case STATE.SLAB:
    case STATE.ROLL:
      this._zoom_end = [x, y];
      break;
    case STATE.PAN:
      this._pan_end = [x, y];
      break;
    case STATE.PAN_ZOOM:
      if (dist == null) { return; } // should not happen
      this._pan_end = [x, y];
      this._pinch_end = dist;
      break;
  }
};

// returned coordinates can be used for atom picking
Controls.prototype.stop = function stop () {
  var ret = null;
  if (this._state === STATE.PAN && !this._panned) { ret = this._pan_end; }
  this._state = STATE.NONE;
  this._rotate_start.copy(this._rotate_end);
  this._pinch_start = this._pinch_end;
  this._pan_start[0] = this._pan_end[0];
  this._pan_start[1] = this._pan_end[1];
  return ret;
};

// cam_up (if set) must be orthogonal to the view
Controls.prototype.go_to = function go_to (targ /*:Vector3*/, cam_pos /*:?Vector3*/, cam_up /*:?Vector3*/,
      steps /*:?number*/) {
  if ((!targ || targ.distanceToSquared(this._target) < 0.001) &&
      (!cam_pos || cam_pos.distanceToSquared(this._camera.position) < 0.1) &&
      (!cam_up || cam_up.distanceToSquared(this._camera.up) < 0.1)) {
    return;
  }
  this._state = STATE.GO;
  steps = (steps || 60) / auto_speed;
  var alphas = [];
  var prev_pos = 0;
  for (var i = 1; i <= steps; ++i) {
    var pos = i / steps;
    // quadratic easing
    pos = pos < 0.5 ? 2 * pos * pos : -2 * pos * (pos-2) - 1;
    alphas.push((pos - prev_pos) / (1 - prev_pos));
    prev_pos = pos;
  }
  this._go_func = function () {
    var a = alphas.shift();
    if (targ) {
      // unspecified cam_pos - _camera stays in the same distance to _target
      if (!cam_pos) { this._camera.position.sub(this._target); }
      this._target.lerp(targ, a);
      if (!cam_pos) { this._camera.position.add(this._target); }
    }
    if (cam_pos) { this._camera.position.lerp(cam_pos, a); }
    if (cam_up) { this._camera.up.lerp(cam_up, a); }
    if (alphas.length === 0) {
      this._state = STATE.NONE;
      this._go_func = null;
    }
  };
};

// @flow

/*::
 import type {AtomT, Model} from './model.js'
 import type {Mesh} from './fromthree.js'

 type ColorScheme = {
   name: string,
   bg: number,
   fg: number,
   [name:string]: number | number[],
 };
 type Num2 = [number, number]
 type Num3 = [number, number, number];
 */

var ColorSchemes /*:ColorScheme[]*/ = [ // Viewer.prototype.ColorSchemes
  { // generally mimicks Coot
    name: 'coot dark',
    bg: 0x000000,
    fg: 0xFFFFFF,
    map_den: 0x3362B2,
    map_pos: 0x298029,
    map_neg: 0x8B2E2E,
    center: 0xC997B0,
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
    def: 0xa0a0a0, // default atom color
  },
  // scheme made of "solarized" colors (http://ethanschoonover.com/solarized):
  // base03  base02  base01  base00  base0   base1   base2   base3
  // #002b36 #073642 #586e75 #657b83 #839496 #93a1a1 #eee8d5 #fdf6e3
  // yellow  orange  red     magenta violet  blue    cyan    green
  // #b58900 #cb4b16 #dc322f #d33682 #6c71c4 #268bd2 #2aa198 #859900
  {
    name: 'solarized dark',
    bg: 0x002b36,
    fg: 0xfdf6e3,
    map_den: 0x268bd2,
    map_pos: 0x859900,
    map_neg: 0xd33682,
    center: 0xfdf6e3,
    H: 0x586e75,
    C: 0x93a1a1,
    N: 0x6c71c4,
    O: 0xcb4b16,
    S: 0xb58900,
    def: 0xeee8d5,
  },
  {
    name: 'solarized light',
    bg: 0xfdf6e3,
    fg: 0x002b36,
    map_den: 0x268bd2,
    map_pos: 0x859900,
    map_neg: 0xd33682,
    center: 0x002b36,
    H: 0x93a1a1,
    C: 0x586e75,
    N: 0x6c71c4,
    O: 0xcb4b16,
    S: 0xb58900,
    def: 0x073642,
  },
  { // like in Coot after Edit > Background Color > White
    name: 'coot light',
    bg: 0xFFFFFF,
    fg: 0x000000,
    map_den: 0x3362B2,
    map_pos: 0x298029,
    map_neg: 0x8B2E2E,
    center: 0xC7C769,
    H: 0x999999,
    C: 0xA96464,
    N: 0x1C51B3,
    O: 0xC33869,
    S: 0x9E7B3D,
    def: 0x808080,
  } ];

var INIT_HUD_TEXT = 'This is UglyMol not Coot. ' +
  '<a href="#" onclick="V.toggle_help(); return false;">H shows help.</a>';

// options handled by select_next()

var COLOR_PROPS = ['element', 'B-factor', 'occupancy', 'index', 'chain'];
var RENDER_STYLES = ['lines', 'trace', 'ribbon', 'ball&stick'];
var LIGAND_STYLES = ['ball&stick', 'lines'];
var WATER_STYLES = ['cross', 'dot', 'invisible'];
var MAP_STYLES = ['marching cubes', 'squarish' ];
var LINE_STYLES = ['normal', 'simplistic'];
var LABEL_FONTS = ['bold 14px', '14px', '16px', 'bold 16px'];

function rainbow_value(v/*:number*/, vmin/*:number*/, vmax/*:number*/) {
  var c = new Color(0xe0e0e0);
  if (vmin < vmax) {
    var ratio = (v - vmin) / (vmax - vmin);
    var hue = (240 - (240 * ratio)) / 360;
    c.setHSL(hue, 1.0, 0.5);
  }
  return c;
}

function color_by(prop, atoms /*:AtomT[]*/, elem_colors, hue_shift) {
  var color_func;
  var last_atom = atoms[atoms.length-1];
  if (prop === 'index') {
    color_func = function (atom) {
      return rainbow_value(atom.i_seq, 0, last_atom.i_seq);
    };
  } else if (prop === 'B-factor') {
    var vmin = Infinity;
    var vmax = -Infinity;
    for (var i = 0; i < atoms.length; i++) {
      var v = atoms[i].b;
      if (v > vmax) { vmax = v; }
      if (v < vmin) { vmin = v; }
    }
    //console.log('B-factors in [' + vmin + ', ' + vmax + ']');
    color_func = function (atom) {
      return rainbow_value(atom.b, vmin, vmax);
    };
  } else if (prop === 'occupancy') {
    color_func = function (atom) {
      return rainbow_value(atom.occ, 0, 1);
    };
  } else if (prop === 'chain') {
    color_func = function (atom) {
      return rainbow_value(atom.chain_index, 0, last_atom.chain_index);
    };
  } else { // element
    if (hue_shift === 0) {
      color_func = function (atom) {
        return elem_colors[atom.element] || elem_colors.def;
      };
    } else {
      var c_hsl = elem_colors['C'].getHSL();
      var c_col = new Color(0, 0, 0);
      c_col.setHSL(c_hsl.h + hue_shift, c_hsl.s, c_hsl.l);
      color_func = function (atom) {
        var el = atom.element;
        return el === 'C' ? c_col : (elem_colors[el] || elem_colors.def);
      };
    }
  }
  return atoms.map(color_func);
}

function scale_by_height(value, size) { // for scaling bond_line
  return value * size[1] / 700;
}

var MapBag = function MapBag(map/*:ElMap*/, config/*:Object*/, is_diff_map/*:boolean*/) {
  this.map = map;
  this.name = '';
  this.isolevel = is_diff_map ? 3.0 : config.default_isolevel;
  this.visible = true;
  this.types = is_diff_map ? ['map_pos', 'map_neg'] : ['map_den'];
  this.block_ctr = new Vector3(Infinity, 0, 0);
  this.el_objects = []; // three.js objects
};


var ModelBag = function ModelBag(model/*:Model*/, config/*:Object*/, win_size/*:Num2*/) {
  this.model = model;
  this.label = '(model #' + ++ModelBag.ctor_counter + ')';
  this.visible = true;
  this.hue_shift = 0;
  this.conf = config;
  this.win_size = win_size;
  this.objects = []; // list of three.js objects
  this.atom_array = [];
};

ModelBag.prototype.get_visible_atoms = function get_visible_atoms () {
  var atoms = this.model.atoms;
  if (this.conf.hydrogens || !this.model.has_hydrogens) {
    return atoms;
  }
  // with filter() it's twice slower (on Node 4.2)
  //return atoms.filter(function(a) { return a.element !== 'H'; });
  var non_h = [];
  for (var i = 0, list = atoms; i < list.length; i += 1) {
    var atom = list[i];

      if (atom.element !== 'H') { non_h.push(atom); }
  }
  return non_h;
};

ModelBag.prototype.add_bonds = function add_bonds (polymers/*:boolean*/, ligands/*:boolean*/, ball_size/*:?number*/) {
  var visible_atoms = this.get_visible_atoms();
  var colors = color_by(this.conf.color_prop, visible_atoms,
                          this.conf.colors, this.hue_shift);
  var vertex_arr /*:Vector3[]*/ = [];
  var color_arr = [];
  var sphere_arr = [];
  var sphere_color_arr = [];
  var hydrogens = this.conf.hydrogens;
  for (var i = 0; i < visible_atoms.length; i++) {
    var atom = visible_atoms[i];
    var color = colors[i];
    if (!(atom.is_ligand ? ligands : polymers)) { continue; }
    if (atom.is_water() && this.conf.water_style === 'invisible') { continue; }
    if (atom.bonds.length === 0 && ball_size == null) { // nonbonded - cross
      if (!atom.is_water() || this.conf.water_style === 'cross') {
        addXyzCross(vertex_arr, atom.xyz, 0.7);
        for (var n = 0; n < 6; n++) {
          color_arr.push(color);
        }
      }
    } else { // bonded, draw lines
      for (var j = 0; j < atom.bonds.length; j++) {
        var other = this.model.atoms[atom.bonds[j]];
        if (!hydrogens && other.element === 'H') { continue; }
        // Coot show X-H bonds as thinner lines in a single color.
        // Here we keep it simple and render such bonds like all others.
        var mid = atom.midpoint(other);
        vertex_arr.push(atom.xyz, mid);
        color_arr.push(color, color);
      }
    }
    sphere_arr.push(atom);
    sphere_color_arr.push(color);
  }

  if (ball_size != null) {
    if (vertex_arr.length !== 0) {
      this.objects.push(makeSticks(vertex_arr, color_arr, ball_size / 2));
    }
    if (sphere_arr.length !== 0) {
      this.objects.push(makeBalls(sphere_arr, sphere_color_arr, ball_size));
    }
  } else if (vertex_arr.length !== 0) {
    var linewidth = scale_by_height(this.conf.bond_line, this.win_size);
    var material = makeLineMaterial({
      linewidth: linewidth,
      win_size: this.win_size,
      segments: true,
    });
    this.objects.push(makeLineSegments(material, vertex_arr, color_arr));
    if (this.conf.line_style !== 'simplistic') {
      // wheels (discs) as round caps
      this.objects.push(makeWheels(sphere_arr, sphere_color_arr, linewidth));
    }
  }

  sphere_arr.forEach(function (v) { this.atom_array.push(v); }, this);
};

ModelBag.prototype.add_trace = function add_trace () {
  var segments = this.model.extract_trace();
  var visible_atoms = [].concat.apply([], segments);
  var colors = color_by(this.conf.color_prop, visible_atoms,
                          this.conf.colors, this.hue_shift);
  var material = makeLineMaterial({
    linewidth: scale_by_height(this.conf.bond_line, this.win_size),
    win_size: this.win_size,
  });
  var k = 0;
  for (var i$1 = 0, list$1 = segments; i$1 < list$1.length; i$1 += 1) {
    var seg = list$1[i$1];

      var color_slice = colors.slice(k, k + seg.length);
    k += seg.length;
    var pos = [];
    for (var i = 0, list = seg; i < list.length; i += 1) {
      var atom = list[i];

        pos.push(atom.xyz);
    }
    var line = makeLine(material, pos, color_slice);
    this.objects.push(line);
  }
  this.atom_array = visible_atoms;
};

ModelBag.prototype.add_ribbon = function add_ribbon (smoothness/*:number*/) {
  var segments = this.model.extract_trace();
  var res_map = this.model.get_residues();
  var visible_atoms = [].concat.apply([], segments);
  var colors = color_by(this.conf.color_prop, visible_atoms,
                          this.conf.colors, this.hue_shift);
  var k = 0;
  for (var i$1 = 0, list$1 = segments; i$1 < list$1.length; i$1 += 1) {
    var seg = list$1[i$1];

      var tangents = [];
    var last = [0, 0, 0];
    for (var i = 0, list = seg; i < list.length; i += 1) {
      var atom = list[i];

        var residue = res_map[atom.resid()];
      var tang = this.model.calculate_tangent_vector(residue);
      // untwisting (usually applies to beta-strands)
      if (tang[0]*last[0] + tang[1]*last[1] + tang[2]*last[2] < 0) {
        tang[0] = -tang[0];
        tang[1] = -tang[1];
        tang[2] = -tang[2];
      }
      tangents.push(tang);
      last = tang;
    }
    var color_slice = colors.slice(k, k + seg.length);
    k += seg.length;
    var obj = makeRibbon(seg, color_slice, tangents, smoothness);
    this.objects.push(obj);
  }
};

ModelBag.ctor_counter = 0;

function vec3_to_fixed(vec, n) {
  return [vec.x.toFixed(n), vec.y.toFixed(n), vec.z.toFixed(n)];
}

// for two-finger touch events
function touch_info(evt/*:TouchEvent*/) {
  var touches = evt.touches;
  var dx = touches[0].pageX - touches[1].pageX;
  var dy = touches[0].pageY - touches[1].pageY;
  return {pageX: (touches[0].pageX + touches[1].pageX) / 2,
          pageY: (touches[0].pageY + touches[1].pageY) / 2,
          dist: Math.sqrt(dx * dx + dy * dy)};
}

// makes sense only for full-window viewer
function parse_url_fragment() {
  var ret = {};
  if (typeof window === 'undefined') { return ret; }
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


var Viewer = function Viewer(options /*: {[key: string]: any}*/) {
  // rendered objects
  this.model_bags = [];
  this.map_bags = [];
  this.decor = { cell_box: null, selection: null, zoom_grid: makeGrid(),
                 mark: null };
  this.labels = {};
  this.nav = null;
  this.xhr_headers = {};

  this.config = {
    bond_line: 4.0, // ~ to height, like in Coot (see scale_by_height())
    map_line: 1.25,// for any height
    map_radius: 10.0,
    max_map_radius: 40,
    default_isolevel: 1.5,
    center_cube_size: 0.1,
    map_style: MAP_STYLES[0],
    render_style: RENDER_STYLES[0],
    ligand_style: LIGAND_STYLES[0],
    water_style: WATER_STYLES[0],
    color_prop: COLOR_PROPS[0],
    line_style: LINE_STYLES[0],
    label_font: LABEL_FONTS[0],
    colors: this.ColorSchemes[0],
    hydrogens: false,
    ball_size: 0.4,
  };

  // options of the constructor overwrite default values of the config
  for (var i = 0, list = options; i < list.length; i += 1) {
    var o = list[i];

    if (o in this.config) {
      this.config[o] = options[o];
    }
  }

  this.set_colors();
  this.window_size = [1, 1]; // it will be set in resize()
  this.window_offset = [0, 0];

  this.last_ctr = new Vector3(Infinity, 0, 0);
  this.selected = {bag: null, atom: null};
  this.dbl_click_callback = this.toggle_label;
  this.scene = new Scene();
  this.scene.fog = new Fog(this.config.colors.bg, 0, 1);
  this.scene.add(new AmbientLight(0xffffff));
  this.default_camera_pos = [0, 0, 100];
  if (options.share_view) {
    this.target = options.share_view.target;
    this.camera = options.share_view.camera;
    this.controls = options.share_view.controls;
    this.tied_viewer = options.share_view;
    this.tied_viewer.tied_viewer = this; // not GC friendly
  } else {
    this.target = new Vector3(0, 0, 0);
    this.camera = new OrthographicCamera();
    this.camera.position.fromArray(this.default_camera_pos);
    this.controls = new Controls(this.camera, this.target);
  }
  this.set_common_key_bindings();
  if (this.constructor === Viewer) { this.set_real_space_key_bindings(); }
  if (typeof document === 'undefined') { return; }// for testing on node

  function get_elem(name) {
    if (options[name] === null) { return null; }
    return document.getElementById(options[name] || name);
  }
  this.hud_el = get_elem('hud');

  try {
    this.renderer = new WebGLRenderer({antialias: true});
  } catch (e) {
    this.hud('No WebGL in your browser?', 'ERR');
    this.renderer = null;
    return;
  }

  this.container = get_elem('viewer');
  this.help_el = get_elem('help');
  if (this.hud_el) {
    if (this.hud_el.innerHTML === '') { this.hud_el.innerHTML = INIT_HUD_TEXT; }
    this.initial_hud_html = this.hud_el.innerHTML;
  }

  if (this.container == null) { return; } // can be null in headless tests
  this.renderer.setClearColor(this.config.colors.bg, 1);
  this.renderer.setPixelRatio(window.devicePixelRatio);
  this.resize();
  this.camera.zoom = this.camera.right / 35.0;// arbitrary choice
  this.update_camera();
  var el = this.renderer.domElement;
  // $FlowFixMe: flow can't figure out that this.container != null
  this.container.appendChild(el);
  if (options.focusable) {
    el.tabIndex = 0;
  }
  this.decor.zoom_grid.visible = false;
  this.scene.add(this.decor.zoom_grid);

  window.addEventListener('resize', this.resize.bind(this));
  var keydown_el = (options.focusable ? el : window);
  keydown_el.addEventListener('keydown', this.keydown.bind(this));
  el.addEventListener('contextmenu', function (e) { e.preventDefault(); });
  el.addEventListener('mousewheel', this.mousewheel.bind(this));
  el.addEventListener('MozMousePixelScroll', this.mousewheel.bind(this));
  el.addEventListener('mousedown', this.mousedown.bind(this));
  el.addEventListener('touchstart', this.touchstart.bind(this));
  el.addEventListener('touchmove', this.touchmove.bind(this));
  el.addEventListener('touchend', this.touchend.bind(this));
  el.addEventListener('touchcancel', this.touchend.bind(this));
  el.addEventListener('dblclick', this.dblclick.bind(this));

  var self = this;

  this.mousemove = function (event/*:MouseEvent*/) {
    event.preventDefault();
    //event.stopPropagation();
    self.controls.move(self.relX(event), self.relY(event));
  };

  this.mouseup = function (event/*:MouseEvent*/) {
    event.preventDefault();
    event.stopPropagation();
    document.removeEventListener('mousemove', self.mousemove);
    document.removeEventListener('mouseup', self.mouseup);
    self.decor.zoom_grid.visible = false;
    var not_panned = self.controls.stop();
    // special case - centering on atoms after action 'pan' with no shift
    if (not_panned) {
      var pick = self.pick_atom(not_panned, self.camera);
      if (pick != null) {
        self.select_atom(pick, {steps: 60});
      }
    }
    self.redraw_maps();
  };

  this.scheduled = false;
  this.request_render();
};

Viewer.prototype.pick_atom = function pick_atom (coords/*:Num2*/, camera/*:OrthographicCamera*/) {
  var pick = null;
  for (var i$1 = 0, list$1 = this.model_bags; i$1 < list$1.length; i$1 += 1) {
    var bag = list$1[i$1];

      if (!bag.visible) { continue; }
    var z = (camera.near + camera.far) / (camera.near - camera.far);
    var ray = new Ray();
    ray.origin.set(coords[0], coords[1], z).unproject(camera);
    ray.direction.set(0, 0, -1).transformDirection(camera.matrixWorld);
    var near = camera.near;
    // '0.15' b/c the furthest 15% is hardly visible in the fog
    var far = camera.far - 0.15 * (camera.far - camera.near);
    /*
    // previous version - line-based search
    let intersects = [];
    for (const object of bag.objects) {
      if (object.visible === false) continue;
      if (object.userData.bond_lines) {
        line_raycast(object, {ray, near, far, precision: 0.3}, intersects);
      }
    }
    ...
    if (intersects.length > 0) {
      intersects.sort(function (x) { return x.dist2 || Infinity; });
      const p = intersects[0].point;
      const atom = bag.model.get_nearest_atom(p.x, p.y, p.z);
      if (atom != null) {
        return {bag, atom};
      }
    }
    */
    // search directly atom array ignoring matrixWorld
    var vec = new Vector3();
    // required picking precision: 0.35A at zoom 50, 0.27A @z30, 0.44 @z80
    var precision2 = 0.35 * 0.35 * 0.02 * camera.zoom;
    for (var i = 0, list = bag.atom_array; i < list.length; i += 1) {
      var atom = list[i];

        vec.set(atom.xyz[0] - ray.origin.x,
              atom.xyz[1] - ray.origin.y,
              atom.xyz[2] - ray.origin.z);
      var distance = vec.dot(ray.direction);
      if (distance < 0 || distance < near || distance > far) { continue; }
      var diff2 = vec.addScaledVector(ray.direction, -distance).lengthSq();
      if (diff2 > precision2) { continue; }
      if (pick == null || distance < pick.distance) {
        pick = {bag: bag, atom: atom, distance: distance};
      }
    }
  }
  return pick;
};

Viewer.prototype.set_colors = function set_colors (scheme/*:?number|string|ColorScheme*/) {
  function to_col(x) { return new Color(x); }
  if (scheme == null) {
    scheme = this.config.colors;
  } else if (typeof scheme === 'number') {
    scheme = this.ColorSchemes[scheme % this.ColorSchemes.length];
  } else if (typeof scheme === 'string') {
    for (var i = 0, list = this.ColorSchemes; i < list.length; i += 1) {
      var sc = list[i];

        if (sc.name === scheme) {
        scheme = sc;
        break;
      }
    }
    throw Error('Unknown color scheme.');
  }
  if (typeof scheme.bg === 'number') {
    for (var key in scheme) {
      if (key !== 'name') {
        scheme[key] = scheme[key] instanceof Array ? scheme[key].map(to_col)
                                                   : to_col(scheme[key]);
      }
    }
  }
  this.decor.zoom_grid.color_value.set(scheme.fg);
  this.config.config = scheme;
  this.redraw_all();
};

// relative position on canvas in normalized device coordinates [-1, +1]
Viewer.prototype.relX = function relX (evt/*:{pageX: number}*/) {
  return 2 * (evt.pageX - this.window_offset[0]) / this.window_size[0] - 1;
};

Viewer.prototype.relY = function relY (evt/*:{pageY: number}*/) {
  return 1 - 2 * (evt.pageY - this.window_offset[1]) / this.window_size[1];
};

Viewer.prototype.hud = function hud (text/*:?string*/, type/*:?string*/) {
  if (typeof document === 'undefined') { return; }// for testing on node
  var el = this.hud_el;
  if (el) {
    if (text != null) {
      if (type === 'HTML') {
        el.innerHTML = text;
      } else {
        el.textContent = text;
      }
    } else {
      el.innerHTML = this.initial_hud_html;
    }
    var err = (type === 'ERR');
    el.style.backgroundColor = (err ? '#b00' : '');
    if (err && text) { console.log('ERR: ' + text); }
  } else {
    console.log('hud:', text);
  }
};

Viewer.prototype.redraw_center = function redraw_center (force/*:?boolean*/) {
  var size = this.config.center_cube_size;
  if (force ||
      this.target.distanceToSquared(this.last_ctr) > 0.01 * size * size) {
    this.last_ctr.copy(this.target);
    if (this.decor.mark) {
      this.scene.remove(this.decor.mark);
    }
    this.decor.mark = makeCube(size, this.target, {
      color: this.config.colors.center,
      linewidth: 2,
      win_size: this.window_size,
    });
    this.scene.add(this.decor.mark);
  }
};

Viewer.prototype.redraw_maps = function redraw_maps (force/*:?boolean*/) {
  this.redraw_center(force);
  var r = this.config.map_radius;
  for (var i = 0, list = this.map_bags; i < list.length; i += 1) {
    var map_bag = list[i];

      if (force || this.target.distanceToSquared(map_bag.block_ctr) > r/100) {
      this.redraw_map(map_bag);
    }
  }
};

Viewer.prototype.remove_and_dispose = function remove_and_dispose (obj/*:Object*/) {
  this.scene.remove(obj);
  if (obj.geometry) { obj.geometry.dispose(); }
  if (obj.material) {
    if (obj.material.uniforms && obj.material.uniforms.map) {
      obj.material.uniforms.map.value.dispose();
    }
    obj.material.dispose();
  }
  for (var i = 0, list = obj.children; i < list.length; i += 1) {
    var o = list[i];

      this.remove_and_dispose(o);
  }
};

Viewer.prototype.clear_el_objects = function clear_el_objects (map_bag/*:MapBag*/) {
  for (var i = 0, list = map_bag.el_objects; i < list.length; i += 1) {
    var o = list[i];

      this.remove_and_dispose(o);
  }
  map_bag.el_objects = [];
};

Viewer.prototype.clear_model_objects = function clear_model_objects (model_bag/*:ModelBag*/) {
  for (var i = 0, list = model_bag.objects; i < list.length; i += 1) {
    var o = list[i];

      this.remove_and_dispose(o);
  }
  model_bag.objects = [];
};

Viewer.prototype.has_frag_depth = function has_frag_depth () {
  return this.renderer && this.renderer.extensions.get('EXT_frag_depth');
};

Viewer.prototype.set_model_objects = function set_model_objects (model_bag/*:ModelBag*/) {
  model_bag.objects = [];
  model_bag.atom_array = [];
  var ligand_balls = null;
  if (model_bag.conf.ligand_style === 'ball&stick' && this.has_frag_depth()) {
    ligand_balls = this.config.ball_size;
  }
  switch (model_bag.conf.render_style) {
    case 'lines':
      if (ligand_balls === null) {
        model_bag.add_bonds(true, true);
      } else {
        model_bag.add_bonds(true, false);
        model_bag.add_bonds(false, true, ligand_balls);
      }
      break;
    case 'ball&stick':
      if (!this.has_frag_depth()) {
        this.hud('Ball-and-stick rendering is not working in this browser' +
                 '\ndue to lack of suppport for EXT_frag_depth', 'ERR');
        return;
      }
      if (ligand_balls === null) {
        model_bag.add_bonds(true, false, this.config.ball_size);
        model_bag.add_bonds(false, true);
      } else {
        model_bag.add_bonds(true, true, this.config.ball_size);
      }
      break;
    case 'trace':
      model_bag.add_trace();
      model_bag.add_bonds(false, true, ligand_balls);
      break;
    case 'ribbon':
      model_bag.add_ribbon(8);
      model_bag.add_bonds(false, true, ligand_balls);
      break;
  }
  for (var i = 0, list = model_bag.objects; i < list.length; i += 1) {
    var o = list[i];

      this.scene.add(o);
  }
};

// Add/remove label if `show` is specified, toggle otherwise.
Viewer.prototype.toggle_label = function toggle_label (pick/*:{bag:?ModelBag, atom:?AtomT}*/, show/*:?boolean*/) {
  if (pick.atom == null) { return; }
  var text = pick.atom.short_label();
  var uid = text; // we assume that the labels inside one model are unique
  var is_shown = (uid in this.labels);
  if (show === undefined) { show = !is_shown; }
  if (show) {
    if (is_shown) { return; }
    if (pick.atom == null) { return; } // silly flow
    var atom_style = pick.atom.is_ligand ? 'ligand_style' : 'render_style';
    var balls = pick.bag && pick.bag.conf[atom_style] === 'ball&stick';
    var label = makeLabel(text, {
      pos: pick.atom.xyz,
      font: this.config.label_font,
      color: '#' + this.config.colors.fg.getHexString(),
      win_size: this.window_size,
      z_shift: balls ? this.config.ball_size + 0.1 : 0.2,
    });
    if (!label) { return; }
    if (pick.bag == null) { return; }
    this.labels[uid] = { o: label, bag: pick.bag };
    this.scene.add(label);
  } else {
    if (!is_shown) { return; }
    this.remove_and_dispose(this.labels[uid].o);
    delete this.labels[uid];
  }
};

Viewer.prototype.redraw_labels = function redraw_labels () {
  for (var uid in this.labels) { // eslint-disable-line guard-for-in
    var text = uid;
    this.labels[uid].o.remake(text, {
      font: this.config.label_font,
      color: '#' + this.config.colors.fg.getHexString(),
    });
  }
};

Viewer.prototype.toggle_map_visibility = function toggle_map_visibility (map_bag/*:MapBag*/) {
  if (typeof map_bag === 'number') {
    map_bag = this.map_bags[map_bag];
  }
  map_bag.visible = !map_bag.visible;
  this.redraw_map(map_bag);
  this.request_render();
};

Viewer.prototype.redraw_map = function redraw_map (map_bag/*:MapBag*/) {
  this.clear_el_objects(map_bag);
  if (map_bag.visible) {
    map_bag.map.block.clear();
    this.add_el_objects(map_bag);
  }
};

Viewer.prototype.toggle_model_visibility = function toggle_model_visibility (model_bag/*:?ModelBag*/, visible/*:?boolean*/) {
  model_bag = model_bag || this.selected.bag;
  if (model_bag == null) { return; }
  model_bag.visible = visible == null ? !model_bag.visible : visible;
  this.redraw_model(model_bag);
  this.request_render();
};

Viewer.prototype.redraw_model = function redraw_model (model_bag/*:ModelBag*/) {
  this.clear_model_objects(model_bag);
  if (model_bag.visible) {
    this.set_model_objects(model_bag);
  }
};

Viewer.prototype.redraw_models = function redraw_models () {
  for (var i = 0, list = this.model_bags; i < list.length; i += 1) {
    var model_bag = list[i];

      this.redraw_model(model_bag);
  }
};

Viewer.prototype.add_el_objects = function add_el_objects (map_bag/*:MapBag*/) {
  if (!map_bag.visible || this.config.map_radius <= 0) { return; }
  if (map_bag.map.block.empty()) {
    var t = this.target;
    map_bag.block_ctr.copy(t);
    map_bag.map.extract_block(this.config.map_radius, [t.x, t.y, t.z]);
  }
  for (var i = 0, list = map_bag.types; i < list.length; i += 1) {
    var mtype = list[i];

      var isolevel = (mtype === 'map_neg' ? -1 : 1) * map_bag.isolevel;
    var iso = map_bag.map.isomesh_in_block(isolevel, this.config.map_style);
    if (iso == null) { continue; }
    var obj = makeChickenWire(iso, {
      color: this.config.colors[mtype],
      linewidth: this.config.map_line,
    });
    map_bag.el_objects.push(obj);
    this.scene.add(obj);
  }
};

Viewer.prototype.change_isolevel_by = function change_isolevel_by (map_idx/*:number*/, delta/*:number*/) {
  if (map_idx >= this.map_bags.length) { return; }
  var map_bag = this.map_bags[map_idx];
  map_bag.isolevel += delta;
  //TODO: move slow part into update()
  this.clear_el_objects(map_bag);
  this.add_el_objects(map_bag);
  var abs_level = map_bag.map.abs_level(map_bag.isolevel);
  var abs_text = abs_level.toFixed(4);
  var tied = this.tied_viewer;
  if (tied && map_idx < tied.map_bags.length) {
    var tied_bag = tied.map_bags[map_idx];
    // Should we tie by sigma or absolute level? Now it's sigma.
    tied_bag.isolevel = map_bag.isolevel;
    abs_text += ' / ' + tied_bag.map.abs_level(tied_bag.isolevel).toFixed(4);
    tied.clear_el_objects(tied_bag);
    tied.add_el_objects(tied_bag);
  }
  this.hud('map ' + (map_idx+1) + ' level =  ' + abs_text + ' ' +
           map_bag.map.unit + ' (' + map_bag.isolevel.toFixed(2) + ' rmsd)');
};

Viewer.prototype.change_map_radius = function change_map_radius (delta/*:number*/) {
  var rmax = this.config.max_map_radius;
  var cf = this.config;
  cf.map_radius = Math.min(Math.max(cf.map_radius + delta, 0), rmax);
  cf.map_radius = Math.round(cf.map_radius * 1e9) / 1e9;
  var info = 'map "radius": ' + cf.map_radius;
  if (cf.map_radius === rmax) { info += ' (max)'; }
  else if (cf.map_radius === 0) { info += ' (hidden maps)'; }
  if (this.map_bags.length === 0) { info += '\nNB: no map is loaded.'; }
  this.hud(info);
  this.redraw_maps(true);
};

Viewer.prototype.change_slab_width_by = function change_slab_width_by (delta/*:number*/) {
  var slab_width = this.controls.slab_width;
  slab_width[0] = Math.max(slab_width[0] + delta, 0.01);
  slab_width[1] = Math.max(slab_width[1] + delta, 0.01);
  this.update_camera();
  var final_width = this.camera.far - this.camera.near;
  this.hud('clip width: ' + final_width.toPrecision(3));
};

Viewer.prototype.change_zoom_by_factor = function change_zoom_by_factor (mult/*:number*/) {
  this.camera.zoom *= mult;
  this.update_camera();
  this.hud('zoom: ' + this.camera.zoom.toPrecision(3));
};

Viewer.prototype.change_bond_line = function change_bond_line (delta/*:number*/) {
  this.config.bond_line = Math.max(this.config.bond_line + delta, 0.1);
  this.redraw_models();
  this.hud('bond width: ' + scale_by_height(this.config.bond_line,
                                            this.window_size).toFixed(1));
};

Viewer.prototype.change_map_line = function change_map_line (delta/*:number*/) {
  this.config.map_line = Math.max(this.config.map_line + delta, 0.1);
  this.redraw_maps(true);
  this.hud('wireframe width: ' + this.config.map_line.toFixed(1));
};

Viewer.prototype.toggle_full_screen = function toggle_full_screen () {
  var d = document;
  // $FlowFixMe: Property mozFullScreenElement is missing in Document
  if (d.fullscreenElement || d.mozFullScreenElement ||
      // $FlowFixMe: Property webkitExitFullscreen is missing in Document
      d.webkitFullscreenElement || d.msFullscreenElement) {
    // $FlowFixMe: Property webkitExitFullscreen is missing in Document
    var ex = d.exitFullscreen || d.webkitExitFullscreen ||
    // $FlowFixMe: property `msExitFullscreen` not found in document
             d.mozCancelFullScreen || d.msExitFullscreen;
    // $FlowFixMe: cannot call property `exitFullscreen` of unknown type
    if (ex) { ex.call(d); }
  } else {
    var el = this.container;
    if (!el) { return; }
    // $FlowFixMe: Property webkitRequestFullscreen is missing in HTMLElement
    var req = el.requestFullscreen || el.webkitRequestFullscreen ||
    // $FlowFixMe: property `msRequestFullscreen` not found in HTMLElement
              el.mozRequestFullScreen || el.msRequestFullscreen;
    if (req) { req.call(el); }
  }
};

Viewer.prototype.toggle_cell_box = function toggle_cell_box () {
  if (this.decor.cell_box) {
    this.scene.remove(this.decor.cell_box);
    this.decor.cell_box = null;
  } else {
    var uc_func = this.get_cell_box_func();
    if (uc_func) {
      this.decor.cell_box = makeRgbBox(uc_func, this.config.colors.fg);
      this.scene.add(this.decor.cell_box);
    }
  }
};

Viewer.prototype.get_cell_box_func = function get_cell_box_func () /*:?Function*/ {
  var uc = null;
  if (this.selected.bag != null) {
    uc = this.selected.bag.model.unit_cell;
  }
  // note: model may not have unit cell
  if (uc == null && this.map_bags.length > 0) {
    uc = this.map_bags[0].map.unit_cell;
  }
  return uc && uc.orthogonalize.bind(uc);
};

Viewer.prototype.shift_clip = function shift_clip (delta/*:number*/) {
  var eye = this.camera.position.clone().sub(this.target);
  eye.multiplyScalar(delta / eye.length());
  this.target.add(eye);
  this.camera.position.add(eye);
  this.update_camera();
  this.redraw_maps();
  this.hud('clip shifted by [' + vec3_to_fixed(eye, 2).join(' ') + ']');
};

Viewer.prototype.go_to_nearest_Ca = function go_to_nearest_Ca () {
  var t = this.target;
  var bag = this.selected.bag;
  if (bag == null) { return; }
  var atom = bag.model.get_nearest_atom(t.x, t.y, t.z, 'CA');
  if (atom != null) {
    this.select_atom({bag: bag, atom: atom}, {steps: 30});
  } else {
    this.hud('no nearby CA');
  }
};

Viewer.prototype.toggle_inactive_models = function toggle_inactive_models () {
  var n = this.model_bags.length;
  if (n < 2) {
    this.hud((n == 0 ? 'No' : 'Only one') + ' model is loaded. ' +
             '"V" is for working with multiple models.');
    return;
  }
  var show_all = !this.model_bags.every(function (m) { return m.visible; });
  for (var i = 0, list = this.model_bags; i < list.length; i += 1) {
    var model_bag = list[i];

      var show = show_all || model_bag === this.selected.bag;
    this.toggle_model_visibility(model_bag, show);
  }
  this.hud(show_all ? 'All models visible' : 'Inactive models hidden');
};

Viewer.prototype.permalink = function permalink () {
  if (typeof window === 'undefined') { return; }
  var xyz_prec = Math.round(-Math.log10(0.001));
  window.location.hash =
    '#xyz=' + vec3_to_fixed(this.target, xyz_prec).join(',') +
    '&eye=' + vec3_to_fixed(this.camera.position, 1).join(',') +
    '&zoom=' + this.camera.zoom.toFixed(0);
  this.hud('copy URL from the location bar');
};

Viewer.prototype.redraw_all = function redraw_all () {
  if (!this.renderer) { return; }
  this.scene.fog.color = this.config.colors.bg;
  if (this.renderer) { this.renderer.setClearColor(this.config.colors.bg, 1); }
  this.redraw_models();
  this.redraw_maps(true);
  this.redraw_labels();
};

Viewer.prototype.toggle_help = function toggle_help () {
  var el = this.help_el;
  if (!el) { return; }
  el.style.display = el.style.display === 'block' ? 'none' : 'block';
  if (el.innerHTML === '') {
    el.innerHTML = [this.MOUSE_HELP, this.KEYBOARD_HELP,
                    this.ABOUT_HELP].join('\n\n');
  }
};

Viewer.prototype.select_next = function select_next (info/*:string*/, key/*:string*/,
            options/*:Array<string|ColorScheme>*/, back/*:boolean*/) {
  var old_idx = options.indexOf(this.config[key]);
  var len = options.length;
  var new_idx = (old_idx + (back ? len - 1 : 1)) % len;
  this.config[key] = options[new_idx];
  var html = info + ':';
  for (var i = 0; i < len; i++) {
    var tag = (i === new_idx ? 'u' : 's');
    var opt_name = typeof options[i] === 'string' ? options[i]
                                                    : options[i].name;
    html += ' <' + tag + '>' + opt_name + '</' + tag + '>';
  }
  this.hud(html, 'HTML');
};

Viewer.prototype.keydown = function keydown (evt/*:KeyboardEvent*/) {
  if (evt.ctrlKey) { return; }
  var action = this.key_bindings[evt.keyCode];
  if (action) {
    (action.bind(this))(evt);
  } else {
    if (action === false) { evt.preventDefault(); }
    if (this.help_el) { this.hud('Nothing here. Press H for help.'); }
  }
  this.request_render();
};

Viewer.prototype.set_common_key_bindings = function set_common_key_bindings () {
  var kb = new Array(256);
  // b
  kb[66] = function (evt) {
    this.select_next('color scheme', 'colors', this.ColorSchemes,
                     evt.shiftKey);
    this.set_colors();
  };
  // c
  kb[67] = function (evt) {
    this.select_next('coloring by', 'color_prop', COLOR_PROPS, evt.shiftKey);
    this.redraw_models();
  };
  // e
  kb[69] = function toggle_fog() {
    var fog = this.scene.fog;
    var has_fog = (fog.far === 1);
    fog.far = (has_fog ? 1e9 : 1);
    this.hud((has_fog ? 'dis': 'en') + 'able fog');
    this.redraw_all();
  };
  // h
  kb[72] = this.toggle_help;
  // i
  kb[73] = function (evt) {
    this.hud('toggled spinning');
    this.controls.toggle_auto(evt.shiftKey);
  };
  // k
  kb[75] = function () {
    this.hud('toggled rocking');
    this.controls.toggle_auto(0.0);
  };
  // m
  kb[77] = function (evt) {
    this.change_zoom_by_factor(evt.shiftKey ? 1.2 : 1.03);
  };
  // n
  kb[78] = function (evt) {
    this.change_zoom_by_factor(1 / (evt.shiftKey ? 1.2 : 1.03));
  };
  // q
  kb[81] = function (evt) {
    this.select_next('label font', 'label_font', LABEL_FONTS, evt.shiftKey);
    this.redraw_labels();
  };
  // r
  kb[82] = function (evt) {
    if (evt.shiftKey) {
      this.hud('redraw!');
      this.redraw_all();
    } else {
      this.hud('recentered');
      this.recenter();
    }
  };
  // w
  kb[87] = function (evt) {
    this.select_next('map style', 'map_style', MAP_STYLES, evt.shiftKey);
    this.redraw_maps(true);
  };
  // add, equals/firefox, equal sign
  kb[107] = kb[61] = kb[187] = function (evt) {
    this.change_isolevel_by(evt.shiftKey ? 1 : 0, 0.1);
  };
  // subtract, minus/firefox, dash
  kb[109] = kb[173] = kb[189] = function (evt) {
    this.change_isolevel_by(evt.shiftKey ? 1 : 0, -0.1);
  };
  // [
  kb[219] = function () { this.change_map_radius(-2); };
  // ]
  kb[221] = function () { this.change_map_radius(2); };
  // \ (backslash)
  kb[220] = function (evt) {
    this.select_next('bond lines', 'line_style', LINE_STYLES, evt.shiftKey);
    this.redraw_models();
  };
  // shift, ctrl, alt, altgr
  kb[16] = kb[17] = kb[18] = kb[225] = function () {};
  // slash, single quote
  kb[191] = kb[222] = false;// -> preventDefault()

  this.key_bindings = kb;
};

Viewer.prototype.set_real_space_key_bindings = function set_real_space_key_bindings () {
  var kb = this.key_bindings;
  // Home
  kb[36] = function (evt) {
    evt.shiftKey ? this.change_map_line(0.1) : this.change_bond_line(0.2);
  };
  // End
  kb[35] = function (evt) {
    evt.shiftKey ? this.change_map_line(-0.1) : this.change_bond_line(-0.2);
  };
  // Space
  kb[32] = function (evt) { this.center_next_residue(evt.shiftKey); };
  // d
  kb[68] = function () { this.change_slab_width_by(-0.1); };
  // f
  kb[70] = function (evt) {
    evt.shiftKey ? this.toggle_full_screen() : this.change_slab_width_by(0.1);
  };
  // l
  kb[76] = function (evt) {
    this.select_next('ligands as', 'ligand_style', LIGAND_STYLES,
                     evt.shiftKey);
    this.redraw_models();
  };
  // p
  kb[80] = function (evt) {
    evt.shiftKey ? this.permalink() : this.go_to_nearest_Ca();
  };
  // s
  kb[83] = function (evt) {
    this.select_next('rendering as', 'render_style', RENDER_STYLES,
                     evt.shiftKey);
    this.redraw_models();
  };
  // t
  kb[84] = function (evt) {
    this.select_next('waters as', 'water_style', WATER_STYLES, evt.shiftKey);
    this.redraw_models();
  };
  // u
  kb[85] = function () {
    this.hud('toggled unit cell box');
    this.toggle_cell_box();
  };
  // v
  kb[86] = function () { this.toggle_inactive_models(); };
  // y
  kb[89] = function (evt) {
    this.config.hydrogens = !this.config.hydrogens;
    this.hud((this.config.hydrogens ? 'show' : 'hide') +
             ' hydrogens (if any)');
    this.redraw_models();
  };
  // comma
  kb[188] = function (evt) { if (evt.shiftKey) { this.shift_clip(1); } };
  // period
  kb[190] = function (evt) { if (evt.shiftKey) { this.shift_clip(-1); } };
};

Viewer.prototype.mousedown = function mousedown (event/*:MouseEvent*/) {
  //event.preventDefault(); // default involves setting focus, which we need
  event.stopPropagation();
  document.addEventListener('mouseup', this.mouseup);
  document.addEventListener('mousemove', this.mousemove);
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
      this.decor.zoom_grid.visible = true;
      state = STATE.ZOOM;
    }
  }
  this.controls.start(state, this.relX(event), this.relY(event));
  this.request_render();
};

Viewer.prototype.dblclick = function dblclick (event/*:MouseEvent*/) {
  if (event.button !== 0) { return; }
  if (this.decor.selection) {
    this.remove_and_dispose(this.decor.selection);
    this.decor.selection = null;
  }
  var mouse = [this.relX(event), this.relY(event)];
  var pick = this.pick_atom(mouse, this.camera);
  if (pick) {
    var atom = pick.atom;
    this.hud(pick.bag.label + ' ' + atom.long_label());
    this.dbl_click_callback(pick);
    var color = this.config.colors[atom.element] || this.config.colors.def;
    var size = 2.5 * scale_by_height(this.config.bond_line,
                                       this.window_size);
    this.decor.selection = makeWheels([atom], [color], size);
    this.scene.add(this.decor.selection);
  } else {
    this.hud();
  }
  this.request_render();
};

Viewer.prototype.touchstart = function touchstart (event/*:TouchEvent*/) {
  var touches = event.touches;
  if (touches.length === 1) {
    this.controls.start(STATE.ROTATE,
                        this.relX(touches[0]), this.relY(touches[0]));
  } else { // for now using only two touches
    var info = touch_info(event);
    this.controls.start(STATE.PAN_ZOOM,
                        this.relX(info), this.relY(info), info.dist);
  }
  this.request_render();
};

Viewer.prototype.touchmove = function touchmove (event/*:TouchEvent*/) {
  event.preventDefault();
  event.stopPropagation();
  var touches = event.touches;
  if (touches.length === 1) {
    this.controls.move(this.relX(touches[0]), this.relY(touches[0]));
  } else { // for now using only two touches
    var info = touch_info(event);
    this.controls.move(this.relX(info), this.relY(info), info.dist);
  }
};

Viewer.prototype.touchend = function touchend (/*event*/) {
  this.controls.stop();
  this.redraw_maps();
};

// $FlowFixMe TODO: wheel()+WheelEvent are more standard
Viewer.prototype.mousewheel = function mousewheel (evt/*:MouseWheelEvent*/) {
  evt.preventDefault();
  evt.stopPropagation();
  // evt.wheelDelta for WebKit, evt.detail for Firefox
  var delta = evt.wheelDelta || -2 * (evt.detail || 0);
  this.mousewheel_action(delta, evt);
  this.request_render();
};

Viewer.prototype.mousewheel_action = function mousewheel_action (delta/*:number*/, evt/*:WheelEvent*/) {
  this.change_isolevel_by(evt.shiftKey ? 1 : 0, 0.0005 * delta);
};

Viewer.prototype.resize = function resize (/*evt*/) {
  var el = this.container;
  if (el == null) { return; }
  var width = el.clientWidth;
  var height = el.clientHeight;
  this.window_offset[0] = el.offsetLeft;
  this.window_offset[1] = el.offsetTop;
  this.camera.left = -width;
  this.camera.right = width;
  this.camera.top = height;
  this.camera.bottom = -height;
  this.camera.updateProjectionMatrix();
  this.renderer.setSize(width, height);
  if (width !== this.window_size[0] || height !== this.window_size[1]) {
    this.window_size[0] = width;
    this.window_size[1] = height;
    this.redraw_models(); // b/c bond_line is scaled by height
  }
  this.request_render();
};

// If xyz set recenter on it looking toward the model center.
// Otherwise recenter on the model center looking along the z axis.
Viewer.prototype.recenter = function recenter (xyz/*:?Num3*/, cam/*:?Num3*/, steps/*:?number*/) {
  var bag = this.selected.bag;
  var new_up = new Vector3(0, 1, 0);
  var eye;
  if (xyz != null && cam == null && bag != null) {
    // look from specified point toward the center of the molecule,
    // i.e. shift camera away from the molecule center.
    var mc = bag.model.get_center();
    eye = new Vector3(xyz[0] - mc[0], xyz[1] - mc[1], xyz[2] - mc[2]);
    eye.setLength(100);
    xyz = new Vector3(xyz[0], xyz[1], xyz[2]);
    cam = eye.clone().add(xyz);
  } else {
    if (xyz == null) {
      if (bag != null) {
        xyz = bag.model.get_center();
      } else {
        var uc_func = this.get_cell_box_func();
        xyz = uc_func ? uc_func([0.5, 0.5, 0.5]) : [0, 0, 0];
      }
    }
    xyz = new Vector3(xyz[0], xyz[1], xyz[2]);
    if (cam != null) {
      cam = new Vector3(cam[0], cam[1], cam[2]);
      eye = cam.clone().sub(xyz);
      new_up.copy(this.camera.up); // preserve the up direction
    } else {
      var dc = this.default_camera_pos;
      cam = new Vector3(xyz.x + dc[0], xyz.y + dc[1], xyz.z + dc[2]);
    }
  }
  if (eye != null) {
    new_up.projectOnPlane(eye);
    if (new_up.lengthSq() < 0.0001) { new_up.x += 1; }
    new_up.normalize();
  }
  this.controls.go_to(xyz, cam, new_up, steps);
};

Viewer.prototype.center_next_residue = function center_next_residue (back/*:boolean*/) {
  var bag = this.selected.bag;
  if (bag == null) { return; }
  var atom = bag.model.next_residue(this.selected.atom, back);
  if (atom != null) {
    this.select_atom({bag: bag, atom: atom}, {steps: 30});
  }
};

Viewer.prototype.select_atom = function select_atom (pick/*:{bag:ModelBag, atom:AtomT}*/, options) {
    if ( options === void 0 ) options/*:Object*/={};

  this.hud('-> ' + pick.bag.label + ' ' + pick.atom.long_label());
  var xyz = pick.atom.xyz;
  this.controls.go_to(new Vector3(xyz[0], xyz[1], xyz[2]),
                      null, null, options.steps);
  this.toggle_label(this.selected, false);
  this.selected = {bag: pick.bag, atom: pick.atom}; // not ...=pick b/c flow
  this.toggle_label(this.selected, true);
};

Viewer.prototype.update_camera = function update_camera () {
  var dxyz = this.camera.position.distanceTo(this.target);
  var w = this.controls.slab_width;
  var scale = w[2] || this.camera.zoom;
  this.camera.near = dxyz * (1 - w[0] / scale);
  this.camera.far = dxyz * (1 + w[1] / scale);
  this.camera.updateProjectionMatrix();
};

// The main loop. Running when a mouse button is pressed or when the view
// is moving (and run once more after the mouse button is released).
// It is also triggered by keydown events.
Viewer.prototype.render = function render () {
  this.scheduled = true;
  if (this.renderer === null) { return; }
  if (this.controls.update()) {
    this.update_camera();
  }
  var tied = this.tied_viewer;
  if (!this.controls.is_going()) {
    this.redraw_maps();
    if (tied && !tied.scheduled) { tied.redraw_maps(); }
  }
  this.renderer.render(this.scene, this.camera);
  if (tied && !tied.scheduled) { tied.renderer.render(tied.scene, tied.camera); }
  if (this.nav) {
    this.nav.renderer.render(this.nav.scene, this.camera);
  }
  this.scheduled = false;
  if (this.controls.is_moving()) {
    this.request_render();
  }
};

Viewer.prototype.request_render = function request_render () {
  if (typeof window !== 'undefined' && !this.scheduled) {
    this.scheduled = true;
    window.requestAnimationFrame(this.render.bind(this));
  }
};

Viewer.prototype.add_model = function add_model (model/*:Model*/, options) {
    if ( options === void 0 ) options/*:Object*/={};

  var model_bag = new ModelBag(model, this.config, this.window_size);
  model_bag.hue_shift = options.hue_shift || 0.06 * this.model_bags.length;
  this.model_bags.push(model_bag);
  this.set_model_objects(model_bag);
  this.request_render();
};

Viewer.prototype.add_map = function add_map (map/*:ElMap*/, is_diff_map/*:boolean*/) {
  //map.show_debug_info();
  var map_bag = new MapBag(map, this.config, is_diff_map);
  this.map_bags.push(map_bag);
  this.add_el_objects(map_bag);
  this.request_render();
};

Viewer.prototype.load_file = function load_file (url/*:string*/, options/*:{[id:string]: mixed}*/,
          callback/*:Function*/) {
  if (this.renderer === null) { return; }// no WebGL detected
  var req = new XMLHttpRequest();
  req.open('GET', url, true);
  if (options.binary) {
    req.responseType = 'arraybuffer';
  } else {
    // http://stackoverflow.com/questions/7374911/
    req.overrideMimeType('text/plain');
  }
  var self = this;
  Object.keys(this.xhr_headers).forEach(function (name) {
    req.setRequestHeader(name, self.xhr_headers[name]);
  });
  req.onreadystatechange = function () {
    if (req.readyState === 4) {
      // chrome --allow-file-access-from-files gives status 0
      if (req.status === 200 || (req.status === 0 && req.response !== null &&
                                                     req.response !== '')) {
        try {
          callback(req);
        } catch (e) {
          self.hud('Error: ' + e.message + '\nwhen processing ' + url, 'ERR');
        }
      } else {
        self.hud('Failed to fetch ' + url, 'ERR');
      }
    }
  };
  if (options.progress) {
    req.addEventListener('progress', function (evt /*:ProgressEvent*/) {
      if (evt.lengthComputable && evt.loaded && evt.total) {
        var fn = url.split('/').pop();
        self.hud('loading ' + fn + ' ... ' + ((evt.loaded / 1024) | 0) +
                 ' / ' + ((evt.total / 1024) | 0) + ' kB');
        if (evt.loaded === evt.total) { self.hud(); } // clear progress message
      }
    });
  }
  try {
    req.send(null);
  } catch (e) {
    self.hud('loading ' + url + ' failed:\n' + e, 'ERR');
  }
};

Viewer.prototype.set_dropzone = function set_dropzone (zone/*:Object*/, callback/*:Function*/) {
  var self = this;
  zone.addEventListener('dragover', function (e) {
    e.stopPropagation();
    e.preventDefault();
    e.dataTransfer.dropEffect = 'copy';
    self.hud('ready for file drop...');
  });
  zone.addEventListener('drop', function (e) {
    e.stopPropagation();
    e.preventDefault();
    var names = [];
    for (var i = 0, list = e.dataTransfer.files; i < list.length; i += 1) {
      var file = list[i];

        try {
        callback(file);
      } catch (e) {
        self.hud('Loading ' + file.name + ' failed.\n' + e.message, 'ERR');
        return;
      }
      names.push(file.name);
    }
    self.hud('loading ' + names.join(', '));
  });
};

// for use with set_dropzone
Viewer.prototype.pick_pdb_and_map = function pick_pdb_and_map (file/*:File*/) {
  var self = this;
  var reader = new FileReader();
  if (/\.(pdb|ent)$/.test(file.name)) {
    reader.onload = function (evt/*:any*/) {
      self.load_pdb_from_text(evt.target.result);
      self.recenter();
    };
    reader.readAsText(file);
  } else if (/\.(map|ccp4|mrc|dsn6|omap)$/.test(file.name)) {
    var map_format = /\.(dsn6|omap)$/.test(file.name) ? 'dsn6' : 'ccp4';
    reader.onloadend = function (evt/*:any*/) {
      if (evt.target.readyState == 2) {
        self.load_map_from_buffer(evt.target.result, {format: map_format});
        if (self.model_bags.length === 0 && self.map_bags.length === 1) {
          self.recenter();
        }
      }
    };
    reader.readAsArrayBuffer(file);
  } else {
    throw Error('Unknown file extension. ' +
                'Use: pdb, ent, ccp4, mrc, map, dsn6 or omap.');
  }
};

Viewer.prototype.set_view = function set_view (options/*:?Object*/) {
  var frag = parse_url_fragment();
  if (frag.zoom) { this.camera.zoom = frag.zoom; }
  this.recenter(frag.xyz || (options && options.center), frag.eye, 1);
};

// Load molecular model from PDB file and centers the view
Viewer.prototype.load_pdb_from_text = function load_pdb_from_text (text/*:string*/) {
  var len = this.model_bags.length;
  var models = modelsFromPDB(text);
  for (var i = 0, list = models; i < list.length; i += 1) {
    var model = list[i];

      this.add_model(model);
  }
  this.selected.bag = this.model_bags[len];
};

Viewer.prototype.load_pdb = function load_pdb (url/*:string*/, options/*:?Object*/, callback/*:?Function*/) {
  var self = this;
  this.load_file(url, {binary: false, progress: true}, function (req) {
    self.load_pdb_from_text(req.responseText);
    if (options == null || !options.stay) { self.set_view(options); }
    if (callback) { callback(); }
  });
};

Viewer.prototype.load_map = function load_map (url/*:?string*/, options/*:Object*/, callback/*:?Function*/) {
  if (url == null) {
    if (callback) { callback(); }
    return;
  }
  if (options.format !== 'ccp4' && options.format !== 'dsn6') {
    throw Error('Unknown map format.');
  }
  var self = this;
  this.load_file(url, {binary: true, progress: true}, function (req) {
    self.load_map_from_buffer(req.response, options);
    if (callback) { callback(); }
  });
};

Viewer.prototype.load_map_from_buffer = function load_map_from_buffer (buffer/*:ArrayBuffer*/, options/*:Object*/) {
  var map = new ElMap();
  if (options.format === 'dsn6') {
    map.from_dsn6(buffer);
  } else {
    map.from_ccp4(buffer, true);
  }
  this.add_map(map, options.diff_map);
};

// Load a normal map and a difference map.
// To show the first map ASAP we do not download both maps in parallel.
Viewer.prototype.load_maps = function load_maps (url1/*:string*/, url2/*:string*/,
          options/*:Object*/, callback/*:?Function*/) {
  var format = options.format || 'ccp4';
  var self = this;
  this.load_map(url1, {diff_map: false, format: format}, function () {
    self.load_map(url2, {diff_map: true, format: format}, callback);
  });
};

// Load a model (PDB), normal map and a difference map - in this order.
Viewer.prototype.load_pdb_and_maps = function load_pdb_and_maps (pdb/*:string*/, map1/*:string*/, map2/*:string*/,
                  options/*:Object*/, callback/*:?Function*/) {
  var self = this;
  this.load_pdb(pdb, options, function () {
    self.load_maps(map1, map2, options, callback);
  });
};

// for backward compatibility:
Viewer.prototype.load_ccp4_maps = function load_ccp4_maps (url1/*:string*/, url2/*:string*/, callback/*:?Function*/) {
  this.load_maps(url1, url2, {format: 'ccp4'}, callback);
};
Viewer.prototype.load_pdb_and_ccp4_maps = function load_pdb_and_ccp4_maps (pdb/*:string*/, map1/*:string*/, map2/*:string*/,
                       callback/*:?Function*/) {
  this.load_pdb_and_maps(pdb, map1, map2, {format: 'ccp4'}, callback);
};

// pdb_id here should be lowercase ('1abc')
Viewer.prototype.load_from_pdbe = function load_from_pdbe (pdb_id/*:string*/, callback/*:?Function*/) {
  var id = pdb_id.toLowerCase();
  this.load_pdb_and_maps(
    'https://www.ebi.ac.uk/pdbe/entry-files/pdb' + id + '.ent',
    'https://www.ebi.ac.uk/pdbe/coordinates/files/' + id + '.ccp4',
    'https://www.ebi.ac.uk/pdbe/coordinates/files/' + id + '_diff.ccp4',
    {format: 'ccp4'}, callback);
};
Viewer.prototype.load_from_rcsb = function load_from_rcsb (pdb_id/*:string*/, callback/*:?Function*/) {
  var id = pdb_id.toLowerCase();
  this.load_pdb_and_maps(
    'https://files.rcsb.org/download/' + id + '.pdb',
    'https://edmaps.rcsb.org/maps/' + id + '_2fofc.dsn6',
    'https://edmaps.rcsb.org/maps/' + id + '_fofc.dsn6',
    {format: 'dsn6'}, callback);
};

Viewer.prototype.MOUSE_HELP = [
  '<b>mouse:</b>',
  'Left = rotate',
  'Middle or Ctrl+Left = pan',
  'Right = zoom',
  'Ctrl+Right = clipping',
  'Ctrl+Shift+Right = roll',
  'Wheel = σ level',
  'Shift+Wheel = diff map σ' ].join('\n');

Viewer.prototype.KEYBOARD_HELP = [
  '<b>keyboard:</b>',
  'H = toggle help',
  'S = general style',
  'L = ligand style',
  'T = water style',
  'C = coloring',
  'B = bg color',
  'E = toggle fog',
  'Q = label font',
  '+/- = sigma level',
  ']/[ = map radius',
  'D/F = clip width',
  '&lt;/> = move clip',
  'M/N = zoom',
  'U = unitcell box',
  'Y = hydrogens',
  'V = inactive models',
  'R = center view',
  'W = wireframe style',
  'I = spin',
  'K = rock',
  'Home/End = bond width',
  '\\ = bond caps',
  'P = nearest Cα',
  'Shift+P = permalink',
  '(Shift+)space = next res.',
  'Shift+F = full screen' ].join('\n');

Viewer.prototype.ABOUT_HELP =
  '&nbsp; <a href="https://uglymol.github.io">uglymol</a> ' +
  // $FlowFixMe: Cannot resolve name VERSION.
  (typeof VERSION === 'string' ? VERSION : 'dev'); // eslint-disable-line

Viewer.prototype.ColorSchemes = ColorSchemes;

// @flow


// options handled by Viewer#select_next()
var SPOT_SEL = ['all', 'unindexed', '#1']; //extended when needed
var SHOW_AXES = ['two', 'three', 'none'];
var SPOT_SHAPES = ['wheel', 'square'];

// Modified ElMap for handling output of dials.rs_mapper.
// rs_mapper outputs map in ccp4 format, but we need to rescale it,
// shift it so the box is centered at 0,0,0,
// and the translational symmetry doesn't apply.
var ReciprocalSpaceMap = /*@__PURE__*/(function (ElMap) {
  function ReciprocalSpaceMap(buf /*:ArrayBuffer*/) {
    ElMap.call(this);
    this.box_size = [1, 1, 1];
    this.from_ccp4(buf, false);
    if (this.unit_cell == null) { return; }
    // unit of the map from dials.rs_mapper is (100A)^-1, we scale it to A^-1
    // We assume the "unit cell" is cubic -- as it is in rs_mapper.
    var par = this.unit_cell.parameters;
    this.box_size = [par[0]/ 100, par[1] / 100, par[2] / 100];
    this.unit_cell = null;
  }

  if ( ElMap ) ReciprocalSpaceMap.__proto__ = ElMap;
  ReciprocalSpaceMap.prototype = Object.create( ElMap && ElMap.prototype );
  ReciprocalSpaceMap.prototype.constructor = ReciprocalSpaceMap;

  ReciprocalSpaceMap.prototype.extract_block = function extract_block (radius/*:number*/, center/*:[number,number,number]*/) {
    var grid = this.grid;
    if (grid == null) { return; }
    var b = this.box_size;
    var lo_bounds = [];
    var hi_bounds = [];
    for (var n = 0; n < 3; n++) {
      var lo = Math.floor(grid.dim[n] * ((center[n] - radius) / b[n] + 0.5));
      var hi = Math.floor(grid.dim[n] * ((center[n] + radius) / b[n] + 0.5));
      lo = Math.min(Math.max(0, lo), grid.dim[n] - 1);
      hi = Math.min(Math.max(0, hi), grid.dim[n] - 1);
      if (lo === hi) { return; }
      lo_bounds.push(lo);
      hi_bounds.push(hi);
    }

    var points = [];
    var values = [];
    for (var i = lo_bounds[0]; i <= hi_bounds[0]; i++) {
      for (var j = lo_bounds[1]; j <= hi_bounds[1]; j++) {
        for (var k = lo_bounds[2]; k <= hi_bounds[2]; k++) {
          points.push([(i / grid.dim[0] - 0.5) * b[0],
                       (j / grid.dim[1] - 0.5) * b[1],
                       (k / grid.dim[2] - 0.5) * b[2]]);
          var index = grid.grid2index_unchecked(i, j, k);
          values.push(grid.values[index]);
        }
      }
    }

    var size = [hi_bounds[0] - lo_bounds[0] + 1,
                  hi_bounds[1] - lo_bounds[1] + 1,
                  hi_bounds[2] - lo_bounds[2] + 1];
    this.block.set(points, values, size);
  };

  return ReciprocalSpaceMap;
}(ElMap));

ReciprocalSpaceMap.prototype.unit = '';

function find_max_dist(pos) {
  var max_sq = 0;
  for (var i = 0; i < pos.length; i += 3) {
    var n = 3 * i;
    var sq = pos[n]*pos[n] + pos[n+1]*pos[n+1] + pos[n+2]*pos[n+2];
    if (sq > max_sq) { max_sq = sq; }
  }
  return Math.sqrt(max_sq);
}

function max_val(arr) {
  var max = -Infinity;
  for (var i = 0; i < arr.length; i++) {
    if (arr[i] > max) { max = arr[i]; }
  }
  return max;
}

function parse_csv(text) {
  var lines = text.split('\n').filter(function (line) {
    return line.length > 0 && line[0] !== '#';
  });
  var pos = new Float32Array(lines.length * 3);
  var lattice_ids = [];
  for (var i = 0; i < lines.length; i++) {
    var nums = lines[i].split(',').map(Number);
    for (var j = 0; j < 3; j++) {
      pos[3*i+j] = nums[j];
    }
    lattice_ids.push(nums[3]);
  }
  return { pos: pos, lattice_ids: lattice_ids };
}

function minus_ones(n) {
  var a = [];
  for (var i = 0; i < n; i++) { a.push(-1); }
  return a;
}

function parse_json(text) {
  var d = JSON.parse(text);
  var n = d.rlp.length;
  var pos;
  if (n > 0 && d.rlp[0] instanceof Array) { // deprecated format
    pos = new Float32Array(3*n);
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < 3; j++) {
        pos[3*i+j] = d.rlp[i][j];
      }
    }
  } else { // flat array - new format
    pos = new Float32Array(d.rlp);
  }
  var lattice_ids = d.experiment_id || minus_ones(n);
  return { pos: pos, lattice_ids: lattice_ids };
}

var point_vert = "\nattribute vec3 color;\nattribute float group;\nuniform float show_only;\nuniform float r2_max;\nuniform float r2_min;\nuniform float size;\nvarying vec3 vcolor;\nvoid main() {\n  vcolor = color;\n  float r2 = dot(position, position);\n  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);\n  if (r2 < r2_min || r2 >= r2_max || (show_only != -2.0 && show_only != group))\n    gl_Position.x = 2.0;\n  gl_PointSize = size;\n}";

var round_point_frag = "\n" + fog_pars_fragment + "\nvarying vec3 vcolor;\nvoid main() {\n  // not sure how reliable is such rounding of points\n  vec2 diff = gl_PointCoord - vec2(0.5, 0.5);\n  float dist_sq = 4.0 * dot(diff, diff);\n  if (dist_sq >= 1.0) discard;\n  float alpha = 1.0 - dist_sq * dist_sq * dist_sq;\n  gl_FragColor = vec4(vcolor, alpha);\n" + fog_end_fragment + "\n}";

var square_point_frag = "\n" + fog_pars_fragment + "\nvarying vec3 vcolor;\nvoid main() {\n  gl_FragColor = vec4(vcolor, 1.0);\n" + fog_end_fragment + "\n}";


var ReciprocalViewer = /*@__PURE__*/(function (Viewer) {
  function ReciprocalViewer(options/*:{[key:string]: any}*/) {
    Viewer.call(this, options);
    this.default_camera_pos = [100, 0, 0];
    this.axes = null;
    this.points = null;
    this.max_dist = -1;
    this.d_min = -1;
    this.d_max_inv = 0;
    this.data = {};
    this.config.show_only = SPOT_SEL[0];
    this.config.show_axes = SHOW_AXES[0];
    this.config.spot_shape = SPOT_SHAPES[0];
    this.config.center_cube_size = 0.001;
    this.set_reciprocal_key_bindings();
    if (typeof document !== 'undefined') {
      this.set_dropzone(this.renderer.domElement,
                        this.file_drop_callback.bind(this));
    }
    this.point_material = new ShaderMaterial({
      uniforms: makeUniforms({
        size: 3,
        show_only: -2,
        r2_max: 100,
        r2_min: 0,
      }),
      vertexShader: point_vert,
      fragmentShader: round_point_frag,
      fog: true,
      transparent: true,
      type: 'um_point',
    });
  }

  if ( Viewer ) ReciprocalViewer.__proto__ = Viewer;
  ReciprocalViewer.prototype = Object.create( Viewer && Viewer.prototype );
  ReciprocalViewer.prototype.constructor = ReciprocalViewer;

  ReciprocalViewer.prototype.set_reciprocal_key_bindings = function set_reciprocal_key_bindings () {
    var kb = this.key_bindings;
    // a
    kb[65] = function (evt) {
      this.select_next('axes', 'show_axes', SHOW_AXES, evt.shiftKey);
      this.set_axes();
    };
    // d
    kb[68] = function () { this.change_slab_width_by(-0.01); };
    // f
    kb[70] = function (evt) {
      evt.shiftKey ? this.toggle_full_screen()
                   : this.change_slab_width_by(0.01);
    };
    // p
    kb[80] = function (evt) { this.permalink(); };
    // s
    kb[83] = function (evt) {
      this.select_next('spot shape', 'spot_shape', SPOT_SHAPES, evt.shiftKey);
      if (this.config.spot_shape === 'wheel') {
        this.point_material.fragmentShader = round_point_frag;
      } else {
        this.point_material.fragmentShader = square_point_frag;
      }
      this.point_material.needsUpdate = true;
    };
    // u
    kb[85] = function () {
      if (this.map_bags.length === 0) {
        this.hud('Reciprocal-space density map not loaded.');
        return;
      }
      this.hud('toggled map box');
      this.toggle_cell_box();
    };
    // v
    kb[86] = function (evt) {
      this.select_next('show', 'show_only', SPOT_SEL, evt.shiftKey);
      var idx = SPOT_SEL.indexOf(this.config.show_only);
      this.point_material.uniforms.show_only.value = idx - 2;
    };
    // x
    kb[88] = function (evt) {
      evt.shiftKey ? this.change_map_line(0.1) : this.change_point_size(0.5);
    };
    // z
    kb[90] = function (evt) {
      evt.shiftKey ? this.change_map_line(-0.1) : this.change_point_size(-0.5);
    };
    // comma
    kb[188] = function (evt) { if (evt.shiftKey) { this.shift_clip(0.1); } };
    // period
    kb[190] = function (evt) { if (evt.shiftKey) { this.shift_clip(-0.1); } };
    // <-
    kb[37] = function () { this.change_dmin(0.05); };
    // ->
    kb[39] = function () { this.change_dmin(-0.05); };
    // up arrow
    kb[38] = function () { this.change_dmax(0.025); };
    // down arrow
    kb[40] = function () { this.change_dmax(-0.025); };
    // add, equals/firefox, equal sign
    kb[107] = kb[61] = kb[187] = function (evt) {
      this.change_isolevel_by(0, 0.01);
    };
    // subtract, minus/firefox, dash
    kb[109] = kb[173] = kb[189] = function (evt) {
      this.change_isolevel_by(0, -0.01);
    };
    // [
    kb[219] = function () { this.change_map_radius(-0.001); };
    // ]
    kb[221] = function () { this.change_map_radius(0.001); };
  };

  ReciprocalViewer.prototype.file_drop_callback = function file_drop_callback (file/*:File*/) {
    var self = this;
    var reader = new FileReader();
    if (/\.(map|ccp4)$/.test(file.name)) {
      reader.onloadend = function (evt/*:any*/) {
        if (evt.target.readyState == 2) {
          self.load_map_from_ab(evt.target.result);
        }
      };
      reader.readAsArrayBuffer(file);
    } else {
      reader.onload = function (evt/*:any*/) {
        self.load_from_string(evt.target.result, {});
      };
      reader.readAsText(file);
    }
  };

  ReciprocalViewer.prototype.load_data = function load_data (url/*:string*/, options) {
    if ( options === void 0 ) options/*:Object*/ = {};

    var self = this;
    this.load_file(url, {binary: false, progress: true}, function (req) {
      var ok = self.load_from_string(req.responseText, options);
      if (ok && options.callback) { options.callback(); }
    });
  };

  ReciprocalViewer.prototype.load_from_string = function load_from_string (text/*:string*/, options/*:Object*/) {
    if (text[0] === '{') {
      this.data = parse_json(text);
    } else if (text[0] === '#') {
      this.data = parse_csv(text);
    } else {
      this.hud('Unrecognized file type.');
      return false;
    }
    this.max_dist = find_max_dist(this.data.pos);
    this.d_min = 1 / this.max_dist;
    var last_group = max_val(this.data.lattice_ids);
    SPOT_SEL.splice(3);
    for (var i = 1; i <= last_group; i++) {
      SPOT_SEL.push('#' + (i + 1));
    }
    this.set_axes();
    this.set_points(this.data);
    this.camera.zoom = 0.5 * (this.camera.top - this.camera.bottom);
    // default scale is set to 100 - same as default_camera_pos
    var d = 1.01 * this.max_dist;
    this.controls.slab_width = [d, d, 100];
    this.set_view(options);
    this.hud('Loaded ' + this.data.pos.length + ' spots.');
    return true;
  };

  ReciprocalViewer.prototype.load_map_from_ab = function load_map_from_ab (buffer/*:ArrayBuffer*/) {
    if (this.map_bags.length > 0) {
      this.clear_el_objects(this.map_bags.pop());
    }
    var map = new ReciprocalSpaceMap(buffer);
    if (map == null) { return; }
    var map_range = map.box_size[0] / 2;
    this.config.map_radius = Math.round(map_range / 2 * 100) / 100;
    this.config.max_map_radius = Math.round(1.5 * map_range * 100) / 100;
    this.config.default_isolevel = 2.0;
    this.add_map(map, false);
    var map_dmin = 1 / map_range;
    var msg = 'Loaded density map (' + map_dmin.toFixed(2) + 'Å).\n';
    if (this.points !== null && map_dmin > this.d_min) {
      msg += 'Adjusted spot clipping. ';
      this.change_dmin(map_dmin - this.d_min);
    }
    this.hud(msg + 'Use +/- to change the isolevel.');
  };

  ReciprocalViewer.prototype.set_axes = function set_axes () {
    if (this.axes != null) {
      this.remove_and_dispose(this.axes);
      this.axes = null;
    }
    if (this.config.show_axes === 'none') { return; }
    var axis_length = 1.2 * this.max_dist;
    var vertices = [];
    addXyzCross(vertices, [0, 0, 0], axis_length);
    var ca = this.config.colors.axes;
    var colors = [ca[0], ca[0], ca[1], ca[1], ca[2], ca[2]];
    if (this.config.show_axes === 'two') {
      vertices.splice(4);
      colors.splice(4);
    }
    var material = makeLineMaterial({
      win_size: this.window_size,
      linewidth: 3,
      segments: true,
    });
    this.axes = makeLineSegments(material, vertices, colors);
    this.scene.add(this.axes);
  };

  ReciprocalViewer.prototype.set_points = function set_points (data/*:Object*/) {
    if (this.points != null) {
      this.remove_and_dispose(this.points);
      this.points = null;
    }
    if (data == null || data.lattice_ids == null || data.pos == null) { return; }
    var color_arr = new Float32Array(3 * data.lattice_ids.length);
    this.colorize_by_id(color_arr, data.lattice_ids);
    var geometry = new BufferGeometry();
    geometry.addAttribute('position', new BufferAttribute(data.pos, 3));
    geometry.addAttribute('color', new BufferAttribute(color_arr, 3));
    var groups = new Float32Array(data.lattice_ids);
    geometry.addAttribute('group', new BufferAttribute(groups, 1));
    this.points = new Points(geometry, this.point_material);
    this.scene.add(this.points);
    this.request_render();
  };

  ReciprocalViewer.prototype.colorize_by_id = function colorize_by_id (color_arr/*:Float32Array*/, group_id/*:number[]*/) {
    var palette = this.config.colors.lattices;
    for (var i = 0; i < group_id.length; i++) {
      var c = palette[(group_id[i] + 1) % 4];
      color_arr[3*i] = c.r;
      color_arr[3*i+1] = c.g;
      color_arr[3*i+2] = c.b;
    }
  };

  ReciprocalViewer.prototype.mousewheel_action = function mousewheel_action (delta/*:number*/, evt/*:Event*/) {
    this.change_zoom_by_factor(1 + 0.0005 * delta);
  };

  ReciprocalViewer.prototype.change_point_size = function change_point_size (delta/*:number*/) {
    var size = this.point_material.uniforms.size;
    size.value = Math.max(size.value + delta, 0.5);
    this.hud('point size: ' + size.value.toFixed(1));
  };

  ReciprocalViewer.prototype.change_dmin = function change_dmin (delta/*:number*/) {
    this.d_min = Math.max(this.d_min + delta, 0.1);
    var dmax = this.d_max_inv > 0 ? 1 / this.d_max_inv : null;
    if (dmax !== null && this.d_min > dmax) { this.d_min = dmax; }
    this.point_material.uniforms.r2_max.value = 1 / (this.d_min * this.d_min);
    var low_res = dmax !== null ? dmax.toFixed(2) : '∞';
    this.hud('res. limit: ' + low_res + ' - ' + this.d_min.toFixed(2) + 'Å');
  };

  ReciprocalViewer.prototype.change_dmax = function change_dmax (delta/*:number*/) {
    var v = Math.min(this.d_max_inv + delta, 1 / this.d_min);
    if (v < 1e-6) { v = 0; }
    this.d_max_inv = v;
    this.point_material.uniforms.r2_min.value = v * v;
    var low_res = v > 0 ? (1 / v).toFixed(2) : '∞';
    this.hud('res. limit: ' + low_res + ' - ' + this.d_min.toFixed(2) + 'Å');
  };

  ReciprocalViewer.prototype.redraw_models = function redraw_models () {
    this.set_points(this.data);
  };

  ReciprocalViewer.prototype.get_cell_box_func = function get_cell_box_func () {
    if (this.map_bags.length === 0) { return null; }
    // $FlowFixMe: here the map is ReciprocalSpaceMap not ElMap
    var a = this.map_bags[0].map.box_size;
    return function (xyz/*:[number,number,number]*/) {
      return [(xyz[0]-0.5) * a[0], (xyz[1]-0.5) * a[1], (xyz[2]-0.5) * a[2]];
    };
  };

  return ReciprocalViewer;
}(Viewer));

ReciprocalViewer.prototype.ColorSchemes = [
  {
    name: 'solarized dark',
    bg: 0x002b36,
    fg: 0xfdf6e3,
    map_den: 0xeee8d5,
    center: 0xfdf6e3,
    lattices: [0xdc322f, 0x2aa198, 0x268bd2, 0x859900,
               0xd33682, 0xb58900, 0x6c71c4, 0xcb4b16],
    axes: [0xffaaaa, 0xaaffaa, 0xaaaaff],
  },
  {
    name: 'solarized light',
    bg: 0xfdf6e3,
    fg: 0x002b36,
    map_den: 0x073642,
    center: 0x002b36,
    lattices: [0xdc322f, 0x2aa198, 0x268bd2, 0x859900,
               0xd33682, 0xb58900, 0x6c71c4, 0xcb4b16],
    axes: [0xffaaaa, 0xaaffaa, 0xaaaaff],
  } ];

ReciprocalViewer.prototype.KEYBOARD_HELP = [
  '<b>keyboard:</b>',
  'H = toggle help',
  'V = show (un)indexed',
  'A = toggle axes',
  'U = toggle map box',
  'B = bg color',
  'E = toggle fog',
  'M/N = zoom',
  'D/F = clip width',
  '&lt;/> = move clip',
  'R = center view',
  'Z/X = point size',
  'S = point shape',
  'Shift+P = permalink',
  'Shift+F = full screen',
  '←/→ = max resol.',
  '↑/↓ = min resol.',
  '+/- = map level' ].join('\n');

ReciprocalViewer.prototype.MOUSE_HELP =
    Viewer.prototype.MOUSE_HELP.split('\n').slice(0, -2).join('\n');

// @flow

/*::
 import type {Viewer} from './viewer.js'
 */

function log_timing(t0/*:number*/, text/*:string*/) {
  console.log(text + ': ' + (performance.now() - t0).toFixed(2) + ' ms.');
}

function add_map_from_mtz(viewer, mtz, map_data, is_diff/*:boolean*/) {
  var map = new ElMap();
  var mc = mtz.cell;
  map.unit_cell = new UnitCell(mc.a, mc.b, mc.c, mc.alpha, mc.beta, mc.gamma);
  map.stats.rms = mtz.rmsd;
  map.grid = new GridArray([mtz.nx, mtz.ny, mtz.nz]);
  map.grid.values.set(map_data);
  viewer.add_map(map, is_diff);
}

function load_maps_from_mtz_buffer(viewer/*:Viewer*/, mtz/*:Object*/,
                                   labels/*:?string[]*/) {
  if (labels != null) {
    for (var n = 0; n < labels.length; n += 2) {
      if (labels[n] === '') { continue; }
      var t0 = performance.now();
      var map_data = mtz.calculate_map_from_labels(labels[n], labels[n+1]);
      log_timing(t0, 'map ' + mtz.nx + 'x' + mtz.ny + 'x' + mtz.nz +
                     ' calculated in');
      if (map_data == null) {
        viewer.hud(mtz.last_error, 'ERR');
        continue;
      }
      var is_diff = (n % 4 == 2);
      add_map_from_mtz(viewer, mtz, map_data, is_diff);
    }
  } else {  // use default labels
    for (var nmap = 0; nmap < 2; ++nmap) {
      var is_diff$1 = (nmap == 1);
      var t0$1 = performance.now();
      var map_data$1 = mtz.calculate_map(is_diff$1);
      log_timing(t0$1, 'map ' + mtz.nx + 'x' + mtz.ny + 'x' + mtz.nz +
                     ' calculated in');
      if (map_data$1 != null) {
        add_map_from_mtz(viewer, mtz, map_data$1, is_diff$1);
      }
    }
  }
  mtz.delete();
}

function load_maps_from_mtz(Gemmi/*:Object*/, viewer/*:Viewer*/, url/*:string*/,
                            labels/*:?string[]*/, callback/*:?Function*/) {
  viewer.load_file(url, {binary: true, progress: true}, function (req) {
    var t0 = performance.now();
    try {
      var mtz = Gemmi.readMtz(req.response);
      load_maps_from_mtz_buffer(viewer, mtz, labels);
    } catch (e) {
      viewer.hud(e.message, 'ERR');
      return;
    }
    log_timing(t0, 'load_maps_from_mtz');
    if (callback) { callback(); }
  });
}

function set_pdb_and_mtz_dropzone(Gemmi/*:Object*/, viewer/*:Viewer*/,
                                  zone/*:Object*/) {
  viewer.set_dropzone(zone, function (file) {
    if (/\.mtz$/.test(file.name)) {
      var reader = new FileReader();
      reader.onloadend = function (evt/*:any*/) {
        if (evt.target.readyState == 2) {
          var t0 = performance.now();
          try {
            var mtz = Gemmi.readMtz(evt.target.result);
            load_maps_from_mtz_buffer(viewer, mtz);
          } catch (e) {
            viewer.hud(e.message, 'ERR');
            return;
          }
          log_timing(t0, 'mtz -> maps');
          if (viewer.model_bags.length === 0 && viewer.map_bags.length <= 2) {
            viewer.recenter();
          }
        }
      };
      reader.readAsArrayBuffer(file);
    } else {
      viewer.pick_pdb_and_map(file);
    }
  });
}

exports.AmbientLight = AmbientLight;
exports.Block = Block;
exports.BufferAttribute = BufferAttribute;
exports.BufferGeometry = BufferGeometry;
exports.CatmullRomCurve3 = CatmullRomCurve3;
exports.Color = Color;
exports.Controls = Controls;
exports.ElMap = ElMap;
exports.Fog = Fog;
exports.GridArray = GridArray;
exports.Line = Line;
exports.LineSegments = LineSegments;
exports.Matrix4 = Matrix4;
exports.Mesh = Mesh;
exports.Model = Model;
exports.Object3D = Object3D;
exports.OrthographicCamera = OrthographicCamera;
exports.Points = Points;
exports.Quaternion = Quaternion;
exports.Ray = Ray;
exports.ReciprocalSpaceMap = ReciprocalSpaceMap;
exports.ReciprocalViewer = ReciprocalViewer;
exports.STATE = STATE;
exports.Scene = Scene;
exports.ShaderMaterial = ShaderMaterial;
exports.Texture = Texture;
exports.TriangleStripDrawMode = TriangleStripDrawMode;
exports.UnitCell = UnitCell;
exports.Vector3 = Vector3;
exports.Viewer = Viewer;
exports.WebGLRenderer = WebGLRenderer;
exports.addXyzCross = addXyzCross;
exports.fog_end_fragment = fog_end_fragment;
exports.fog_pars_fragment = fog_pars_fragment;
exports.line_raycast = line_raycast;
exports.load_maps_from_mtz = load_maps_from_mtz;
exports.load_maps_from_mtz_buffer = load_maps_from_mtz_buffer;
exports.makeBalls = makeBalls;
exports.makeChickenWire = makeChickenWire;
exports.makeCube = makeCube;
exports.makeGrid = makeGrid;
exports.makeLabel = makeLabel;
exports.makeLine = makeLine;
exports.makeLineMaterial = makeLineMaterial;
exports.makeLineSegments = makeLineSegments;
exports.makeLines = makeLines;
exports.makeMultiColorLines = makeMultiColorLines;
exports.makeRgbBox = makeRgbBox;
exports.makeRibbon = makeRibbon;
exports.makeSticks = makeSticks;
exports.makeUniforms = makeUniforms;
exports.makeWheels = makeWheels;
exports.modelsFromPDB = modelsFromPDB;
exports.set_pdb_and_mtz_dropzone = set_pdb_and_mtz_dropzone;

Object.defineProperty(exports, '__esModule', { value: true });

})));
