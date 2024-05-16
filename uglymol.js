/*!
 * UglyMol v0.7.2. Macromolecular Viewer for Crystallographers.
 * Copyright 2014 Nat Echols
 * Copyright 2016 Diamond Light Source Ltd
 * Copyright 2016 Marcin Wojdyr
 * Released under the MIT License.
 */
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
typeof define === 'function' && define.amd ? define(['exports'], factory) :
(global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.UM = {}));
})(this, (function (exports) { 'use strict';

var VERSION = exports.VERSION = '0.7.2';


var UnitCell = function UnitCell(a, b, c,
            alpha, beta, gamma) {
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

UnitCell.prototype.fractionalize = function fractionalize (xyz) {
  return multiply(xyz, this.frac);
};

UnitCell.prototype.orthogonalize = function orthogonalize (xyz) {
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

var AMINO_ACIDS = [
  'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
  'LYS', 'MET', 'MSE', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'UNK' ];
var NUCLEIC_ACIDS = [
  'DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U', 'rA', 'rC', 'rG', 'rU',
  'Ar', 'Cr', 'Gr', 'Ur' ];

var NOT_LIGANDS = ['HOH'].concat(AMINO_ACIDS, NUCLEIC_ACIDS);

function modelsFromPDB(pdb_string) {
  var models = [new Model()];
  var pdb_tail = models[0].from_pdb(pdb_string.split('\n'));
  while (pdb_tail != null) {
    var model = new Model();
    pdb_tail = model.from_pdb(pdb_tail);
    if (model.atoms.length > 0) { models.push(model); }
  }
  return models;
}

function modelsFromGemmi(gemmi, buffer, name) {
  var st = gemmi.read_structure(buffer, name);
  var cell = st.cell;  // TODO: check if a copy of cell is created here
  var models = [];
  for (var i_model = 0; i_model < st.length; ++i_model) {
    var model = st.at(i_model);
    var m = new Model();
    m.unit_cell = new UnitCell(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
    var atom_i_seq = 0;
    for (var i_chain = 0; i_chain < model.length; ++i_chain) {
      var chain = model.at(i_chain);
      var chain_name = chain.name;
      for (var i_res = 0; i_res < chain.length; ++i_res) {
        var res = chain.at(i_res);
        var seqid = res.seqid_string;
        var resname = res.name;
        var ent_type = res.entity_type_string;
        var is_ligand = (ent_type === "non-polymer" || ent_type === "branched");
        for (var i_atom = 0; i_atom < res.length; ++i_atom) {
          var atom = res.at(i_atom);
          var new_atom = new Atom();
          new_atom.i_seq = atom_i_seq++;
          new_atom.chain = chain_name;
          new_atom.chain_index = i_chain + 1;
          new_atom.resname = resname;
          new_atom.seqid = seqid;
          new_atom.name = atom.name;
          new_atom.altloc = atom.altloc === 0 ? '' : String.fromCharCode(atom.altloc);
          new_atom.xyz = atom.pos;
          new_atom.occ = atom.occ;
          new_atom.b = atom.b_iso;
          new_atom.element = atom.element_uname;
          new_atom.is_ligand = is_ligand;
          m.atoms.push(new_atom);
        }
      }
    }
    m.calculate_bounds();
    m.calculate_connectivity();
    models.push(m);
  }
  st.delete();
  //console.log("[after modelsFromGemmi] wasm mem:", gemmi.HEAPU8.length / 1024, "kb");
  return models;
}

var Model = function Model() {
  this.atoms = [];
  this.unit_cell = null;
  this.has_hydrogens = false;
  this.lower_bound = [0, 0, 0];
  this.upper_bound = [0, 0, 0];
  this.residue_map = null;
  this.cubes = null;
};

Model.prototype.from_pdb = function from_pdb (pdb_lines) {
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

Model.prototype.next_residue = function next_residue (atom, backward) {
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
Model.prototype.calculate_tangent_vector = function calculate_tangent_vector (residue) {
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

Model.prototype.get_nearest_atom = function get_nearest_atom (x, y, z, atom_name) {
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
  this.name = '';
  this.altloc = '';
  this.resname = '';
  this.chain = '';
  this.chain_index = -1;
  this.seqid = '';
  this.xyz = [0, 0, 0];
  this.occ = 1.0;
  this.b = 0;
  this.element = '';
  this.i_seq = -1;
  this.is_ligand = null;
  this.bonds = [];
};

// http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
Atom.prototype.from_pdb_line = function from_pdb_line (pdb_line) {
  if (pdb_line.length < 66) {
    throw Error('ATOM or HETATM record is too short: ' + pdb_line);
  }
  var rec_type = pdb_line.substring(0, 6);
  if (rec_type !== 'HETATM' && rec_type !== 'ATOM  ') {
    throw Error('Wrong record type: ' + rec_type);
  }
  this.name = pdb_line.substring(12, 16).trim();
  this.altloc = pdb_line.substring(16, 17).trim();
  this.resname = pdb_line.substring(17, 20).trim();
  this.chain = pdb_line.substring(20, 22).trim();
  this.seqid = pdb_line.substring(22, 27).trim();
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

Atom.prototype.distance_sq = function distance_sq (other) {
  var dx = this.xyz[0] - other.xyz[0];
  var dy = this.xyz[1] - other.xyz[1];
  var dz = this.xyz[2] - other.xyz[2];
  return dx*dx + dy*dy + dz*dz;
};

Atom.prototype.distance = function distance (other) {
  return Math.sqrt(this.distance_sq(other));
};

Atom.prototype.midpoint = function midpoint (other) {
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

Atom.prototype.is_same_conformer = function is_same_conformer (other) {
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

Atom.prototype.is_bonded_to = function is_bonded_to (other) {
  var MAX_DIST = 2.2 * 2.2;
  if (!this.is_same_conformer(other)) { return false; }
  var dxyz2 = this.distance_sq(other);
  if (dxyz2 > MAX_DIST) { return false; }
  if (this.element === 'H' && other.element === 'H') { return false; }
  return dxyz2 <= this.bond_radius() * other.bond_radius();
};

Atom.prototype.resid = function resid () {
  return this.seqid + '/' + this.chain;
};

Atom.prototype.long_label = function long_label () {
  var a = this;// eslint-disable-line @typescript-eslint/no-this-alias
  return a.name + ' /' + a.seqid + ' ' + a.resname + '/' + a.chain +
         ' - occ: ' + a.occ.toFixed(2) + ' bf: ' + a.b.toFixed(2) +
         ' ele: ' + a.element + ' pos: (' + a.xyz[0].toFixed(2) + ',' +
         a.xyz[1].toFixed(2) + ',' + a.xyz[2].toFixed(2) + ')';
};

Atom.prototype.short_label = function short_label () {
  var a = this;// eslint-disable-line @typescript-eslint/no-this-alias
  return a.name + ' /' + a.seqid + ' ' + a.resname + '/' + a.chain;
};


// Partition atoms into boxes for quick neighbor searching.
var Cubicles = function Cubicles(atoms, box_length,
            lower_bound, upper_bound) {
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

Cubicles.prototype.find_box_id = function find_box_id (x, y, z) {
  var xstep = Math.floor((x - this.lower_bound[0]) / this.box_length);
  var ystep = Math.floor((y - this.lower_bound[1]) / this.box_length);
  var zstep = Math.floor((z - this.lower_bound[2]) / this.box_length);
  var box_id = (zstep * this.ydim + ystep) * this.xdim + xstep;
  if (box_id < 0 || box_id >= this.boxes.length) { throw Error('Ups!'); }
  return box_id;
};

Cubicles.prototype.get_nearby_atoms = function get_nearby_atoms (box_id) {
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

var Block = function Block() {
  this._points = null;
  this._values = null;
  this._size = [0, 0, 0];
};

Block.prototype.set = function set (points, values, size) {
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

Block.prototype.empty = function empty (){
  return this._values === null;
};

Block.prototype.isosurface = function isosurface (isolevel, method) {
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


function marchingCubes(dims,
                       values,
                       points,
                       isolevel,
                       method) {
  var snap = (method === 'snapped MC');
  var seg_table = (method === 'squarish' ? segTable2 : segTable);
  var vlist = new Array(12);
  var vert_offsets = calculateVertOffsets(dims);
  var vertex_values = new Float32Array(8);
  var vertex_points = new Array(8);
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

function modulo(a, b) {
  var reminder = a % b;
  return reminder >= 0 ? reminder : reminder + b;
}

var GridArray = function GridArray(dim) {
  this.dim = dim; // dimensions of the grid for the entire unit cell
  this.values = new Float32Array(dim[0] * dim[1] * dim[2]);
};

GridArray.prototype.grid2index = function grid2index (i, j, k) {
  i = modulo(i, this.dim[0]);
  j = modulo(j, this.dim[1]);
  k = modulo(k, this.dim[2]);
  return this.dim[2] * (this.dim[1] * i + j) + k;
};

GridArray.prototype.grid2index_unchecked = function grid2index_unchecked (i, j, k) {
  return this.dim[2] * (this.dim[1] * i + j) + k;
};

GridArray.prototype.grid2frac = function grid2frac (i, j, k) {
  return [i / this.dim[0], j / this.dim[1], k / this.dim[2]];
};

// return grid coordinates (rounded down) for the given fractional coordinates
GridArray.prototype.frac2grid = function frac2grid (xyz) {
  // at one point "| 0" here made extract_block() 40% faster on V8 3.14,
  // but I don't see any effect now
  return [Math.floor(xyz[0] * this.dim[0]) | 0,
          Math.floor(xyz[1] * this.dim[1]) | 0,
          Math.floor(xyz[2] * this.dim[2]) | 0];
};

GridArray.prototype.set_grid_value = function set_grid_value (i, j, k, value) {
  var idx = this.grid2index(i, j, k);
  this.values[idx] = value;
};

GridArray.prototype.get_grid_value = function get_grid_value (i, j, k) {
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

ElMap.prototype.abs_level = function abs_level (sigma) {
  return sigma * this.stats.rms + this.stats.mean;
};

// http://www.ccp4.ac.uk/html/maplib.html#description
// eslint-disable-next-line complexity
ElMap.prototype.from_ccp4 = function from_ccp4 (buf, expand_symmetry) {
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
ElMap.prototype.from_dsn6 = function from_dsn6 (buf) {
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
ElMap.prototype.extract_block = function extract_block (radius, center) {
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

ElMap.prototype.isomesh_in_block = function isomesh_in_block (sigma, method) {
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
        var pos = {x: 0, y: 1, z: 2}[m[1]] ;
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

// Copyright 2010-2023 Three.js Authors
// SPDX-License-Identifier: MIT

var _lut = [
  '00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '0a', '0b', '0c',
  '0d', '0e', '0f', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
  '1a', '1b', '1c', '1d', '1e', '1f', '20', '21', '22', '23', '24', '25', '26',
  '27', '28', '29', '2a', '2b', '2c', '2d', '2e', '2f', '30', '31', '32', '33',
  '34', '35', '36', '37', '38', '39', '3a', '3b', '3c', '3d', '3e', '3f', '40',
  '41', '42', '43', '44', '45', '46', '47', '48', '49', '4a', '4b', '4c', '4d',
  '4e', '4f', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '5a',
  '5b', '5c', '5d', '5e', '5f', '60', '61', '62', '63', '64', '65', '66', '67',
  '68', '69', '6a', '6b', '6c', '6d', '6e', '6f', '70', '71', '72', '73', '74',
  '75', '76', '77', '78', '79', '7a', '7b', '7c', '7d', '7e', '7f', '80', '81',
  '82', '83', '84', '85', '86', '87', '88', '89', '8a', '8b', '8c', '8d', '8e',
  '8f', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '9a', '9b',
  '9c', '9d', '9e', '9f', 'a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8',
  'a9', 'aa', 'ab', 'ac', 'ad', 'ae', 'af', 'b0', 'b1', 'b2', 'b3', 'b4', 'b5',
  'b6', 'b7', 'b8', 'b9', 'ba', 'bb', 'bc', 'bd', 'be', 'bf', 'c0', 'c1', 'c2',
  'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'ca', 'cb', 'cc', 'cd', 'ce', 'cf',
  'd0', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'da', 'db', 'dc',
  'dd', 'de', 'df', 'e0', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9',
  'ea', 'eb', 'ec', 'ed', 'ee', 'ef', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6',
  'f7', 'f8', 'f9', 'fa', 'fb', 'fc', 'fd', 'fe', 'ff'
];

// http://stackoverflow.com/questions/105034/how-to-create-a-guid-uuid-in-javascript/21963136#21963136
function generateUUID() {
  var d0 = (Math.random() * 0xffffffff) | 0;
  var d1 = (Math.random() * 0xffffffff) | 0;
  var d2 = (Math.random() * 0xffffffff) | 0;
  var d3 = (Math.random() * 0xffffffff) | 0;
  var uuid =
    _lut[d0 & 0xff] + _lut[(d0 >> 8) & 0xff] +
    _lut[(d0 >> 16) & 0xff] + _lut[(d0 >> 24) & 0xff] + '-' +
    _lut[d1 & 0xff] + _lut[(d1 >> 8) & 0xff] + '-' +
    _lut[((d1 >> 16) & 0x0f) | 0x40] + _lut[(d1 >> 24) & 0xff] + '-' +
    _lut[(d2 & 0x3f) | 0x80] + _lut[(d2 >> 8) & 0xff] + '-' +
    _lut[(d2 >> 16) & 0xff] + _lut[(d2 >> 24) & 0xff] + _lut[d3 & 0xff] +
    _lut[(d3 >> 8) & 0xff] + _lut[(d3 >> 16) & 0xff] + _lut[(d3 >> 24) & 0xff];
  // .toLowerCase() here flattens concatenated strings to save heap memory space.
  return uuid.toLowerCase();
}

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

// compute euclidean modulo of m % n
// https://en.wikipedia.org/wiki/Modulo_operation
function euclideanModulo(n, m) {
  return ((n % m) + m) % m;
}


var Quaternion = function Quaternion(x, y, z, w) {
  if ( x === void 0 ) x = 0;
  if ( y === void 0 ) y = 0;
  if ( z === void 0 ) z = 0;
  if ( w === void 0 ) w = 1;

  this._x = x;
  this._y = y;
  this._z = z;
  this._w = w;
};

var prototypeAccessors$1 = { x: { configurable: true },y: { configurable: true },z: { configurable: true },w: { configurable: true } };

prototypeAccessors$1.x.get = function () {
  return this._x;
};

prototypeAccessors$1.y.get = function () {
  return this._y;
};

prototypeAccessors$1.z.get = function () {
  return this._z;
};

prototypeAccessors$1.w.get = function () {
  return this._w;
};

Quaternion.prototype.setFromAxisAngle = function setFromAxisAngle (axis, angle) {
  // http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
  // assumes axis is normalized
  var halfAngle = angle / 2, s = Math.sin(halfAngle);

  this._x = axis.x * s;
  this._y = axis.y * s;
  this._z = axis.z * s;
  this._w = Math.cos(halfAngle);

  return this;
};

Quaternion.prototype.setFromRotationMatrix = function setFromRotationMatrix (m) {
  // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
  // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

  var te = m.elements,

    m11 = te[0], m12 = te[4], m13 = te[8],
    m21 = te[1], m22 = te[5], m23 = te[9],
    m31 = te[2], m32 = te[6], m33 = te[10],

    trace = m11 + m22 + m33;

  if (trace > 0) {
    var s = 0.5 / Math.sqrt(trace + 1.0);

    this._w = 0.25 / s;
    this._x = (m32 - m23) * s;
    this._y = (m13 - m31) * s;
    this._z = (m21 - m12) * s;
  } else if (m11 > m22 && m11 > m33) {
    var s$1 = 2.0 * Math.sqrt(1.0 + m11 - m22 - m33);

    this._w = (m32 - m23) / s$1;
    this._x = 0.25 * s$1;
    this._y = (m12 + m21) / s$1;
    this._z = (m13 + m31) / s$1;
  } else if (m22 > m33) {
    var s$2 = 2.0 * Math.sqrt(1.0 + m22 - m11 - m33);

    this._w = (m13 - m31) / s$2;
    this._x = (m12 + m21) / s$2;
    this._y = 0.25 * s$2;
    this._z = (m23 + m32) / s$2;
  } else {
    var s$3 = 2.0 * Math.sqrt(1.0 + m33 - m11 - m22);

    this._w = (m21 - m12) / s$3;
    this._x = (m13 + m31) / s$3;
    this._y = (m23 + m32) / s$3;
    this._z = 0.25 * s$3;
  }

  return this;
};

Quaternion.prototype.setFromUnitVectors = function setFromUnitVectors (vFrom, vTo) {
  // assumes direction vectors vFrom and vTo are normalized
  var r = vFrom.dot(vTo) + 1;
  if (r < Number.EPSILON) {
    // vFrom and vTo point in opposite directions
    r = 0;
    if (Math.abs(vFrom.x) > Math.abs(vFrom.z)) {
      this._x = -vFrom.y;
      this._y = vFrom.x;
      this._z = 0;
      this._w = r;
    } else {
      this._x = 0;
      this._y = -vFrom.z;
      this._z = vFrom.y;
      this._w = r;
    }
  } else {
    // crossVectors( vFrom, vTo ); // inlined to avoid cyclic dependency on Vector3
    this._x = vFrom.y * vTo.z - vFrom.z * vTo.y;
    this._y = vFrom.z * vTo.x - vFrom.x * vTo.z;
    this._z = vFrom.x * vTo.y - vFrom.y * vTo.x;
    this._w = r;
  }
  return this.normalize();
};

Quaternion.prototype.length = function length () {
  return Math.sqrt(this._x * this._x + this._y * this._y +
                   this._z * this._z + this._w * this._w);
};

Quaternion.prototype.normalize = function normalize () {
  var l = this.length();
  if (l === 0) {
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
};

Object.defineProperties( Quaternion.prototype, prototypeAccessors$1 );


var Vector3 = function Vector3(x, y, z) {
  if ( x === void 0 ) x = 0;
  if ( y === void 0 ) y = 0;
  if ( z === void 0 ) z = 0;

  Vector3.prototype.isVector3 = true;
  this.x = x;
  this.y = y;
  this.z = z;
};

Vector3.prototype.set = function set (x, y, z) {
  if (z === undefined) { z = this.z; } // sprite.scale.set(x,y)
  this.x = x;
  this.y = y;
  this.z = z;
  return this;
};

Vector3.prototype.clone = function clone () {
  return new this.constructor(this.x, this.y, this.z);
};

Vector3.prototype.copy = function copy (v) {
  this.x = v.x;
  this.y = v.y;
  this.z = v.z;
  return this;
};

Vector3.prototype.add = function add (v) {
  this.x += v.x;
  this.y += v.y;
  this.z += v.z;
  return this;
};

Vector3.prototype.addVectors = function addVectors (a, b) {
  this.x = a.x + b.x;
  this.y = a.y + b.y;
  this.z = a.z + b.z;
  return this;
};

Vector3.prototype.addScaledVector = function addScaledVector (v, s) {
  this.x += v.x * s;
  this.y += v.y * s;
  this.z += v.z * s;
  return this;
};

Vector3.prototype.sub = function sub (v) {
  this.x -= v.x;
  this.y -= v.y;
  this.z -= v.z;
  return this;
};

Vector3.prototype.subVectors = function subVectors (a, b) {
  this.x = a.x - b.x;
  this.y = a.y - b.y;
  this.z = a.z - b.z;
  return this;
};

Vector3.prototype.multiplyScalar = function multiplyScalar (scalar) {
  this.x *= scalar;
  this.y *= scalar;
  this.z *= scalar;
  return this;
};

Vector3.prototype.applyMatrix4 = function applyMatrix4 (m) {
  var x = this.x, y = this.y, z = this.z;
  var e = m.elements;
  var w = 1 / (e[3] * x + e[7] * y + e[11] * z + e[15]);
  this.x = (e[0] * x + e[4] * y + e[8] * z + e[12]) * w;
  this.y = (e[1] * x + e[5] * y + e[9] * z + e[13]) * w;
  this.z = (e[2] * x + e[6] * y + e[10] * z + e[14]) * w;
  return this;
};

Vector3.prototype.applyQuaternion = function applyQuaternion (q) {
  // quaternion q is assumed to have unit length
  var vx = this.x, vy = this.y, vz = this.z;
  var qx = q.x, qy = q.y, qz = q.z, qw = q.w;
  // t = 2 * cross( q.xyz, v );
  var tx = 2 * (qy * vz - qz * vy);
  var ty = 2 * (qz * vx - qx * vz);
  var tz = 2 * (qx * vy - qy * vx);
  // v + q.w * t + cross( q.xyz, t );
  this.x = vx + qw * tx + qy * tz - qz * ty;
  this.y = vy + qw * ty + qz * tx - qx * tz;
  this.z = vz + qw * tz + qx * ty - qy * tx;
  return this;
};

Vector3.prototype.unproject = function unproject (camera) {
  return this.applyMatrix4(camera.projectionMatrixInverse).applyMatrix4(camera.matrixWorld);
};

Vector3.prototype.transformDirection = function transformDirection (m) {
  // input: THREE.Matrix4 affine matrix
  // vector interpreted as a direction

  var x = this.x, y = this.y, z = this.z;
  var e = m.elements;

  this.x = e[0] * x + e[4] * y + e[8] * z;
  this.y = e[1] * x + e[5] * y + e[9] * z;
  this.z = e[2] * x + e[6] * y + e[10] * z;

  return this.normalize();
};

Vector3.prototype.divideScalar = function divideScalar (scalar) {
  return this.multiplyScalar(1 / scalar);
};

Vector3.prototype.dot = function dot (v) {
  return this.x * v.x + this.y * v.y + this.z * v.z;
};

Vector3.prototype.lengthSq = function lengthSq () {
  return this.x * this.x + this.y * this.y + this.z * this.z;
};

Vector3.prototype.length = function length () {
  return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
};

Vector3.prototype.normalize = function normalize () {
  return this.divideScalar(this.length() || 1);
};

Vector3.prototype.setLength = function setLength (length) {
  return this.normalize().multiplyScalar(length);
};

Vector3.prototype.lerp = function lerp (v, alpha) {
  this.x += (v.x - this.x) * alpha;
  this.y += (v.y - this.y) * alpha;
  this.z += (v.z - this.z) * alpha;
  return this;
};

Vector3.prototype.cross = function cross (v) {
  return this.crossVectors(this, v);
};

Vector3.prototype.crossVectors = function crossVectors (a, b) {
  var ax = a.x, ay = a.y, az = a.z;
  var bx = b.x, by = b.y, bz = b.z;
  this.x = ay * bz - az * by;
  this.y = az * bx - ax * bz;
  this.z = ax * by - ay * bx;
  return this;
};

Vector3.prototype.projectOnVector = function projectOnVector (v) {
  var denominator = v.lengthSq();
  if (denominator === 0) { return this.set(0, 0, 0); }
  var scalar = v.dot(this) / denominator;
  return this.copy(v).multiplyScalar(scalar);
};

Vector3.prototype.projectOnPlane = function projectOnPlane (planeNormal) {
  _vector.copy(this).projectOnVector(planeNormal);
  return this.sub(_vector);
};

Vector3.prototype.distanceTo = function distanceTo (v) {
  return Math.sqrt(this.distanceToSquared(v));
};

Vector3.prototype.distanceToSquared = function distanceToSquared (v) {
  var dx = this.x - v.x, dy = this.y - v.y, dz = this.z - v.z;
  return dx * dx + dy * dy + dz * dz;
};

Vector3.prototype.setFromMatrixPosition = function setFromMatrixPosition (m) {
  var e = m.elements;
  this.x = e[12];
  this.y = e[13];
  this.z = e[14];
  return this;
};

Vector3.prototype.setFromMatrixColumn = function setFromMatrixColumn (m, index) {
  return this.fromArray(m.elements, index * 4);
};

Vector3.prototype.equals = function equals (v) {
  return v.x === this.x && v.y === this.y && v.z === this.z;
};

Vector3.prototype.fromArray = function fromArray (array, offset) {
    if ( offset === void 0 ) offset = 0;

  this.x = array[offset];
  this.y = array[offset + 1];
  this.z = array[offset + 2];
  return this;
};
var _vector = /*@__PURE__*/ new Vector3();


var Vector4 = function Vector4(x, y, z, w) {
  if ( x === void 0 ) x = 0;
  if ( y === void 0 ) y = 0;
  if ( z === void 0 ) z = 0;
  if ( w === void 0 ) w = 1;

  this.x = x;
  this.y = y;
  this.z = z;
  this.w = w;
};

Vector4.prototype.set = function set (x, y, z, w) {
  this.x = x;
  this.y = y;
  this.z = z;
  this.w = w;
  return this;
};

Vector4.prototype.copy = function copy (v) {
  this.x = v.x;
  this.y = v.y;
  this.z = v.z;
  this.w = v.w !== undefined ? v.w : 1;
  return this;
};

Vector4.prototype.multiplyScalar = function multiplyScalar (scalar) {
  this.x *= scalar;
  this.y *= scalar;
  this.z *= scalar;
  this.w *= scalar;
  return this;
};

Vector4.prototype.equals = function equals (v) {
  return v.x === this.x && v.y === this.y && v.z === this.z && v.w === this.w;
};


var Matrix4 = function Matrix4() {
  this.elements = [
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1 ];
};

Matrix4.prototype.set = function set (n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44) {
  var te = this.elements;
  te[0] = n11; te[4] = n12; te[8] = n13; te[12] = n14;
  te[1] = n21; te[5] = n22; te[9] = n23; te[13] = n24;
  te[2] = n31; te[6] = n32; te[10] = n33; te[14] = n34;
  te[3] = n41; te[7] = n42; te[11] = n43; te[15] = n44;
  return this;
};

Matrix4.prototype.copy = function copy (m) {
  var te = this.elements;
  var me = m.elements;

  te[0] = me[0];
  te[1] = me[1];
  te[2] = me[2];
  te[3] = me[3];
  te[4] = me[4];
  te[5] = me[5];
  te[6] = me[6];
  te[7] = me[7];
  te[8] = me[8];
  te[9] = me[9];
  te[10] = me[10];
  te[11] = me[11];
  te[12] = me[12];
  te[13] = me[13];
  te[14] = me[14];
  te[15] = me[15];

  return this;
};

Matrix4.prototype.makeRotationFromQuaternion = function makeRotationFromQuaternion (q) {
  return this.compose(_zero, q, _one);
};

Matrix4.prototype.compose = function compose (position, quaternion, scale) {
  var te = this.elements;

  var x = quaternion._x, y = quaternion._y, z = quaternion._z, w = quaternion._w;
  var x2 = x + x, y2 = y + y, z2 = z + z;
  var xx = x * x2, xy = x * y2, xz = x * z2;
  var yy = y * y2, yz = y * z2, zz = z * z2;
  var wx = w * x2, wy = w * y2, wz = w * z2;

  var sx = scale.x, sy = scale.y, sz = scale.z;

  te[0] = (1 - (yy + zz)) * sx;
  te[1] = (xy + wz) * sx;
  te[2] = (xz - wy) * sx;
  te[3] = 0;

  te[4] = (xy - wz) * sy;
  te[5] = (1 - (xx + zz)) * sy;
  te[6] = (yz + wx) * sy;
  te[7] = 0;

  te[8] = (xz + wy) * sz;
  te[9] = (yz - wx) * sz;
  te[10] = (1 - (xx + yy)) * sz;
  te[11] = 0;

  te[12] = position.x;
  te[13] = position.y;
  te[14] = position.z;
  te[15] = 1;

  return this;
};

Matrix4.prototype.lookAt = function lookAt (eye, target, up) {
  var te = this.elements;

  _z.subVectors(eye, target);

  if (_z.lengthSq() === 0) {
    // eye and target are in the same position

    _z.z = 1;
  }

  _z.normalize();
  _x.crossVectors(up, _z);

  if (_x.lengthSq() === 0) {
    // up and z are parallel

    if (Math.abs(up.z) === 1) {
      _z.x += 0.0001;
    } else {
      _z.z += 0.0001;
    }

    _z.normalize();
    _x.crossVectors(up, _z);
  }

  _x.normalize();
  _y.crossVectors(_z, _x);

  te[0] = _x.x; te[4] = _y.x; te[8] = _z.x;
  te[1] = _x.y; te[5] = _y.y; te[9] = _z.y;
  te[2] = _x.z; te[6] = _y.z; te[10] = _z.z;

  return this;
};

Matrix4.prototype.multiplyMatrices = function multiplyMatrices (a, b) {
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
};

Matrix4.prototype.setPosition = function setPosition (x, y, z) {
  var te = this.elements;

  if (x.isVector3) {
    te[12] = x.x;
    te[13] = x.y;
    te[14] = x.z;
  } else {
    te[12] = x;
    te[13] = y;
    te[14] = z;
  }
  return this;
};

Matrix4.prototype.invert = function invert () {
  // based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
  var te = this.elements,
    n11 = te[0], n21 = te[1], n31 = te[2], n41 = te[3],
    n12 = te[4], n22 = te[5], n32 = te[6], n42 = te[7],
    n13 = te[8], n23 = te[9], n33 = te[10], n43 = te[11],
    n14 = te[12], n24 = te[13], n34 = te[14], n44 = te[15],
    t11 = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44,
    t12 = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44,
    t13 = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44,
    t14 = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;

  var det = n11 * t11 + n21 * t12 + n31 * t13 + n41 * t14;

  if (det === 0) { return this.set(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); }

  var detInv = 1 / det;

  te[0] = t11 * detInv;
  te[1] = (n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 + n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44) * detInv;
  te[2] = (n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 - n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44) * detInv;
  te[3] = (n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 + n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43) * detInv;

  te[4] = t12 * detInv;
  te[5] = (n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 - n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44) * detInv;
  te[6] = (n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 + n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44) * detInv;
  te[7] = (n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 - n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43) * detInv;

  te[8] = t13 * detInv;
  te[9] = (n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 + n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44) * detInv;
  te[10] = (n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 - n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44) * detInv;
  te[11] = (n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 + n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43) * detInv;

  te[12] = t14 * detInv;
  te[13] = (n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 - n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34) * detInv;
  te[14] = (n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 + n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34) * detInv;
  te[15] = (n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33) * detInv;

  return this;
};

Matrix4.prototype.scale = function scale (v) {
  var te = this.elements;
  var x = v.x, y = v.y, z = v.z;

  te[0] *= x; te[4] *= y; te[8] *= z;
  te[1] *= x; te[5] *= y; te[9] *= z;
  te[2] *= x; te[6] *= y; te[10] *= z;
  te[3] *= x; te[7] *= y; te[11] *= z;

  return this;
};

Matrix4.prototype.getMaxScaleOnAxis = function getMaxScaleOnAxis () {
  var te = this.elements;

  var scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
  var scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
  var scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];

  return Math.sqrt(Math.max(scaleXSq, scaleYSq, scaleZSq));
};

Matrix4.prototype.makeOrthographic = function makeOrthographic (left, right, top, bottom, near, far) {
  var te = this.elements;
  var w = 1.0 / (right - left);
  var h = 1.0 / (top - bottom);
  var p = 1.0 / (far - near);

  var x = (right + left) * w;
  var y = (top + bottom) * h;
  var z = (far + near) * p;

  te[0] = 2 * w; te[4] = 0; te[8] = 0; te[12] = -x;
  te[1] = 0; te[5] = 2 * h; te[9] = 0; te[13] = -y;
  te[2] = 0; te[6] = 0; te[10] = - 2 * p; te[14] = -z;
  te[3] = 0; te[7] = 0; te[11] = 0; te[15] = 1;

  return this;
};

var _zero = /*@__PURE__*/ new Vector3(0, 0, 0);
var _one = /*@__PURE__*/ new Vector3(1, 1, 1);
var _x = /*@__PURE__*/ new Vector3();
var _y = /*@__PURE__*/ new Vector3();
var _z = /*@__PURE__*/ new Vector3();


function hue2rgb(p, q, t) {
  if (t < 0) { t += 1; }
  if (t > 1) { t -= 1; }
  if (t < 1 / 6) { return p + (q - p) * 6 * t; }
  if (t < 1 / 2) { return q; }
  if (t < 2 / 3) { return p + (q - p) * 6 * (2 / 3 - t); }
  return p;
}

var Color = function Color(r, g, b) {
  this.isColor = true;
  this.r = 1;
  this.g = 1;
  this.b = 1;
  return this.set(r, g, b);
};

Color.prototype.set = function set (r, g, b) {
  if (g === undefined && b === undefined) {
    // r is THREE.Color, hex or string
    var value = r;
    if (value && value.isColor) {
      this.copy(value);
    } else if (typeof value === 'number') {
      this.setHex(value);
    }
  } else {
    this.setRGB(r, g, b);
  }
  return this;
};

Color.prototype.setHex = function setHex (hex) {
  hex = Math.floor(hex);
  this.r = ((hex >> 16) & 255) / 255;
  this.g = ((hex >> 8) & 255) / 255;
  this.b = (hex & 255) / 255;
  return this;
};

Color.prototype.setRGB = function setRGB (r, g, b) {
  this.r = r;
  this.g = g;
  this.b = b;
  return this;
};

Color.prototype.setHSL = function setHSL (h, s, l) {
  // h,s,l ranges are in 0.0 - 1.0
  h = euclideanModulo(h, 1);
  s = clamp(s, 0, 1);
  l = clamp(l, 0, 1);

  if (s === 0) {
    this.r = this.g = this.b = l;
  } else {
    var p = l <= 0.5 ? l * (1 + s) : l + s - l * s;
    var q = 2 * l - p;

    this.r = hue2rgb(q, p, h + 1 / 3);
    this.g = hue2rgb(q, p, h);
    this.b = hue2rgb(q, p, h - 1 / 3);
  }
  return this;
};

Color.prototype.getHSL = function getHSL (target) {
  // h,s,l ranges are in 0.0 - 1.0

  var r = this.r, g = this.g, b = this.b;

  var max = Math.max(r, g, b);
  var min = Math.min(r, g, b);

  var hue, saturation;
  var lightness = (min + max) / 2.0;

  if (min === max) {
    hue = 0;
    saturation = 0;
  } else {
    var delta = max - min;

    saturation = lightness <= 0.5 ? delta / (max + min) : delta / (2 - max - min);

    switch (max) {
      case r:
        hue = (g - b) / delta + (g < b ? 6 : 0);
        break;
      case g:
        hue = (b - r) / delta + 2;
        break;
      case b:
        hue = (r - g) / delta + 4;
        break;
    }

    hue /= 6;
  }
  target.h = hue;
  target.s = saturation;
  target.l = lightness;

  return target;
};

Color.prototype.clone = function clone () {
  return new this.constructor(this.r, this.g, this.b);
};

Color.prototype.copy = function copy (color) {
  this.r = color.r;
  this.g = color.g;
  this.b = color.b;

  return this;
};

Color.prototype.getHex = function getHex () {
  return ( this.r * 255 ) << 16 ^ ( this.g * 255 ) << 8 ^ ( this.b * 255 ) << 0;
};

Color.prototype.getHexString = function getHexString () {
  return ('000000' + this.getHex().toString(16)).slice(-6);
};


//const _vector is already defined above
var _segCenter = /*@__PURE__*/ new Vector3();
var _segDir = /*@__PURE__*/ new Vector3();
var _diff = /*@__PURE__*/ new Vector3();

var Ray = function Ray(origin, direction) {
  if ( origin === void 0 ) origin = new Vector3();
  if ( direction === void 0 ) direction = new Vector3(0, 0, -1);

  this.origin = origin;
  this.direction = direction;
};

Ray.prototype.copy = function copy (ray) {
  this.origin.copy(ray.origin);
  this.direction.copy(ray.direction);
  return this;
};

Ray.prototype.distanceSqToPoint = function distanceSqToPoint (point) {
  var directionDistance = _vector.subVectors(point, this.origin).dot(this.direction);
  // point behind the ray
  if (directionDistance < 0) {
    return this.origin.distanceToSquared(point);
  }
  _vector.copy(this.origin).addScaledVector(this.direction, directionDistance);
  return _vector.distanceToSquared(point);
};

Ray.prototype.distanceSqToSegment = function distanceSqToSegment (v0, v1, optionalPointOnRay, optionalPointOnSegment) {
  // from https://github.com/pmjoniak/GeometricTools/blob/master/GTEngine/Include/Mathematics/GteDistRaySegment.h
  // It returns the min distance between the ray and the segment
  // defined by v0 and v1
  // It can also set two optional targets :
  // - The closest point on the ray
  // - The closest point on the segment

  _segCenter.copy(v0).add(v1).multiplyScalar(0.5);
  _segDir.copy(v1).sub(v0).normalize();
  _diff.copy(this.origin).sub(_segCenter);

  var segExtent = v0.distanceTo(v1) * 0.5;
  var a01 = -this.direction.dot(_segDir);
  var b0 = _diff.dot(this.direction);
  var b1 = -_diff.dot(_segDir);
  var c = _diff.lengthSq();
  var det = Math.abs(1 - a01 * a01);
  var s0, s1, sqrDist, extDet;

  if (det > 0) {
    // The ray and segment are not parallel.

    s0 = a01 * b1 - b0;
    s1 = a01 * b0 - b1;
    extDet = segExtent * det;

    if (s0 >= 0) {
      if (s1 >= -extDet) {
        if (s1 <= extDet) {
          // region 0
          // Minimum at interior points of ray and segment.

          var invDet = 1 / det;
          s0 *= invDet;
          s1 *= invDet;
          sqrDist = s0 * (s0 + a01 * s1 + 2 * b0) + s1 * (a01 * s0 + s1 + 2 * b1) + c;
        } else {
          // region 1

          s1 = segExtent;
          s0 = Math.max(0, -(a01 * s1 + b0));
          sqrDist = -s0 * s0 + s1 * (s1 + 2 * b1) + c;
        }
      } else {
        // region 5

        s1 = -segExtent;
        s0 = Math.max(0, -(a01 * s1 + b0));
        sqrDist = -s0 * s0 + s1 * (s1 + 2 * b1) + c;
      }
    } else {
      if (s1 <= -extDet) {
        // region 4

        s0 = Math.max(0, -(-a01 * segExtent + b0));
        s1 = s0 > 0 ? -segExtent : Math.min(Math.max(-segExtent, -b1), segExtent);
        sqrDist = -s0 * s0 + s1 * (s1 + 2 * b1) + c;
      } else if (s1 <= extDet) {
        // region 3

        s0 = 0;
        s1 = Math.min(Math.max(-segExtent, -b1), segExtent);
        sqrDist = s1 * (s1 + 2 * b1) + c;
      } else {
        // region 2

        s0 = Math.max(0, -(a01 * segExtent + b0));
        s1 = s0 > 0 ? segExtent : Math.min(Math.max(-segExtent, -b1), segExtent);
        sqrDist = -s0 * s0 + s1 * (s1 + 2 * b1) + c;
      }
    }
  } else {
    // Ray and segment are parallel.

    s1 = a01 > 0 ? -segExtent : segExtent;
    s0 = Math.max(0, -(a01 * s1 + b0));
    sqrDist = -s0 * s0 + s1 * (s1 + 2 * b1) + c;
  }

  if (optionalPointOnRay) {
    optionalPointOnRay.copy(this.origin).addScaledVector(this.direction, s0);
  }

  if (optionalPointOnSegment) {
    optionalPointOnSegment.copy(_segCenter).addScaledVector(_segDir, s1);
  }

  return sqrDist;
};

Ray.prototype.applyMatrix4 = function applyMatrix4 (matrix4) {
  this.origin.applyMatrix4(matrix4);
  this.direction.transformDirection(matrix4);
  return this;
};

// Copyright 2010-2023 Three.js Authors
// SPDX-License-Identifier: MIT


// constants.js
var NoBlending = 0;
var NormalBlending = 1;

// core/EventDispatcher.js
var EventDispatcher = function EventDispatcher () {};

EventDispatcher.prototype.addEventListener = function addEventListener (type, listener) {
  if (this._listeners === undefined) { this._listeners = {}; }

  var listeners = this._listeners;

  if (listeners[type] === undefined) {
    listeners[type] = [];
  }

  if (listeners[type].indexOf(listener) === -1) {
    listeners[type].push(listener);
  }
};

EventDispatcher.prototype.removeEventListener = function removeEventListener (type, listener) {
  if (this._listeners === undefined) { return; }

  var listeners = this._listeners;
  var listenerArray = listeners[type];

  if (listenerArray !== undefined) {
    var index = listenerArray.indexOf(listener);

    if (index !== -1) {
      listenerArray.splice(index, 1);
    }
  }
};

EventDispatcher.prototype.dispatchEvent = function dispatchEvent (event) {
  if (this._listeners === undefined) { return; }

  var listeners = this._listeners;
  var listenerArray = listeners[event.type];

  if (listenerArray !== undefined) {
    event.target = this;

    // Make a copy, in case listeners are removed while iterating.
    var array = listenerArray.slice(0);

    for (var i = 0, l = array.length; i < l; i++) {
      array[i].call(this, event);
    }

    event.target = null;
  }
};

// textures/Source.js
var _sourceId = 0;
var Source = function Source(data) {
  if ( data === void 0 ) data = null;

  Object.defineProperty(this, 'id', { value: _sourceId++ });
  this.uuid = generateUUID();
  this.data = data;
  this.dataReady = true;
  this.version = 0;
};

var prototypeAccessors = { needsUpdate: { configurable: true } };
prototypeAccessors.needsUpdate.set = function ( value ) {
  if (value === true) { this.version++; }
};

Object.defineProperties( Source.prototype, prototypeAccessors );

// textures/Texture.js
var _textureId = 0;
var Texture = /*@__PURE__*/(function (EventDispatcher) {
  function Texture(image) {
    EventDispatcher.call(this);
    Object.defineProperty(this, 'id', { value: _textureId++ });
    this.uuid = generateUUID();
    this.name = '';
    this.source = new Source(image);
    this.version = 0;
  }

  if ( EventDispatcher ) Texture.__proto__ = EventDispatcher;
  Texture.prototype = Object.create( EventDispatcher && EventDispatcher.prototype );
  Texture.prototype.constructor = Texture;

  var prototypeAccessors$1 = { image: { configurable: true },needsUpdate: { configurable: true } };

  prototypeAccessors$1.image.get = function () {
    return this.source.data;
  };

  prototypeAccessors$1.image.set = function (value) {
    this.source.data = value;
  };

  Texture.prototype.dispose = function dispose () {
    this.dispatchEvent({ type: 'dispose' });
  };

  prototypeAccessors$1.needsUpdate.set = function (value) {
    if (value === true) {
      this.version++;
      this.source.needsUpdate = true;
    }
  };

  Object.defineProperties( Texture.prototype, prototypeAccessors$1 );

  return Texture;
}(EventDispatcher));


// renderers/webgl/WebGLUniforms.js
/**
 * Uniforms of a program.
 * Those form a tree structure with a special top-level container for the root,
 * which you get by calling 'new WebGLUniforms( gl, program )'.
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
 * .setValue( gl, value, [textures] )
 *
 *              uploads a uniform value(s)
 *      the 'textures' parameter is needed for sampler uniforms
 *
 *
 * Static methods of the top-level container (textures factorizations):
 *
 * .upload( gl, seq, values, textures )
 *
 *              sets uniforms in 'seq' to 'values[id].value'
 *
 * .seqWithValue( seq, values ) : filteredSeq
 *
 *              filters 'seq' entries with corresponding entry in values
 *
 *
 * Methods of the top-level container (textures factorizations):
 *
 * .setValue( gl, name, value, textures )
 *
 *              sets uniform with  name 'name' to 'value'
 *
 * .setOptional( gl, obj, prop )
 *
 *              like .set for an optional property of the object
 *
 */
var emptyTexture = /*@__PURE__*/ new Texture();

// --- Utilities ---

// Float32Array caches used for uploading Matrix uniforms

var mat4array = new Float32Array(16);
var mat3array = new Float32Array(9);
var mat2array = new Float32Array(4);

function arraysEqual(a, b) {
  if (a.length !== b.length) { return false; }

  for (var i = 0, l = a.length; i < l; i++) {
    if (a[i] !== b[i]) { return false; }
  }

  return true;
}

function copyArray(a, b) {
  for (var i = 0, l = b.length; i < l; i++) {
    a[i] = b[i];
  }
}

// --- Setters ---

// Note: Defining these methods externally, because they come in a bunch
// and this way their names minify.

// Single scalar
function setValueV1f(gl, v) {
  var cache = this.cache;
  if (cache[0] === v) { return; }
  gl.uniform1f(this.addr, v);
  cache[0] = v;
}

// Single float vector (from flat array or THREE.VectorN)

function setValueV2f(gl, v) {
  var cache = this.cache;
  if (v.x !== undefined) {
    if (cache[0] !== v.x || cache[1] !== v.y) {
      gl.uniform2f(this.addr, v.x, v.y);
      cache[0] = v.x;
      cache[1] = v.y;
    }
  } else {
    if (arraysEqual(cache, v)) { return; }
    gl.uniform2fv(this.addr, v);
    copyArray(cache, v);
  }
}

function setValueV3f(gl, v) {
  var cache = this.cache;

  if (v.x !== undefined) {
    if (cache[0] !== v.x || cache[1] !== v.y || cache[2] !== v.z) {
      gl.uniform3f(this.addr, v.x, v.y, v.z);
      cache[0] = v.x;
      cache[1] = v.y;
      cache[2] = v.z;
    }
  } else if (v.r !== undefined) {
    if (cache[0] !== v.r || cache[1] !== v.g || cache[2] !== v.b) {
      gl.uniform3f(this.addr, v.r, v.g, v.b);
      cache[0] = v.r;
      cache[1] = v.g;
      cache[2] = v.b;
    }
  } else {
    if (arraysEqual(cache, v)) { return; }
    gl.uniform3fv(this.addr, v);
    copyArray(cache, v);
  }
}

function setValueV4f(gl, v) {
  var cache = this.cache;
  if (v.x !== undefined) {
    if (cache[0] !== v.x || cache[1] !== v.y || cache[2] !== v.z || cache[3] !== v.w) {
      gl.uniform4f(this.addr, v.x, v.y, v.z, v.w);
      cache[0] = v.x;
      cache[1] = v.y;
      cache[2] = v.z;
      cache[3] = v.w;
    }
  } else {
    if (arraysEqual(cache, v)) { return; }
    gl.uniform4fv(this.addr, v);
    copyArray(cache, v);
  }
}

// Single matrix (from flat array or THREE.MatrixN)

function setValueM2(gl, v) {
  var cache = this.cache;
  var elements = v.elements;

  if (elements === undefined) {
    if (arraysEqual(cache, v)) { return; }
    gl.uniformMatrix2fv(this.addr, false, v);
    copyArray(cache, v);
  } else {
    if (arraysEqual(cache, elements)) { return; }
    mat2array.set(elements);
    gl.uniformMatrix2fv(this.addr, false, mat2array);
    copyArray(cache, elements);
  }
}

function setValueM3(gl, v) {
  var cache = this.cache;
  var elements = v.elements;
  if (elements === undefined) {
    if (arraysEqual(cache, v)) { return; }
    gl.uniformMatrix3fv(this.addr, false, v);
    copyArray(cache, v);
  } else {
    if (arraysEqual(cache, elements)) { return; }
    mat3array.set(elements);
    gl.uniformMatrix3fv(this.addr, false, mat3array);
    copyArray(cache, elements);
  }
}

function setValueM4(gl, v) {
  var cache = this.cache;
  var elements = v.elements;
  if (elements === undefined) {
    if (arraysEqual(cache, v)) { return; }
    gl.uniformMatrix4fv(this.addr, false, v);
    copyArray(cache, v);
  } else {
    if (arraysEqual(cache, elements)) { return; }
    mat4array.set(elements);
    gl.uniformMatrix4fv(this.addr, false, mat4array);
    copyArray(cache, elements);
  }
}

// Single integer / boolean

function setValueV1i(gl, v) {
  var cache = this.cache;
  if (cache[0] === v) { return; }
  gl.uniform1i(this.addr, v);
  cache[0] = v;
}

// Single integer / boolean vector (from flat array or THREE.VectorN)

function setValueV2i(gl, v) {
  var cache = this.cache;
  if (v.x !== undefined) {
    if (cache[0] !== v.x || cache[1] !== v.y) {
      gl.uniform2i(this.addr, v.x, v.y);
      cache[0] = v.x;
      cache[1] = v.y;
    }
  } else {
    if (arraysEqual(cache, v)) { return; }
    gl.uniform2iv(this.addr, v);
    copyArray(cache, v);
  }
}

function setValueV3i(gl, v) {
  var cache = this.cache;
  if (v.x !== undefined) {
    if (cache[0] !== v.x || cache[1] !== v.y || cache[2] !== v.z) {
      gl.uniform3i(this.addr, v.x, v.y, v.z);
      cache[0] = v.x;
      cache[1] = v.y;
      cache[2] = v.z;
    }
  } else {
    if (arraysEqual(cache, v)) { return; }
    gl.uniform3iv(this.addr, v);
    copyArray(cache, v);
  }
}

function setValueV4i(gl, v) {
  var cache = this.cache;
  if (v.x !== undefined) {
    if (cache[0] !== v.x || cache[1] !== v.y || cache[2] !== v.z || cache[3] !== v.w) {
      gl.uniform4i(this.addr, v.x, v.y, v.z, v.w);
      cache[0] = v.x;
      cache[1] = v.y;
      cache[2] = v.z;
      cache[3] = v.w;
    }
  } else {
    if (arraysEqual(cache, v)) { return; }
    gl.uniform4iv(this.addr, v);
    copyArray(cache, v);
  }
}

// Single unsigned integer
function setValueV1ui(gl, v) {
  var cache = this.cache;

  if (cache[0] === v) { return; }

  gl.uniform1ui(this.addr, v);

  cache[0] = v;
}

// Single unsigned integer vector (from flat array or THREE.VectorN)

function setValueV2ui(gl, v) {
  var cache = this.cache;
  if (v.x !== undefined) {
    if (cache[0] !== v.x || cache[1] !== v.y) {
      gl.uniform2ui(this.addr, v.x, v.y);

      cache[0] = v.x;
      cache[1] = v.y;
    }
  } else {
    if (arraysEqual(cache, v)) { return; }
    gl.uniform2uiv(this.addr, v);
    copyArray(cache, v);
  }
}

function setValueV3ui(gl, v) {
  var cache = this.cache;
  if (v.x !== undefined) {
    if (cache[0] !== v.x || cache[1] !== v.y || cache[2] !== v.z) {
      gl.uniform3ui(this.addr, v.x, v.y, v.z);
      cache[0] = v.x;
      cache[1] = v.y;
      cache[2] = v.z;
    }
  } else {
    if (arraysEqual(cache, v)) { return; }
    gl.uniform3uiv(this.addr, v);
    copyArray(cache, v);
  }
}

function setValueV4ui(gl, v) {
  var cache = this.cache;
  if (v.x !== undefined) {
    if (cache[0] !== v.x || cache[1] !== v.y || cache[2] !== v.z || cache[3] !== v.w) {
      gl.uniform4ui(this.addr, v.x, v.y, v.z, v.w);

      cache[0] = v.x;
      cache[1] = v.y;
      cache[2] = v.z;
      cache[3] = v.w;
    }
  } else {
    if (arraysEqual(cache, v)) { return; }
    gl.uniform4uiv(this.addr, v);
    copyArray(cache, v);
  }
}

// Single texture (2D / Cube)

function setValueT1(gl, v, textures) {
  var cache = this.cache;
  var unit = textures.allocateTextureUnit();

  if (cache[0] !== unit) {
    gl.uniform1i(this.addr, unit);
    cache[0] = unit;
  }
  var emptyTexture2D = emptyTexture;
  textures.setTexture2D(v || emptyTexture2D, unit);
}

// Helper to pick the right setter for the singular case
function getSingularSetter(type) {
  switch (type) {
    case 0x1406: return setValueV1f; // FLOAT
    case 0x8b50: return setValueV2f; // _VEC2
    case 0x8b51: return setValueV3f; // _VEC3
    case 0x8b52: return setValueV4f; // _VEC4

    case 0x8b5a: return setValueM2; // _MAT2
    case 0x8b5b: return setValueM3; // _MAT3
    case 0x8b5c: return setValueM4; // _MAT4

    case 0x1404: case 0x8b56: return setValueV1i; // INT, BOOL
    case 0x8b53: case 0x8b57: return setValueV2i; // _VEC2
    case 0x8b54: case 0x8b58: return setValueV3i; // _VEC3
    case 0x8b55: case 0x8b59: return setValueV4i; // _VEC4

    case 0x1405: return setValueV1ui; // UINT
    case 0x8dc6: return setValueV2ui; // _VEC2
    case 0x8dc7: return setValueV3ui; // _VEC3
    case 0x8dc8: return setValueV4ui; // _VEC4

    case 0x8b5e: // SAMPLER_2D
    case 0x8d66: // SAMPLER_EXTERNAL_OES
    case 0x8dca: // INT_SAMPLER_2D
    case 0x8dd2: // UNSIGNED_INT_SAMPLER_2D
      return setValueT1;
  }
}

// --- Uniform Classes ---

var SingleUniform = function SingleUniform(id, activeInfo, addr) {
  this.id = id;
  this.addr = addr;
  this.cache = [];
  this.type = activeInfo.type;
  this.setValue = getSingularSetter(activeInfo.type);
  // this.path = activeInfo.name; // DEBUG
};

var StructuredUniform = function StructuredUniform(id) {
  this.id = id;
  this.seq = [];
  this.map = {};
};
StructuredUniform.prototype.setValue = function setValue (gl, value, textures) {
  var seq = this.seq;
  for (var i = 0, n = seq.length; i !== n; ++i) {
    var u = seq[i];
    u.setValue(gl, value[u.id], textures);
  }
};


// --- Top-level ---

// Parser - builds up the property tree from the path strings

var RePathPart = /(\w+)(\])?(\[|\.)?/g;

// extracts
//  - the identifier (member name or array index)
//  - followed by an optional right bracket (found when array index)
//  - followed by an optional left bracket or dot (type of subscript)
//
// Note: These portions can be read in a non-overlapping fashion and
// allow straightforward parsing of the hierarchy that WebGL encodes
// in the uniform names.

function addUniform(container, uniformObject) {
  container.seq.push(uniformObject);
  container.map[uniformObject.id] = uniformObject;
}

function parseUniform(activeInfo, addr, container) {
  var path = activeInfo.name,
    pathLength = path.length;
  // reset RegExp object, because of the early exit of a previous run
  RePathPart.lastIndex = 0;
  while (true) {  // eslint-disable-line no-constant-condition
    var match = RePathPart.exec(path),
      matchEnd = RePathPart.lastIndex;
    var id = match[1];
    var idIsIndex = match[2] === ']',
      subscript = match[3];
    if (idIsIndex) { id = id | 0; } // convert to integer
    if (subscript === undefined || (subscript === '[' && matchEnd + 2 === pathLength)) {
      // bare name or "pure" bottom-level array "[0]" suffix
      if (subscript !== undefined) { throw new TypeError('PureArrayUniform?'); }
      addUniform(container, new SingleUniform(id, activeInfo, addr));
      break;
    } else {
      // step into inner node / create it in case it doesn't exist
      var map = container.map;
      var next = map[id];
      if (next === undefined) {
        next = new StructuredUniform(id);
        addUniform(container, next);
      }
      container = next;
    }
  }
}

// Root Container

var WebGLUniforms = function WebGLUniforms(gl, program) {
  this.seq = [];
  this.map = {};
  var n = gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS);

  for (var i = 0; i < n; ++i) {
    var info = gl.getActiveUniform(program, i),
      addr = gl.getUniformLocation(program, info.name);
    parseUniform(info, addr, this);
  }
};

WebGLUniforms.prototype.setValue = function setValue (gl, name, value, textures) {
  var u = this.map[name];
  if (u !== undefined) { u.setValue(gl, value, textures); }
};

WebGLUniforms.prototype.setOptional = function setOptional (gl, object, name) {
  var v = object[name];
  if (v !== undefined) { this.setValue(gl, name, v); }
};

WebGLUniforms.upload = function upload (gl, seq, values, textures) {
  for (var i = 0, n = seq.length; i !== n; ++i) {
    var u = seq[i],
      v = values[u.id];
    if (v.needsUpdate !== false) {
      // note: always updating when .needsUpdate is undefined
      u.setValue(gl, v.value, textures);
    }
  }
};

WebGLUniforms.seqWithValue = function seqWithValue (seq, values) {
  var r = [];
  for (var i = 0, n = seq.length; i !== n; ++i) {
    var u = seq[i];
    if (u.id in values) { r.push(u); }
  }
  return r;
};



// materials/Material.js
var _materialId = 0;

var Material = /*@__PURE__*/(function (EventDispatcher) {
  function Material() {
    EventDispatcher.call(this);
    this.isMaterial = true;
    Object.defineProperty(this, 'id', { value: _materialId++ });
    this.uuid = generateUUID();
    this.name = '';
    this.type = 'Material';

    this.opacity = 1;
    this.transparent = false;

    this.depthTest = true;
    this.depthWrite = true;

    this.precision = null; // override the renderer's default precision for this material

    this.premultipliedAlpha = false;

    this.visible = true;

    //TODO
    //this.version = 0;
    this._needsUpdate = true;
  }

  if ( EventDispatcher ) Material.__proto__ = EventDispatcher;
  Material.prototype = Object.create( EventDispatcher && EventDispatcher.prototype );
  Material.prototype.constructor = Material;

  var prototypeAccessors$2 = { needsUpdate: { configurable: true } };

  Material.prototype.setValues = function setValues (values) {
    if (values === undefined) { return; }

    for (var key in values) {
      var newValue = values[key];
      if (newValue === undefined) {
        console.warn(("THREE.Material: parameter '" + key + "' has value of undefined."));
        continue;
      }
      var currentValue = this[key];
      if (currentValue === undefined) {
        console.warn(("THREE.Material: '" + key + "' is not a property of THREE." + (this.type) + "."));
        continue;
      }
      if (currentValue && currentValue.isColor) {
        currentValue.set(newValue);
      } else if (currentValue && currentValue.isVector3 && newValue && newValue.isVector3) {
        currentValue.copy(newValue);
      } else {
        this[key] = newValue;
      }
    }
  };

  Material.prototype.dispose = function dispose () {
    this.dispatchEvent({ type: 'dispose' });
  };

  //TODO
  //set needsUpdate(value) {
  //  if (value === true) this.version++;
  //}
  //old:
  prototypeAccessors$2.needsUpdate.get = function () {
    return this._needsUpdate;
  };
  prototypeAccessors$2.needsUpdate.set = function (value) {
    if ( value === true ) { this.update(); }
    this._needsUpdate = value;
  };
  Material.prototype.update = function update () {
    this.dispatchEvent( { type: 'update' } );
  };

  Object.defineProperties( Material.prototype, prototypeAccessors$2 );

  return Material;
}(EventDispatcher));

// materials/ShaderMaterial.js
var ShaderMaterial = /*@__PURE__*/(function (Material) {
  function ShaderMaterial(parameters) {
    Material.call(this);
    this.isShaderMaterial = true;
    this.type = 'ShaderMaterial';
    this.uniforms = {};
    this.vertexShader = '';
    this.fragmentShader = '';
    this.linewidth = 1;
    this.fog = false; // set to use scene fog

    this.extensions = {
      fragDepth: false, // set to use fragment depth values
    };

    this.setValues(parameters);
  }

  if ( Material ) ShaderMaterial.__proto__ = Material;
  ShaderMaterial.prototype = Object.create( Material && Material.prototype );
  ShaderMaterial.prototype.constructor = ShaderMaterial;

  return ShaderMaterial;
}(Material));


// core/Object3D.js
var _object3DId = 0;

var _addedEvent = { type: 'added' };
var _removedEvent = { type: 'removed' };

var Object3D = /*@__PURE__*/(function (EventDispatcher) {
  function Object3D() {
    EventDispatcher.call(this);

    this.isObject3D = true;

    Object.defineProperty(this, 'id', { value: _object3DId++ });

    this.uuid = generateUUID();

    this.name = '';
    this.type = 'Object3D';

    this.parent = null;
    this.children = [];

    this.up = Object3D.DEFAULT_UP.clone();

    var position = new Vector3();
    //const rotation = new Euler();
    var quaternion = new Quaternion();
    var scale = new Vector3(1, 1, 1);

    //function onRotationChange() {
    //  quaternion.setFromEuler(rotation, false);
    //}
    //function onQuaternionChange() {
    //  rotation.setFromQuaternion(quaternion, undefined, false);
    //}
    //rotation._onChange(onRotationChange);
    //quaternion._onChange(onQuaternionChange);

    Object.defineProperties(this, {
      position: {
        configurable: true,
        enumerable: true,
        value: position,
      },
      //rotation: {
      //  configurable: true,
      //  enumerable: true,
      //  value: rotation,
      //},
      quaternion: {
        configurable: true,
        enumerable: true,
        value: quaternion,
      },
      scale: {
        configurable: true,
        enumerable: true,
        value: scale,
      },
      modelViewMatrix: {
        value: new Matrix4(),
      },
      //normalMatrix: {
      //  value: new Matrix3(),
      //},
    });

    this.matrix = new Matrix4();
    this.matrixWorld = new Matrix4();

    this.matrixAutoUpdate = Object3D.DEFAULT_MATRIX_AUTO_UPDATE;

    this.matrixWorldAutoUpdate = Object3D.DEFAULT_MATRIX_WORLD_AUTO_UPDATE; // checked by the renderer
    this.matrixWorldNeedsUpdate = false;

    //this.layers = new Layers();
    this.visible = true;

    //this.castShadow = false;
    //this.receiveShadow = false;

    this.frustumCulled = true;
    this.renderOrder = 0;

    //this.animations = [];

    this.userData = {};
  }

  if ( EventDispatcher ) Object3D.__proto__ = EventDispatcher;
  Object3D.prototype = Object.create( EventDispatcher && EventDispatcher.prototype );
  Object3D.prototype.constructor = Object3D;

  Object3D.prototype.add = function add (object) {
    if (arguments.length > 1) {
      for (var i = 0; i < arguments.length; i++) {
        this.add(arguments[i]);
      }
      return this;
    }
    if (object && object.isObject3D) {
      if (object.parent !== null) {
        object.parent.remove(object);
      }
      object.parent = this;
      this.children.push(object);
      object.dispatchEvent(_addedEvent);
    }
    return this;
  };

  Object3D.prototype.remove = function remove (object) {
    if (arguments.length > 1) {
      for (var i = 0; i < arguments.length; i++) {
        this.remove(arguments[i]);
      }
      return this;
    }

    var index = this.children.indexOf(object);
    if (index !== -1) {
      object.parent = null;
      this.children.splice(index, 1);
      object.dispatchEvent(_removedEvent);
    }
    return this;
  };

  Object3D.prototype.updateMatrix = function updateMatrix () {
    this.matrix.compose(this.position, this.quaternion, this.scale);
    this.matrixWorldNeedsUpdate = true;
  };

  Object3D.prototype.updateMatrixWorld = function updateMatrixWorld (force) {
    if (this.matrixAutoUpdate) { this.updateMatrix(); }

    if (this.matrixWorldNeedsUpdate || force) {
      if (this.parent === null) {
        this.matrixWorld.copy(this.matrix);
      } else {
        this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix);
      }
      this.matrixWorldNeedsUpdate = false;
      force = true;
    }

    // update children
    var children = this.children;
    for (var i = 0, l = children.length; i < l; i++) {
      var child = children[i];
      //if (child.matrixWorldAutoUpdate === true || force === true) {
      child.updateMatrixWorld(force);
      //}
    }
  };

  return Object3D;
}(EventDispatcher));

Object3D.DEFAULT_UP = /*@__PURE__*/ new Vector3(0, 1, 0);
Object3D.DEFAULT_MATRIX_AUTO_UPDATE = true;
Object3D.DEFAULT_MATRIX_WORLD_AUTO_UPDATE = true;



// core/BufferAttribute.js
var BufferAttribute = function BufferAttribute(array, itemSize, normalized) {
  if ( normalized === void 0 ) normalized = false;

  if (Array.isArray(array)) {
    throw new TypeError('BufferAttribute: array should be a Typed Array.');
  }
  this.isBufferAttribute = true;

  this.array = array;
  this.itemSize = itemSize;
  this.count = array !== undefined ? array.length / itemSize : 0;
  this.normalized = normalized;

  // FIXME: new variables
  //this.usage = StaticDrawUsage;
  //this._updateRange = { offset: 0, count: -1 };
  //this.updateRanges = [];
  //this.gpuType = FloatType;
  // FIXME: old variables
  this.dynamic = false;
  this.updateRange = { offset: 0, count: - 1 };
  this.uuid = generateUUID();

  this.version = 0;
};

BufferAttribute.prototype.onUploadCallback = function onUploadCallback () {};


// core/BufferGeometry.js
var _id = 0;

var BufferGeometry = /*@__PURE__*/(function (EventDispatcher) {
  function BufferGeometry() {
    EventDispatcher.call(this);
    this.isBufferGeometry = true;
    Object.defineProperty(this, 'id', { value: _id++ });
    this.uuid = generateUUID();
    this.name = '';
    this.type = 'BufferGeometry';
    this.index = null;
    this.attributes = {};
    this.groups = [];
    this.boundingBox = null;
    this.boundingSphere = null;
    this.drawRange = { start: 0, count: Infinity };
  }

  if ( EventDispatcher ) BufferGeometry.__proto__ = EventDispatcher;
  BufferGeometry.prototype = Object.create( EventDispatcher && EventDispatcher.prototype );
  BufferGeometry.prototype.constructor = BufferGeometry;

  BufferGeometry.prototype.getIndex = function getIndex () {
    return this.index;
  };

  BufferGeometry.prototype.setIndex = function setIndex (index) {
    this.index = index;
  };

  BufferGeometry.prototype.setAttribute = function setAttribute (name, attribute) {
    this.attributes[name] = attribute;
    return this;
  };

  BufferGeometry.prototype.dispose = function dispose () {
    this.dispatchEvent({ type: 'dispose' });
  };

  return BufferGeometry;
}(EventDispatcher));


// objects/Mesh.js
var Mesh = /*@__PURE__*/(function (Object3D) {
  function Mesh(geometry, material) {
    Object3D.call(this);

    this.isMesh = true;

    this.type = 'Mesh';

    if (!geometry) { throw new TypeError('Mesh: geometry not set'); }
    this.geometry = geometry;
    this.material = material;
  }

  if ( Object3D ) Mesh.__proto__ = Object3D;
  Mesh.prototype = Object.create( Object3D && Object3D.prototype );
  Mesh.prototype.constructor = Mesh;

  return Mesh;
}(Object3D));


// cameras/Camera.js
var Camera = /*@__PURE__*/(function (Object3D) {
  function Camera() {
    Object3D.call(this);
    this.isCamera = true;
    this.type = 'Camera';
    this.matrixWorldInverse = new Matrix4();
    this.projectionMatrix = new Matrix4();
    this.projectionMatrixInverse = new Matrix4();
    //this.coordinateSystem = WebGLCoordinateSystem;
  }

  if ( Object3D ) Camera.__proto__ = Object3D;
  Camera.prototype = Object.create( Object3D && Object3D.prototype );
  Camera.prototype.constructor = Camera;

  return Camera;
}(Object3D));

// cameras/OrthographicCamera.js
var OrthographicCamera = /*@__PURE__*/(function (Camera) {
  function OrthographicCamera(left, right, top, bottom, near, far) {
    if ( left === void 0 ) left = -1;
    if ( right === void 0 ) right = 1;
    if ( top === void 0 ) top = 1;
    if ( bottom === void 0 ) bottom = -1;
    if ( near === void 0 ) near = 0.1;
    if ( far === void 0 ) far = 2000;

    Camera.call(this);
    this.type = 'OrthographicCamera';
    this.zoom = 1;
    //this.view = null;
    this.left = left;
    this.right = right;
    this.top = top;
    this.bottom = bottom;
    this.near = near;
    this.far = far;
    this.updateProjectionMatrix();
  }

  if ( Camera ) OrthographicCamera.__proto__ = Camera;
  OrthographicCamera.prototype = Object.create( Camera && Camera.prototype );
  OrthographicCamera.prototype.constructor = OrthographicCamera;

  OrthographicCamera.prototype.updateProjectionMatrix = function updateProjectionMatrix () {
    var dx = (this.right - this.left) / (2 * this.zoom);
    var dy = (this.top - this.bottom) / (2 * this.zoom);
    var cx = (this.right + this.left) / 2;
    var cy = (this.top + this.bottom) / 2;

    var left = cx - dx;
    var right = cx + dx;
    var top = cy + dy;
    var bottom = cy - dy;

    this.projectionMatrix.makeOrthographic(left, right, top, bottom, this.near, this.far);
    this.projectionMatrixInverse.copy(this.projectionMatrix).invert();
  };

  return OrthographicCamera;
}(Camera));


// renderers/webgl/WebGLIndexedBufferRenderer.js (not updated)
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


// renderers/webgl/WebGLBufferRenderer.js (not updated)
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


// renderers/webgl/WebGLShader.js (not updated)
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


// renderers/webgl/WebGLProgram.js (not updated)
var programIdCount = 0;

function generateExtensions( extensions, parameters, rendererExtensions ) {
  extensions = extensions || {};

  var chunks = [
  ( extensions.fragDepth ) && rendererExtensions.get( 'EXT_frag_depth' ) ? '#extension GL_EXT_frag_depth : enable' : '' ];

  return chunks.join( '\n' );
}

function fetchAttributeLocations( gl, program ) {
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
      cachedUniforms = new WebGLUniforms( gl, program );
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


function WebGLPrograms( renderer, capabilities ) {
  var programs = [];

  var parameterNames = [
    'precision',
    'fog', 'useFog',
    'premultipliedAlpha' ];

  this.getParameters = function ( material, fog ) {
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


// renderers/webgl/WebGLGeometries.js (not updated)
function WebGLGeometries( gl, properties ) {
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


// renderers/webgl/WebGLObjects.js (not updated)
function WebGLObjects( gl, properties, info ) {
  var geometries = new WebGLGeometries( gl, properties);

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


// renderers/webgl/WebGLTextures.js (not updated)
function WebGLTextures( _gl, extensions, state, properties, capabilities ) {
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

  var textureUnits = 0;

  function resetTextureUnits() {
    textureUnits = 0;
  }

  function allocateTextureUnit() {
    var textureUnit = textureUnits;
    if (textureUnit >= capabilities.maxTextures) {
      console.warn('WebGLTextures: Trying to use ' + textureUnit +
                   ' texture units while this GPU supports only ' + capabilities.maxTextures);
    }
    textureUnits += 1;
    return textureUnit;
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
  this.resetTextureUnits = resetTextureUnits;
  this.allocateTextureUnit = allocateTextureUnit;
}


// renderers/webgl/WebGLProperties.js (not updated)
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

    clear: function () {
      properties = {};
    }
  };
}


// renderers/webgl/WebGLState.js (not updated)
function WebGLState( gl ) {
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

      reset: function () {
        currentColorClear.set( 0, 0, 0, 1 );
      }
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

  var glVersion = gl.getParameter(gl.VERSION);
  var lineWidthAvailable = false;
  if ( glVersion.indexOf( 'WebGL' ) !== - 1 ) {
    var version = parseFloat( /^WebGL\ (\d)/.exec( glVersion )[1] );
    lineWidthAvailable = ( version >= 1.0 );
  }

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

  function reset() {
    for ( var i = 0; i < enabledAttributes.length; i ++ ) {
      if ( enabledAttributes[ i ] === 1 ) {
        gl.disableVertexAttribArray( i );
        enabledAttributes[ i ] = 0;
      }
    }
    capabilities = {};
    currentTextureSlot = null;
    currentBoundTextures = {};
    currentBlending = null;
    colorBuffer.reset();
    depthBuffer.reset();
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
    reset: reset
  };
}


// renderers/webgl/WebGLCapabilities.js (not updated)
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


// renderers/webgl/WebGLExtensions.js (not updated)
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


// renderers/WebGLRenderer.js (not updated)
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
    _currentMaterialId = - 1,
    _currentGeometryProgram = '',
    _currentCamera = null,

    _currentViewport = new Vector4(),

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
  extensions.get( 'OES_element_index_uint' );

  var capabilities = new WebGLCapabilities( _gl, extensions, parameters );

  var state = new WebGLState( _gl );
  var properties = new WebGLProperties();
  var textures = new WebGLTextures( _gl, extensions, state, properties, capabilities );
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
      renderer.setMode( _gl.TRIANGLES );
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

    if ( scene.matrixWorldAutoUpdate === true ) { scene.updateMatrixWorld(); }

    // update camera matrices and frustum

    if ( camera.parent === null ) { camera.updateMatrixWorld(); }

    camera.matrixWorldInverse.copy( camera.matrixWorld ).invert();

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

    this.setRenderTarget( );

    state.buffers.color.setClear( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha, _premultipliedAlpha );

    if ( this.autoClear || forceClear ) {
      this.clear( this.autoClearColor, this.autoClearDepth );
    }

    // opaque pass (front-to-back order)

    state.setBlending( NoBlending );
    renderObjects( opaqueObjects, scene, camera );

    // transparent pass (back-to-front order)

    renderObjects( transparentObjects, scene, camera );

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
          _vector3.applyMatrix4( _projScreenMatrix );
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
      var material = renderItem.material ;
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
    textures.resetTextureUnits();

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
      p_uniforms.setValue(_gl, 'projectionMatrix', camera.projectionMatrix);

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
        _gl, materialProperties.uniformsList, m_uniforms, textures );
    }


    // common matrices

    p_uniforms.setValue(_gl, 'modelViewMatrix', object.modelViewMatrix);
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

  this.setRenderTarget = function ( ) {
    _currentRenderTarget = null;
    _currentViewport.copy( _viewport ).multiplyScalar( _pixelRatio );
    state.viewport( _currentViewport );
  };
}

// scenes/Fog.js
var Fog = function Fog(color, near, far) {
  if ( near === void 0 ) near = 1;
  if ( far === void 0 ) far = 1000;

  this.isFog = true;
  this.name = '';
  this.color = new Color(color);
  this.near = near;
  this.far = far;
};


// scenes/Scene.js
var Scene = /*@__PURE__*/(function (Object3D) {
  function Scene() {
    Object3D.call(this);
    this.type = 'Scene';
    this.fog = null;
  }

  if ( Object3D ) Scene.__proto__ = Object3D;
  Scene.prototype = Object.create( Object3D && Object3D.prototype );
  Scene.prototype.constructor = Scene;

  return Scene;
}(Object3D));


// objects/Line.js
var Line = /*@__PURE__*/(function (Object3D) {
  function Line(geometry, material) {
    Object3D.call(this);
    this.isLine = true;
    this.type = 'Line';
    this.geometry = geometry;
    this.material = material;
  }

  if ( Object3D ) Line.__proto__ = Object3D;
  Line.prototype = Object.create( Object3D && Object3D.prototype );
  Line.prototype.constructor = Line;

  return Line;
}(Object3D));


// objects/LineSegments.js
var LineSegments = /*@__PURE__*/(function (Line) {
  function LineSegments(geometry, material) {
    Line.call(this, geometry, material);
    this.isLineSegments = true;
    this.type = 'LineSegments';
  }

  if ( Line ) LineSegments.__proto__ = Line;
  LineSegments.prototype = Object.create( Line && Line.prototype );
  LineSegments.prototype.constructor = LineSegments;

  return LineSegments;
}(Line));


// objects/Points.js
var Points = /*@__PURE__*/(function (Object3D) {
  function Points(geometry, material) {
    Object3D.call(this);
    this.isPoints = true;
    this.type = 'Points';
    this.geometry = geometry;
    this.material = material;
  }

  if ( Object3D ) Points.__proto__ = Object3D;
  Points.prototype = Object.create( Object3D && Object3D.prototype );
  Points.prototype.constructor = Points;

  return Points;
}(Object3D));

// kept for compatibility with THREE (lights/AmbientLight.js)
var AmbientLight = function AmbientLight(color) {};

// Copyright 2010-2023 Three.js Authors
// SPDX-License-Identifier: MIT


// from extras/core/Curve.js
var Curve = function Curve() {
  this.type = 'Curve';
  this.arcLengthDivisions = 200;
};

Curve.prototype.getPoints = function getPoints (divisions) {
    if ( divisions === void 0 ) divisions = 5;

  var points = [];
  for (var d = 0; d <= divisions; d++) {
    points.push(this.getPoint(d / divisions));
  }
  return points;
};

// from extras/curves/CatmullRomCurve3.js
/**
 * Centripetal CatmullRom Curve - which is useful for avoiding
 * cusps and self-intersections in non-uniform catmull rom curves.
 * http://www.cemyuksel.com/research/catmullrom_param/catmullrom.pdf
 *
 * curve.type accepts centripetal(default), chordal and catmullrom
 * curve.tension is used for catmullrom which defaults to 0.5
 */

/*
Based on an optimized c++ solution in
 - http://stackoverflow.com/questions/9489736/catmull-rom-curve-with-no-cusps-and-no-self-intersections/
 - http://ideone.com/NoEbVM

This CubicPoly class could be used for reusing some variables and calculations,
but for three.js curve use, it could be possible inlined and flatten into a single function call
which can be placed in CurveUtils.
*/

function CubicPoly() {
  var c0 = 0,
    c1 = 0,
    c2 = 0,
    c3 = 0;

  /*
   * Compute coefficients for a cubic polynomial
   *   p(s) = c0 + c1*s + c2*s^2 + c3*s^3
   * such that
   *   p(0) = x0, p(1) = x1
   *  and
   *   p'(0) = t0, p'(1) = t1.
   */
  function init(x0, x1, t0, t1) {
    c0 = x0;
    c1 = t0;
    c2 = -3 * x0 + 3 * x1 - 2 * t0 - t1;
    c3 = 2 * x0 - 2 * x1 + t0 + t1;
  }

  return {
    initCatmullRom: function (x0, x1, x2, x3, tension) {
      init(x1, x2, tension * (x2 - x0), tension * (x3 - x1));
    },

    initNonuniformCatmullRom: function (x0, x1, x2, x3, dt0, dt1, dt2) {
      // compute tangents when parameterized in [t1,t2]
      var t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1;
      var t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2;

      // rescale tangents for parametrization in [0,1]
      t1 *= dt1;
      t2 *= dt1;

      init(x1, x2, t1, t2);
    },

    calc: function (t) {
      var t2 = t * t;
      var t3 = t2 * t;
      return c0 + c1 * t + c2 * t2 + c3 * t3;
    },
  };
}

//

var tmp = /*@__PURE__*/ new Vector3();
var px = /*@__PURE__*/ new CubicPoly();
var py = /*@__PURE__*/ new CubicPoly();
var pz = /*@__PURE__*/ new CubicPoly();

var CatmullRomCurve3 = /*@__PURE__*/(function (Curve) {
  function CatmullRomCurve3(points, closed, curveType, tension) {
    if ( points === void 0 ) points = [];
    if ( closed === void 0 ) closed = false;
    if ( curveType === void 0 ) curveType = 'centripetal';
    if ( tension === void 0 ) tension = 0.5;

    Curve.call(this);

    this.isCatmullRomCurve3 = true;

    this.type = 'CatmullRomCurve3';

    this.points = points;
    this.closed = closed;
    this.curveType = curveType;
    this.tension = tension;
  }

  if ( Curve ) CatmullRomCurve3.__proto__ = Curve;
  CatmullRomCurve3.prototype = Object.create( Curve && Curve.prototype );
  CatmullRomCurve3.prototype.constructor = CatmullRomCurve3;

  CatmullRomCurve3.prototype.getPoint = function getPoint (t, optionalTarget) {
    if ( optionalTarget === void 0 ) optionalTarget = new Vector3();

    var point = optionalTarget;

    var points = this.points;
    var l = points.length;

    var p = (l - (this.closed ? 0 : 1)) * t;
    var intPoint = Math.floor(p);
    var weight = p - intPoint;

    if (this.closed) {
      intPoint += intPoint > 0 ? 0 : (Math.floor(Math.abs(intPoint) / l) + 1) * l;
    } else if (weight === 0 && intPoint === l - 1) {
      intPoint = l - 2;
      weight = 1;
    }

    var p0, p3; // 4 points (p1 & p2 defined below)

    if (this.closed || intPoint > 0) {
      p0 = points[(intPoint - 1) % l];
    } else {
      // extrapolate first point
      tmp.subVectors(points[0], points[1]).add(points[0]);
      p0 = tmp;
    }

    var p1 = points[intPoint % l];
    var p2 = points[(intPoint + 1) % l];

    if (this.closed || intPoint + 2 < l) {
      p3 = points[(intPoint + 2) % l];
    } else {
      // extrapolate last point
      tmp.subVectors(points[l - 1], points[l - 2]).add(points[l - 1]);
      p3 = tmp;
    }

    if (this.curveType === 'centripetal' || this.curveType === 'chordal') {
      // init Centripetal / Chordal Catmull-Rom
      var pow = this.curveType === 'chordal' ? 0.5 : 0.25;
      var dt0 = Math.pow(p0.distanceToSquared(p1), pow);
      var dt1 = Math.pow(p1.distanceToSquared(p2), pow);
      var dt2 = Math.pow(p2.distanceToSquared(p3), pow);

      // safety check for repeated points
      if (dt1 < 1e-4) { dt1 = 1.0; }
      if (dt0 < 1e-4) { dt0 = dt1; }
      if (dt2 < 1e-4) { dt2 = dt1; }

      px.initNonuniformCatmullRom(p0.x, p1.x, p2.x, p3.x, dt0, dt1, dt2);
      py.initNonuniformCatmullRom(p0.y, p1.y, p2.y, p3.y, dt0, dt1, dt2);
      pz.initNonuniformCatmullRom(p0.z, p1.z, p2.z, p3.z, dt0, dt1, dt2);
    } else if (this.curveType === 'catmullrom') {
      px.initCatmullRom(p0.x, p1.x, p2.x, p3.x, this.tension);
      py.initCatmullRom(p0.y, p1.y, p2.y, p3.y, this.tension);
      pz.initCatmullRom(p0.z, p1.z, p2.z, p3.z, this.tension);
    }

    point.set(px.calc(weight), py.calc(weight), pz.calc(weight));

    return point;
  };

  return CatmullRomCurve3;
}(Curve));

var CUBE_EDGES =
  [[0, 0, 0], [1, 0, 0],
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

function makeColorAttribute(colors) {
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

function makeLines(pos, color, linewidth) {
  var material = new ShaderMaterial({
    uniforms: makeUniforms({vcolor: color}),
    vertexShader: unicolor_vert,
    fragmentShader: unicolor_frag,
    fog: true,
    linewidth: linewidth,
    type: 'um_lines',
  });
  var geometry = new BufferGeometry();
  geometry.setAttribute('position', new BufferAttribute(pos, 3));
  return new LineSegments(geometry, material);
}






function makeCube(size, ctr, options) {
  var pos = new Float32Array(CUBE_EDGES.length * 3);
  for (var i = 0; i < CUBE_EDGES.length; i++) {
    var coor = CUBE_EDGES[i];
    pos[3*i+0] = ctr.x + size * (coor[0] - 0.5);
    pos[3*i+1] = ctr.y + size * (coor[1] - 0.5);
    pos[3*i+2] = ctr.z + size * (coor[2] - 0.5);
  }
  return makeLines(pos, options.color, options.linewidth);
}

function makeMultiColorLines(pos,
                                    colors,
                                    linewidth) {
  var material = new ShaderMaterial({
    uniforms: makeUniforms({}),
    vertexShader: varcolor_vert,
    fragmentShader: varcolor_frag,
    fog: true,
    linewidth: linewidth,
    type: 'um_multicolor_lines',
  });
  var geometry = new BufferGeometry();
  geometry.setAttribute('position', new BufferAttribute(pos, 3));
  geometry.setAttribute('color', makeColorAttribute(colors));
  return new LineSegments(geometry, material);
}

// A cube with 3 edges (for x, y, z axes) colored in red, green and blue.
function makeRgbBox(transform_func, color) {
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

function double_pos(pos) {
  var double_pos = [];
  for (var i = 0; i < pos.length; i++) {
    var v = pos[i];
    double_pos.push(v[0], v[1], v[2]);
    double_pos.push(v[0], v[1], v[2]);
  }
  return double_pos;
}

function double_color(color_arr) {
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


var wide_segments_vert = "\nattribute vec3 color;\nattribute vec3 other;\nattribute float side;\nuniform vec2 win_size;\nuniform float linewidth;\nvarying vec3 vcolor;\n\nvoid main() {\n  vcolor = color;\n  mat4 mat = projectionMatrix * modelViewMatrix;\n  vec2 dir = normalize((mat * vec4(position - other, 0.0)).xy);\n  vec2 normal = vec2(-dir.y, dir.x);\n  gl_Position = mat * vec4(position, 1.0);\n  gl_Position.xy += side * linewidth * normal / win_size;\n}";

function interpolate_vertices(segment, smooth) {
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

function makeUniforms(params) {
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
function makeRibbon(vertices,
                           colors,
                           tangents,
                           smoothness) {
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
  geometry.setAttribute('position', new BufferAttribute(pos, 3));
  geometry.setAttribute('color', makeColorAttribute(color_arr));
  var tan = new Float32Array(tang_arr);
  geometry.setAttribute('tan', new BufferAttribute(tan, 3));
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


function makeChickenWire(data,
                         options) {
  var geom = new BufferGeometry();
  var position = new Float32Array(data.vertices);
  geom.setAttribute('position', new BufferAttribute(position, 3));

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
  geom.setAttribute('position', new BufferAttribute(pos_arr, 3));
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
  return obj;
}


function makeLineMaterial(options) {
  var uniforms = makeUniforms({
    linewidth: options.linewidth,
    win_size: options.win_size,
  });
  return new ShaderMaterial({
    uniforms: uniforms,
    vertexShader: wide_segments_vert,
    fragmentShader: varcolor_frag,
    fog: true,
    type: 'um_line',
  });
}

// vertex_arr and color_arr must be of the same length
function makeLineSegments(material,
                                 vertex_arr,
                                 color_arr) {
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
  geometry.setAttribute('position', new BufferAttribute(position, 3));
  geometry.setAttribute('other', new BufferAttribute(other_vert, 3));
  geometry.setAttribute('side', new BufferAttribute(side, 1));
  if (color_arr != null) {
    var color = double_color(color_arr);
    geometry.setAttribute('color', new BufferAttribute(color, 3));
  }
  geometry.setIndex(make_quad_index_buffer(len/2));

  var mesh = new Mesh(geometry, material);
  //mesh.userData.bond_lines = true;
  return mesh;
}

var wheel_vert = "\nattribute vec3 color;\nuniform float size;\nvarying vec3 vcolor;\nvoid main() {\n  vcolor = color;\n  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);\n  gl_PointSize = size;\n}";

// not sure how portable it is
var wheel_frag = "\n" + fog_pars_fragment + "\nvarying vec3 vcolor;\nvoid main() {\n  vec2 diff = gl_PointCoord - vec2(0.5, 0.5);\n  if (dot(diff, diff) >= 0.25) discard;\n  gl_FragColor = vec4(vcolor, 1.0);\n" + fog_end_fragment + "\n}";

function makeWheels(atom_arr, color_arr, size) {
  var geometry = new BufferGeometry();
  var pos = new Float32Array(atom_arr.length * 3);
  for (var i = 0; i < atom_arr.length; i++) {
    var xyz = atom_arr[i].xyz;
    pos[3*i+0] = xyz[0];
    pos[3*i+1] = xyz[1];
    pos[3*i+2] = xyz[2];
  }
  geometry.setAttribute('position', new BufferAttribute(pos, 3));
  geometry.setAttribute('color', makeColorAttribute(color_arr));
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

function makeSticks(vertex_arr, color_arr, radius) {
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
  geometry.setAttribute('position', new BufferAttribute(position, 3));
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
  geometry.setAttribute('axis', new BufferAttribute(axis, 3));
  geometry.setAttribute('corner', new BufferAttribute(corner, 2));
  var color = double_color(color_arr);
  geometry.setAttribute('color', new BufferAttribute(color, 3));
  geometry.setIndex(make_quad_index_buffer(len/2));

  var mesh = new Mesh(geometry, material);
  //mesh.userData.bond_lines = true;
  return mesh;
}

function makeBalls(atom_arr, color_arr, radius) {
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
  geometry.setAttribute('position', new BufferAttribute(pos, 3));

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
  geometry.setAttribute('corner', new BufferAttribute(corner, 2));

  var colors = new Float32Array(N * 4 * 3);
  for (var i$2 = 0; i$2 < N; i$2++) {
    var col = color_arr[i$2];
    for (var j$1 = 0; j$1 < 4; j$1++) {
      colors[3 * (4*i$2 + j$1) + 0] = col.r;
      colors[3 * (4*i$2 + j$1) + 1] = col.g;
      colors[3 * (4*i$2 + j$1) + 2] = col.b;
    }
  }
  geometry.setAttribute('color', new BufferAttribute(colors, 3));

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

/*
interface LineRaycastOptions {
  precision: number;
  ray: Ray;
  near: number;
  far: number;
}
// based on Line.prototype.raycast(), but skipping duplicated points
const inverseMatrix = new Matrix4();
const ray = new Ray();
export
function line_raycast(mesh: Mesh, options: LineRaycastOptions,
                      intersects: object[]) {
  const precisionSq = options.precision * options.precision;
  inverseMatrix.copy(mesh.matrixWorld).invert();
  ray.copy(options.ray).applyMatrix4(inverseMatrix);
  const vStart = new Vector3();
  const vEnd = new Vector3();
  const interSegment = new Vector3();
  const interRay = new Vector3();
  const positions = mesh.geometry.attributes.position.array;
  for (let i = 0, l = positions.length / 6 - 1; i < l; i += 2) {
    vStart.fromArray(positions, 6 * i);
    vEnd.fromArray(positions, 6 * i + 6);
    const distSq = ray.distanceSqToSegment(vStart, vEnd, interRay, interSegment);
    if (distSq > precisionSq) continue;
    interRay.applyMatrix4(mesh.matrixWorld);
    const distance = options.ray.origin.distanceTo(interRay);
    if (distance < options.near || distance > options.far) continue;
    intersects.push({
      distance: distance,
      point: interSegment.clone().applyMatrix4(mesh.matrixWorld),
      index: i,
      object: mesh,
      line_dist: Math.sqrt(distSq), // extra property, not in Three.js
    });
  }
}
*/

var label_vert = "\nattribute vec2 uvs;\nuniform vec2 canvas_size;\nuniform vec2 win_size;\nuniform float z_shift;\nvarying vec2 vUv;\nvoid main() {\n  vUv = uvs;\n  vec2 rel_offset = vec2(0.02, -0.3);\n  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);\n  gl_Position.xy += (uvs + rel_offset) * 2.0 * canvas_size / win_size;\n  gl_Position.z += z_shift * projectionMatrix[2][2];\n}";

var label_frag = "\n" + fog_pars_fragment + "\nvarying vec2 vUv;\nuniform sampler2D map;\nvoid main() {\n  gl_FragColor = texture2D(map, vUv);\n" + fog_end_fragment + "\n}";

var Label = function Label(text, options) {
  this.texture = new Texture();
  var canvas_size = this.redraw(text, options);
  if (canvas_size === undefined) { return; }

  // Rectangle geometry.
  var geometry = new BufferGeometry();
  var pos = options.pos;
  var position = new Float32Array([].concat(pos, pos, pos, pos));
  var uvs = new Float32Array([0, 1, 1, 1, 0, 0, 1, 0]);
  var indices = new Uint16Array([0, 2, 1, 2, 3, 1]);
  geometry.setIndex(new BufferAttribute(indices, 1));
  geometry.setAttribute('position', new BufferAttribute(position, 3));
  geometry.setAttribute('uvs', new BufferAttribute(uvs, 2));

  var material = new ShaderMaterial({
    uniforms: makeUniforms({map: this.texture,
                            canvas_size: canvas_size,
                            win_size: options.win_size,
                            z_shift: options.z_shift}),
    vertexShader: label_vert,
    fragmentShader: label_frag,
    fog: true,
    type: 'um_label',
  });
  material.transparent = true;
  this.mesh = new Mesh(geometry, material);
};

Label.prototype.redraw = function redraw (text, options) {
  if (typeof document === 'undefined') { return; }// for testing on node
  var canvas = document.createElement('canvas');
  // Canvas size should be 2^N.
  canvas.width = 256;// arbitrary limit, to keep it simple
  canvas.height = 16;// font size
  var context = canvas.getContext('2d');
  if (!context) { return null; }
  context.font = (options.font || 'bold 14px') + ' sans-serif';
  //context.fillStyle = 'green';
  //context.fillRect(0, 0, canvas.width, canvas.height);
  context.textBaseline = 'bottom';
  if (options.color) { context.fillStyle = options.color; }
  context.fillText(text, 0, canvas.height);
  this.texture.image = canvas;
  this.texture.needsUpdate = true;
  return [canvas.width, canvas.height];
};


// Add vertices of a 3d cross (representation of an unbonded atom)
function addXyzCross(vertices, xyz, r) {
  vertices.push([xyz[0]-r, xyz[1], xyz[2]], [xyz[0]+r, xyz[1], xyz[2]]);
  vertices.push([xyz[0], xyz[1]-r, xyz[2]], [xyz[0], xyz[1]+r, xyz[2]]);
  vertices.push([xyz[0], xyz[1], xyz[2]-r], [xyz[0], xyz[1], xyz[2]+r]);
}

// Properties defined with Object.defineProperties() in JS are not understood
// by TypeScript; add them here.
 






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

// object used in computations (to avoid creating and deleting it)
var _m1 = new Matrix4();

var STATE = { NONE: -1, ROTATE: 0, PAN: 1, ZOOM: 2, PAN_ZOOM: 3,
                       SLAB: 4, ROLL: 5, AUTO_ROTATE: 6, GO: 7 };

var auto_speed = 1.0;

// based on three.js/examples/js/controls/OrthographicTrackballControls.js
var Controls = function Controls(camera, target) {
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

Controls.prototype._rotate_camera = function _rotate_camera (eye) {
  var quat = new Quaternion();
  quat.setFromUnitVectors(this._rotate_end, this._rotate_start);
  eye.applyQuaternion(quat);
  this._camera.up.applyQuaternion(quat);
  this._rotate_end.applyQuaternion(quat);
  this._rotate_start.copy(this._rotate_end);
};

Controls.prototype._zoom_camera = function _zoom_camera (eye) {
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

Controls.prototype._pan_camera = function _pan_camera (eye) {
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

Controls.prototype._auto_rotate = function _auto_rotate (eye) {
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

Controls.prototype.toggle_auto = function toggle_auto (param) {
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
  //this._camera.lookAt(this._target);
  _m1.lookAt(this._camera.position, this._target, this._camera.up);
  this._camera.quaternion.setFromRotationMatrix(_m1);
  return changed;
};

Controls.prototype.start = function start (new_state, x, y, dist) {
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

Controls.prototype.move = function move (x, y, dist) {
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
Controls.prototype.go_to = function go_to (targ, cam_pos, cam_up, steps) {
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

var ColorSchemes$1 = {
  // the default scheme that generally mimicks Coot
  'coot dark': {
    bg: new Color(0x000000),
    fg: new Color(0xFFFFFF),
    map_den: new Color(0x3362B2),
    map_pos: new Color(0x298029),
    map_neg: new Color(0x8B2E2E),
    center: new Color(0xC997B0),
    // atoms
    H: new Color(0x858585), // H is normally invisible
    // C, N and O are taken approximately (by color-picker) from coot
    C: new Color(0xb3b300),
    N: new Color(0x7EAAFB),
    O: new Color(0xF24984),
    S: new Color(0x40ff40), // S in coot is too similar to C, here it's greener
    // Coot doesn't define other colors (?)
    MG: new Color(0xc0c0c0),
    P:  new Color(0xffc040),
    CL: new Color(0xa0ff60),
    CA: new Color(0xffffff),
    MN: new Color(0xff90c0),
    FE: new Color(0xa03000),
    NI: new Color(0x00ff80),
    def: new Color(0xa0a0a0), // default atom color
  },

  // scheme made of "solarized" colors (http://ethanschoonover.com/solarized):
  // base03  base02  base01  base00  base0   base1   base2   base3
  // #002b36 #073642 #586e75 #657b83 #839496 #93a1a1 #eee8d5 #fdf6e3
  // yellow  orange  red     magenta violet  blue    cyan    green
  // #b58900 #cb4b16 #dc322f #d33682 #6c71c4 #268bd2 #2aa198 #859900
  'solarized dark': {
    bg: new Color(0x002b36),
    fg: new Color(0xfdf6e3),
    map_den: new Color(0x268bd2),
    map_pos: new Color(0x859900),
    map_neg: new Color(0xd33682),
    center: new Color(0xfdf6e3),
    H: new Color(0x586e75),
    C: new Color(0x93a1a1),
    N: new Color(0x6c71c4),
    O: new Color(0xcb4b16),
    S: new Color(0xb58900),
    def: new Color(0xeee8d5),
  },

  'solarized light': {
    bg: new Color(0xfdf6e3),
    fg: new Color(0x002b36),
    map_den: new Color(0x268bd2),
    map_pos: new Color(0x859900),
    map_neg: new Color(0xd33682),
    center: new Color(0x002b36),
    H: new Color(0x93a1a1),
    C: new Color(0x586e75),
    N: new Color(0x6c71c4),
    O: new Color(0xcb4b16),
    S: new Color(0xb58900),
    def: new Color(0x073642),
  },

  // like in Coot after Edit > Background Color > White
  'coot light': {
    bg: new Color(0xFFFFFF),
    fg: new Color(0x000000),
    map_den: new Color(0x3362B2),
    map_pos: new Color(0x298029),
    map_neg: new Color(0x8B2E2E),
    center: new Color(0xC7C769),
    H: new Color(0x999999),
    C: new Color(0xA96464),
    N: new Color(0x1C51B3),
    O: new Color(0xC33869),
    S: new Color(0x9E7B3D),
    def: new Color(0x808080),
  },
};


var INIT_HUD_TEXT = 'This is UglyMol not Coot. ' +
  '<a href="#" onclick="V.toggle_help(); return false;">H shows help.</a>';

// options handled by select_next()

var COLOR_PROPS = ['element', 'B-factor', 'pLDDT', 'occupancy', 'index', 'chain'];
var RENDER_STYLES = ['lines', 'trace', 'ribbon', 'ball&stick'];
var LIGAND_STYLES = ['ball&stick', 'lines'];
var WATER_STYLES = ['cross', 'dot', 'invisible'];
var MAP_STYLES = ['marching cubes', 'squarish' ];
var LINE_STYLES = ['normal', 'simplistic'];
var LABEL_FONTS = ['bold 14px', '14px', '16px', 'bold 16px'];

function rainbow_value(v, vmin, vmax) {
  var c = new Color(0xe0e0e0);
  if (vmin < vmax) {
    var ratio = (v - vmin) / (vmax - vmin);
    var hue = (240 - (240 * ratio)) / 360;
    c.setHSL(hue, 1.0, 0.5);
  }
  return c;
}

function color_by(prop, atoms, elem_colors,
                  hue_shift) {
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
  } else if (prop === 'pLDDT') {
    var steps = [90, 70, 50];
    var colors = [
      new Color(0x0053d6), // dark blue
      new Color(0x65cbf3), // light blue
      new Color(0xffdb13), // yellow
      new Color(0xff7d45)  // orange
    ];
    color_func = function (atom) {
      var i = 0;
      while (i < 3 && atom.b < steps[i]) {
        ++i;
      }
      return colors[i];
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
      var c_hsl = { h: 0, s: 0, l: 0 };
      elem_colors['C'].getHSL(c_hsl);
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

var MapBag = function MapBag(map, config, is_diff_map) {
  this.map = map;
  this.name = '';
  this.isolevel = is_diff_map ? 3.0 : config.default_isolevel;
  this.visible = true;
  this.types = is_diff_map ? ['map_pos', 'map_neg'] : ['map_den'];
  this.block_ctr = new Vector3(Infinity, 0, 0);
  this.el_objects = []; // three.js objects
};

var ModelBag = function ModelBag(model, config, win_size) {
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

ModelBag.prototype.add_bonds = function add_bonds (polymers, ligands, ball_size) {
  var visible_atoms = this.get_visible_atoms();
  var colors = color_by(this.conf.color_prop, visible_atoms,
                          this.conf.colors, this.hue_shift);
  var vertex_arr = [];
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
  var vertex_arr = [];
  var color_arr = [];
  var k = 0;
  for (var i$1 = 0, list = segments; i$1 < list.length; i$1 += 1) {
    var seg = list[i$1];

      for (var i = 1; i < seg.length; ++i) {
      vertex_arr.push(seg[i-1].xyz, seg[i].xyz);
      color_arr.push(colors[k+i-1], colors[k+i]);
    }
    k += seg.length;
  }
  var linewidth = scale_by_height(this.conf.bond_line, this.win_size);
  var material = makeLineMaterial({
    linewidth: linewidth,
    win_size: this.win_size,
  });
  this.objects.push(makeLineSegments(material, vertex_arr, color_arr));
  if (this.conf.line_style !== 'simplistic') {
    // wheels (discs) as round caps
    this.objects.push(makeWheels(visible_atoms, colors, linewidth));
  }
  this.atom_array = visible_atoms;
};

ModelBag.prototype.add_ribbon = function add_ribbon (smoothness) {
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
function touch_info(evt) {
  var touches = evt.touches;
  var dx = touches[0].pageX - touches[1].pageX;
  var dy = touches[0].pageY - touches[1].pageY;
  return {pageX: (touches[0].pageX + touches[1].pageX) / 2,
          pageY: (touches[0].pageY + touches[1].pageY) / 2,
          dist: Math.sqrt(dx * dx + dy * dy)};
}

// makes sense only for full-window viewer
function parse_url_fragment() {
  var ret  = {};
  if (typeof window === 'undefined') { return ret; }
  var params = window.location.hash.substr(1).split('&');
  for (var i = 0; i < params.length; i++) {
    var kv = params[i].split('=');
    var key = kv[0];
    var val = kv[1];
    if (key === 'xyz' || key === 'eye') {
      ret[key] = val.split(',').map(Number);
    } else if (key === 'zoom') {
      ret[key] = Number(val);
    } else {
      ret[key] = val;
    }
  }
  return ret;
}


var Viewer = function Viewer(options) {
  // rendered objects
  this.model_bags = [];
  this.map_bags = [];
  this.decor = {
    cell_box: null,
    selection: null,
    zoom_grid: makeGrid(),
    mark: null,
  };
  this.labels = {};
  //this.nav = null;
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
    color_scheme: 'coot dark',
    // `colors` is assigned in set_colors()
    hydrogens: false,
    ball_size: 0.4,
  };

  // options of the constructor overwrite default values of the config
  for (var i = 0, list = Object.keys(options); i < list.length; i += 1) {
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
  this.scene.add(new AmbientLight());
  this.default_camera_pos = [0, 0, 100];
  if (options.share_view) {
    this.target = options.share_view.target;
    this.camera = options.share_view.camera;
    this.controls = options.share_view.controls;
    this.tied_viewer = options.share_view;
    this.tied_viewer.tied_viewer = this; // not GC friendly
  } else {
    this.target = new Vector3(0, 0, 0);
    this.camera = new OrthographicCamera() ;
    this.camera.position.fromArray(this.default_camera_pos);
    this.controls = new Controls(this.camera, this.target);
  }
  this.set_common_key_bindings();
  if (this.constructor === Viewer) { this.set_real_space_key_bindings(); }

  function get_elem(name) {
    if (options[name] === null || typeof document === 'undefined') { return null; }
    return document.getElementById(options[name] || name);
  }
  this.hud_el = get_elem('hud');
  this.container = get_elem('viewer');
  this.help_el = get_elem('help');
  if (this.hud_el) {
    if (this.hud_el.innerHTML === '') { this.hud_el.innerHTML = INIT_HUD_TEXT; }
    this.initial_hud_html = this.hud_el.innerHTML;
  }

  try {
    this.renderer = new WebGLRenderer({antialias: true});
  } catch (e) {
    this.hud('No WebGL in your browser?', 'ERR');
    this.renderer = null;
    return;
  }

  if (this.container == null) { return; } // can be null in headless tests
  this.renderer.setClearColor(this.config.colors.bg, 1);
  this.renderer.setPixelRatio(window.devicePixelRatio);
  this.resize();
  this.camera.zoom = this.camera.right / 35.0;// arbitrary choice
  this.update_camera();
  var el = this.renderer.domElement;
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
  el.addEventListener('wheel', this.wheel.bind(this));
  el.addEventListener('mousedown', this.mousedown.bind(this));
  el.addEventListener('touchstart', this.touchstart.bind(this));
  el.addEventListener('touchmove', this.touchmove.bind(this));
  el.addEventListener('touchend', this.touchend.bind(this));
  el.addEventListener('touchcancel', this.touchend.bind(this));
  el.addEventListener('dblclick', this.dblclick.bind(this));

  var self = this;

  this.mousemove = function (event) {
    event.preventDefault();
    //event.stopPropagation();
    self.controls.move(self.relX(event), self.relY(event));
  };

  this.mouseup = function (event) {
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

Viewer.prototype.pick_atom = function pick_atom (coords, camera) {
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

Viewer.prototype.set_colors = function set_colors () {
  var scheme = this.ColorSchemes[this.config.color_scheme];
  if (!scheme) { throw Error('Unknown color scheme.'); }
  this.decor.zoom_grid.material.uniforms.ucolor.value.set(scheme.fg);
  this.config.colors = scheme;
  this.redraw_all();
};

// relative position on canvas in normalized device coordinates [-1, +1]
Viewer.prototype.relX = function relX (evt) {
  return 2 * (evt.pageX - this.window_offset[0]) / this.window_size[0] - 1;
};

Viewer.prototype.relY = function relY (evt) {
  return 1 - 2 * (evt.pageY - this.window_offset[1]) / this.window_size[1];
};

Viewer.prototype.hud = function hud (text, type) {
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

Viewer.prototype.redraw_center = function redraw_center (force) {
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
    });
    this.scene.add(this.decor.mark);
  }
};

Viewer.prototype.redraw_maps = function redraw_maps (force) {
  this.redraw_center(force);
  var r = this.config.map_radius;
  for (var i = 0, list = this.map_bags; i < list.length; i += 1) {
    var map_bag = list[i];

      if (force || this.target.distanceToSquared(map_bag.block_ctr) > r/100) {
      this.redraw_map(map_bag);
    }
  }
};

Viewer.prototype.remove_and_dispose = function remove_and_dispose (obj) {
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

Viewer.prototype.clear_el_objects = function clear_el_objects (map_bag) {
  for (var i = 0, list = map_bag.el_objects; i < list.length; i += 1) {
    var o = list[i];

      this.remove_and_dispose(o);
  }
  map_bag.el_objects = [];
};

Viewer.prototype.clear_model_objects = function clear_model_objects (model_bag) {
  for (var i = 0, list = model_bag.objects; i < list.length; i += 1) {
    var o = list[i];

      this.remove_and_dispose(o);
  }
  model_bag.objects = [];
};

Viewer.prototype.has_frag_depth = function has_frag_depth () {
  return this.renderer && this.renderer.extensions.get('EXT_frag_depth');
};

Viewer.prototype.set_model_objects = function set_model_objects (model_bag) {
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
Viewer.prototype.toggle_label = function toggle_label (pick, show) {
  if (pick.atom == null) { return; }
  var text = pick.atom.short_label();
  var uid = text; // we assume that the labels inside one model are unique
  var is_shown = (uid in this.labels);
  if (show === undefined) { show = !is_shown; }
  if (show) {
    if (is_shown) { return; }
    var atom_style = pick.atom.is_ligand ? 'ligand_style' : 'render_style';
    var balls = pick.bag && pick.bag.conf[atom_style] === 'ball&stick';
    var label = new Label(text, {
      pos: pick.atom.xyz,
      font: this.config.label_font,
      color: '#' + this.config.colors.fg.getHexString(),
      win_size: this.window_size,
      z_shift: balls ? this.config.ball_size + 0.1 : 0.2,
    });
    if (pick.bag == null || label.mesh == null) { return; }
    this.labels[uid] = { o: label, bag: pick.bag };
    this.scene.add(label.mesh);
  } else {
    if (!is_shown) { return; }
    this.remove_and_dispose(this.labels[uid].o.mesh);
    delete this.labels[uid];
  }
};

Viewer.prototype.redraw_labels = function redraw_labels () {
  for (var uid in this.labels) { // eslint-disable-line guard-for-in
    var text = uid;
    this.labels[uid].o.redraw(text, {
      font: this.config.label_font,
      color: '#' + this.config.colors.fg.getHexString(),
    });
  }
};

Viewer.prototype.toggle_map_visibility = function toggle_map_visibility (map_bag) {
  if (typeof map_bag === 'number') {
    map_bag = this.map_bags[map_bag];
  }
  map_bag.visible = !map_bag.visible;
  this.redraw_map(map_bag);
  this.request_render();
};

Viewer.prototype.redraw_map = function redraw_map (map_bag) {
  this.clear_el_objects(map_bag);
  if (map_bag.visible) {
    map_bag.map.block.clear();
    this.add_el_objects(map_bag);
  }
};

Viewer.prototype.toggle_model_visibility = function toggle_model_visibility (model_bag, visible) {
  model_bag = model_bag || this.selected.bag;
  if (model_bag == null) { return; }
  model_bag.visible = visible == null ? !model_bag.visible : visible;
  this.redraw_model(model_bag);
  this.request_render();
};

Viewer.prototype.redraw_model = function redraw_model (model_bag) {
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

Viewer.prototype.add_el_objects = function add_el_objects (map_bag) {
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

Viewer.prototype.change_isolevel_by = function change_isolevel_by (map_idx, delta) {
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

Viewer.prototype.change_map_radius = function change_map_radius (delta) {
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

Viewer.prototype.change_slab_width_by = function change_slab_width_by (delta) {
  var slab_width = this.controls.slab_width;
  slab_width[0] = Math.max(slab_width[0] + delta, 0.01);
  slab_width[1] = Math.max(slab_width[1] + delta, 0.01);
  this.update_camera();
  var final_width = this.camera.far - this.camera.near;
  this.hud('clip width: ' + final_width.toPrecision(3));
};

Viewer.prototype.change_zoom_by_factor = function change_zoom_by_factor (mult) {
  this.camera.zoom *= mult;
  this.update_camera();
  this.hud('zoom: ' + this.camera.zoom.toPrecision(3));
};

Viewer.prototype.change_bond_line = function change_bond_line (delta) {
  this.config.bond_line = Math.max(this.config.bond_line + delta, 0.1);
  this.redraw_models();
  this.hud('bond width: ' + scale_by_height(this.config.bond_line,
                                            this.window_size).toFixed(1));
};

Viewer.prototype.change_map_line = function change_map_line (delta) {
  this.config.map_line = Math.max(this.config.map_line + delta, 0.1);
  this.redraw_maps(true);
  this.hud('wireframe width: ' + this.config.map_line.toFixed(1));
};

Viewer.prototype.toggle_full_screen = function toggle_full_screen () {
  var d = document;
  // @ts-expect-error no mozFullScreenElement
  if (d.fullscreenElement || d.mozFullScreenElement ||
      // @ts-expect-error no msFullscreenElement
      d.webkitFullscreenElement || d.msFullscreenElement) {
    // @ts-expect-error no webkitExitFullscreen
    var ex = d.exitFullscreen || d.webkitExitFullscreen ||
               // @ts-expect-error no msExitFullscreen
               d.mozCancelFullScreen || d.msExitFullscreen;
    if (ex) { ex.call(d); }
  } else {
    var el = this.container;
    if (!el) { return; }
    // @ts-expect-error no webkitRequestFullscreen
    var req = el.requestFullscreen || el.webkitRequestFullscreen ||
                // @ts-expect-error no msRequestFullscreen
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

Viewer.prototype.get_cell_box_func = function get_cell_box_func () {
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

Viewer.prototype.shift_clip = function shift_clip (delta) {
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

Viewer.prototype.select_next = function select_next (info, key, options, back) {
  var old_idx = options.indexOf(this.config[key]);
  var len = options.length;
  var new_idx = (old_idx + (back ? len - 1 : 1)) % len;
  this.config[key] = options[new_idx];
  var html = info + ':';
  for (var i = 0; i < len; i++) {
    var tag = (i === new_idx ? 'u' : 's');
    html += ' <' + tag + '>' + options[i] + '</' + tag + '>';
  }
  this.hud(html, 'HTML');
};

Viewer.prototype.keydown = function keydown (evt) {
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
  kb[66] = function ( evt) {
    var schemes = Object.keys(this.ColorSchemes);
    this.select_next('color scheme', 'color_scheme', schemes, evt.shiftKey);
    this.set_colors();
  };
  // c
  kb[67] = function ( evt) {
    this.select_next('coloring by', 'color_prop', COLOR_PROPS, evt.shiftKey);
    this.redraw_models();
  };
  // e
  kb[69] = function () {
    var fog = this.scene.fog;
    var has_fog = (fog.far === 1);
    fog.far = (has_fog ? 1e9 : 1);
    this.hud((has_fog ? 'dis': 'en') + 'able fog');
    this.redraw_all();
  };
  // h
  kb[72] = this.toggle_help;
  // i
  kb[73] = function ( evt) {
    this.hud('toggled spinning');
    this.controls.toggle_auto(evt.shiftKey);
  };
  // k
  kb[75] = function () {
    this.hud('toggled rocking');
    this.controls.toggle_auto(0.0);
  };
  // m
  kb[77] = function ( evt) {
    this.change_zoom_by_factor(evt.shiftKey ? 1.2 : 1.03);
  };
  // n
  kb[78] = function ( evt) {
    this.change_zoom_by_factor(1 / (evt.shiftKey ? 1.2 : 1.03));
  };
  // q
  kb[81] = function ( evt) {
    this.select_next('label font', 'label_font', LABEL_FONTS, evt.shiftKey);
    this.redraw_labels();
  };
  // r
  kb[82] = function ( evt) {
    if (evt.shiftKey) {
      this.hud('redraw!');
      this.redraw_all();
    } else {
      this.hud('recentered');
      this.recenter();
    }
  };
  // w
  kb[87] = function ( evt) {
    this.select_next('map style', 'map_style', MAP_STYLES, evt.shiftKey);
    this.redraw_maps(true);
  };
  // add, equals/firefox, equal sign
  kb[107] = kb[61] = kb[187] = function ( evt) {
    this.change_isolevel_by(evt.shiftKey ? 1 : 0, 0.1);
  };
  // subtract, minus/firefox, dash
  kb[109] = kb[173] = kb[189] = function ( evt) {
    this.change_isolevel_by(evt.shiftKey ? 1 : 0, -0.1);
  };
  // [
  kb[219] = function () { this.change_map_radius(-2); };
  // ]
  kb[221] = function () { this.change_map_radius(2); };
  // \ (backslash)
  kb[220] = function ( evt) {
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
  kb[36] = function ( evt) {
    evt.shiftKey ? this.change_map_line(0.1) : this.change_bond_line(0.2);
  };
  // End
  kb[35] = function ( evt) {
    evt.shiftKey ? this.change_map_line(-0.1) : this.change_bond_line(-0.2);
  };
  // Space
  kb[32] = function ( evt) {
    this.center_next_residue(evt.shiftKey);
  };
  // d
  kb[68] = function () {
    this.change_slab_width_by(-0.1);
  };
  // f
  kb[70] = function ( evt) {
    evt.shiftKey ? this.toggle_full_screen() : this.change_slab_width_by(0.1);
  };
  // l
  kb[76] = function ( evt) {
    this.select_next('ligands as', 'ligand_style', LIGAND_STYLES, evt.shiftKey);
    this.redraw_models();
  };
  // p
  kb[80] = function ( evt) {
    evt.shiftKey ? this.permalink() : this.go_to_nearest_Ca();
  };
  // s
  kb[83] = function ( evt) {
    this.select_next('rendering as', 'render_style', RENDER_STYLES, evt.shiftKey);
    this.redraw_models();
  };
  // t
  kb[84] = function ( evt) {
    this.select_next('waters as', 'water_style', WATER_STYLES, evt.shiftKey);
    this.redraw_models();
  };
  // u
  kb[85] = function () {
    this.hud('toggled unit cell box');
    this.toggle_cell_box();
  };
  // v
  kb[86] = function () {
    this.toggle_inactive_models();
  };
  // y
  kb[89] = function () {
    this.config.hydrogens = !this.config.hydrogens;
    this.hud((this.config.hydrogens ? 'show' : 'hide') +
             ' hydrogens (if any)');
    this.redraw_models();
  };
  // comma
  kb[188] = function ( evt) {
    if (evt.shiftKey) { this.shift_clip(1); }
  };
  // period
  kb[190] = function ( evt) {
    if (evt.shiftKey) { this.shift_clip(-1); }
  };
};

Viewer.prototype.mousedown = function mousedown (event) {
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

Viewer.prototype.dblclick = function dblclick (event) {
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

Viewer.prototype.touchstart = function touchstart (event) {
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

Viewer.prototype.touchmove = function touchmove (event) {
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

Viewer.prototype.wheel = function wheel (evt) {
  evt.preventDefault();
  evt.stopPropagation();
  this.mousewheel_action(evt.deltaY, evt);
  this.request_render();
};

// overrided in ReciprocalViewer
Viewer.prototype.mousewheel_action = function mousewheel_action (delta, evt) {
  var map_idx = evt.shiftKey ? 1 : 0;
  this.change_isolevel_by(map_idx, 0.0005 * delta);
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
Viewer.prototype.recenter = function recenter (xyz, cam, steps) {
  var bag = this.selected.bag;
  var new_up = new Vector3(0, 1, 0);
  var vec_cam;
  var vec_xyz;
  var eye;
  if (xyz != null && cam == null && bag != null) {
    // look from specified point toward the center of the molecule,
    // i.e. shift camera away from the molecule center.
    var mc = bag.model.get_center();
    eye = new Vector3(xyz[0] - mc[0], xyz[1] - mc[1], xyz[2] - mc[2]);
    eye.setLength(100);
    vec_xyz = new Vector3(xyz[0], xyz[1], xyz[2]);
    vec_cam = eye.clone().add(vec_xyz);
  } else {
    if (xyz == null) {
      if (bag != null) {
        xyz = bag.model.get_center();
      } else {
        var uc_func = this.get_cell_box_func();
        xyz = uc_func ? uc_func([0.5, 0.5, 0.5]) : [0, 0, 0];
      }
    }
    vec_xyz = new Vector3(xyz[0], xyz[1], xyz[2]);
    if (cam != null) {
      vec_cam = new Vector3(cam[0], cam[1], cam[2]);
      eye = vec_cam.clone().sub(vec_xyz);
      new_up.copy(this.camera.up); // preserve the up direction
    } else {
      var dc = this.default_camera_pos;
      vec_cam = new Vector3(xyz[0] + dc[0], xyz[1] + dc[1], xyz[2] + dc[2]);
    }
  }
  if (eye != null) {
    new_up.projectOnPlane(eye);
    if (new_up.lengthSq() < 0.0001) { new_up.x += 1; }
    new_up.normalize();
  }
  this.controls.go_to(vec_xyz, vec_cam, new_up, steps);
};

Viewer.prototype.center_next_residue = function center_next_residue (back) {
  var bag = this.selected.bag;
  if (bag == null) { return; }
  var atom = bag.model.next_residue(this.selected.atom, back);
  if (atom != null) {
    this.select_atom({bag: bag, atom: atom}, {steps: 30});
  }
};

Viewer.prototype.select_atom = function select_atom (pick, options) {
    if ( options === void 0 ) options={};

  this.hud('-> ' + pick.bag.label + ' ' + pick.atom.long_label());
  var xyz = pick.atom.xyz;
  this.controls.go_to(new Vector3(xyz[0], xyz[1], xyz[2]),
                      null, null, options.steps);
  this.toggle_label(this.selected, false);
  this.selected = pick;
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
  //if (this.nav) {
  //this.nav.renderer.render(this.nav.scene, this.camera);
  //}
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

Viewer.prototype.add_model = function add_model (model, options) {
    if ( options === void 0 ) options={};

  var model_bag = new ModelBag(model, this.config, this.window_size);
  model_bag.hue_shift = options.hue_shift || 0.06 * this.model_bags.length;
  this.model_bags.push(model_bag);
  this.set_model_objects(model_bag);
  this.request_render();
};

Viewer.prototype.add_map = function add_map (map, is_diff_map) {
  //map.show_debug_info();
  var map_bag = new MapBag(map, this.config, is_diff_map);
  this.map_bags.push(map_bag);
  this.add_el_objects(map_bag);
  this.request_render();
};

Viewer.prototype.load_file = function load_file (url, options,
          callback ) {
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
    req.addEventListener('progress', function (evt) {
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

Viewer.prototype.set_dropzone = function set_dropzone (zone, callback) {
  var self = this;
  zone.addEventListener('dragover', function (e) {
    e.stopPropagation();
    e.preventDefault();
    if (e.dataTransfer != null) { e.dataTransfer.dropEffect = 'copy'; }
    self.hud('ready for file drop...');
  });
  zone.addEventListener('drop', function (e) {
    e.stopPropagation();
    e.preventDefault();
    if (e.dataTransfer == null) { return; }
    var names = [];
    for (var i = 0; i < e.dataTransfer.files.length; i++) {
      var file = e.dataTransfer.files.item(i);
      try {
        self.hud('Loading ' + file.name);
        callback(file);
      } catch (e$1) {
        self.hud('Loading ' + file.name + ' failed.\n' + e$1.message, 'ERR');
        return;
      }
      names.push(file.name);
    }
    self.hud('loaded ' + names.join(', '));
  });
};

// for use with set_dropzone
Viewer.prototype.pick_pdb_and_map = function pick_pdb_and_map (file) {
  var self = this;
  var reader = new FileReader();
  if (/\.(pdb|ent)$/.test(file.name)) {
    reader.onload = function (evt) {
      self.load_pdb_from_text(evt.target.result );
      self.recenter();
    };
    reader.readAsText(file);
  } else if (/\.(map|ccp4|mrc|dsn6|omap)$/.test(file.name)) {
    var map_format = /\.(dsn6|omap)$/.test(file.name) ? 'dsn6' : 'ccp4';
    reader.onloadend = function (evt) {
      if (evt.target != null && evt.target.readyState == 2) {
        self.load_map_from_buffer(evt.target.result ,
                                  {format: map_format});
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

Viewer.prototype.set_view = function set_view (options) {
  var frag = parse_url_fragment();
  if (frag.zoom) { this.camera.zoom = frag.zoom; }
  this.recenter(frag.xyz || (options && options.center), frag.eye, 1);
};

// Load molecular model from PDB file and centers the view
Viewer.prototype.load_pdb_from_text = function load_pdb_from_text (text) {
  var len = this.model_bags.length;
  var models = modelsFromPDB(text);
  for (var i = 0, list = models; i < list.length; i += 1) {
    var model = list[i];

      this.add_model(model);
  }
  this.selected.bag = this.model_bags[len];
};

Viewer.prototype.load_structure_from_buffer = function load_structure_from_buffer (gemmi, buffer, name) {
  var len = this.model_bags.length;
  var models = modelsFromGemmi(gemmi, buffer, name);
  for (var i = 0, list = models; i < list.length; i += 1) {
    var model = list[i];

      this.add_model(model);
  }
  this.selected.bag = this.model_bags[len];
};

Viewer.prototype.load_pdb = function load_pdb (url, options,
         callback) {
  var self = this;
  var gemmi = options && options.gemmi;
  this.load_file(url, {binary: !!gemmi, progress: true}, function (req) {
    var t0 = performance.now();
    if (gemmi) {
      self.load_structure_from_buffer(gemmi, req.response, url);
    } else {
      self.load_pdb_from_text(req.responseText);
    }
    console.log('coordinate file processed in', (performance.now() - t0).toFixed(2),
                gemmi ? 'ms (using gemmi)': 'ms');
    if (options == null || !options.stay) { self.set_view(options); }
    if (callback) { callback(); }
  });
};

Viewer.prototype.load_map = function load_map (url, options,
         callback) {
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

Viewer.prototype.load_map_from_buffer = function load_map_from_buffer (buffer, options) {
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
Viewer.prototype.load_maps = function load_maps (url1, url2,
          options, callback) {
  var format = options.format || 'ccp4';
  var self = this;
  this.load_map(url1, {diff_map: false, format: format}, function () {
    self.load_map(url2, {diff_map: true, format: format}, callback);
  });
};

// Load a model (PDB), normal map and a difference map - in this order.
Viewer.prototype.load_pdb_and_maps = function load_pdb_and_maps (pdb, map1, map2,
                  options, callback) {
  var self = this;
  this.load_pdb(pdb, options, function () {
    self.load_maps(map1, map2, options, callback);
  });
};

// for backward compatibility:
Viewer.prototype.load_ccp4_maps = function load_ccp4_maps (url1, url2, callback) {
  this.load_maps(url1, url2, {format: 'ccp4'}, callback);
};
Viewer.prototype.load_pdb_and_ccp4_maps = function load_pdb_and_ccp4_maps (pdb, map1, map2,
                       callback) {
  this.load_pdb_and_maps(pdb, map1, map2, {format: 'ccp4'}, callback);
};

// pdb_id here should be lowercase ('1abc')
Viewer.prototype.load_from_pdbe = function load_from_pdbe (pdb_id, callback) {
  var id = pdb_id.toLowerCase();
  this.load_pdb_and_maps(
    'https://www.ebi.ac.uk/pdbe/entry-files/pdb' + id + '.ent',
    'https://www.ebi.ac.uk/pdbe/coordinates/files/' + id + '.ccp4',
    'https://www.ebi.ac.uk/pdbe/coordinates/files/' + id + '_diff.ccp4',
    {format: 'ccp4'}, callback);
};
Viewer.prototype.load_from_rcsb = function load_from_rcsb (pdb_id, callback) {
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
  // @ts-expect-error Cannot find name 'VERSION'
  (typeof VERSION === 'string' ? VERSION : 'dev'); // eslint-disable-line

Viewer.prototype.ColorSchemes = ColorSchemes$1;

function to_col(num) { return new Color(num); }

var ColorSchemes = {
  'solarized dark': {
    bg: new Color(0x002b36),
    fg: new Color(0xfdf6e3),
    map_den: new Color(0xeee8d5),
    center: new Color(0xfdf6e3),
    lattices: [0xdc322f, 0x2aa198, 0x268bd2, 0x859900,
               0xd33682, 0xb58900, 0x6c71c4, 0xcb4b16].map(to_col),
    axes: [0xffaaaa, 0xaaffaa, 0xaaaaff].map(to_col),
  },
  'solarized light': {
    bg: new Color(0xfdf6e3),
    fg: new Color(0x002b36),
    map_den: new Color(0x073642),
    center: new Color(0x002b36),
    lattices: [0xdc322f, 0x2aa198, 0x268bd2, 0x859900,
               0xd33682, 0xb58900, 0x6c71c4, 0xcb4b16].map(to_col),
    axes: [0xffaaaa, 0xaaffaa, 0xaaaaff].map(to_col),
  },
};

// options handled by Viewer#select_next()
var SPOT_SEL = ['all', 'unindexed', '#1']; //extended when needed
var SHOW_AXES = ['two', 'three', 'none'];
var SPOT_SHAPES = ['wheel', 'square'];

// Modified ElMap for handling output of dials.rs_mapper.
// rs_mapper outputs map in ccp4 format, but we need to rescale it,
// shift it so the box is centered at 0,0,0,
// and the translational symmetry doesn't apply.
var ReciprocalSpaceMap = /*@__PURE__*/(function (ElMap) {
  function ReciprocalSpaceMap(buf) {
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

  ReciprocalSpaceMap.prototype.extract_block = function extract_block (radius, center) {
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
  function ReciprocalViewer(options) {
    if ( options === void 0 ) options = {};

    options.color_scheme = 'solarized dark';
    Viewer.call(this, options);
    this.default_camera_pos = [100, 0, 0];
    this.axes = null;
    this.points = null;
    this.max_dist = -1;
    this.d_min = -1;
    this.d_max_inv = 0;
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
    kb[80] = function () { this.permalink(); };
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
    kb[107] = kb[61] = kb[187] = function () {
      this.change_isolevel_by(0, 0.01);
    };
    // subtract, minus/firefox, dash
    kb[109] = kb[173] = kb[189] = function () {
      this.change_isolevel_by(0, -0.01);
    };
    // [
    kb[219] = function () { this.change_map_radius(-0.001); };
    // ]
    kb[221] = function () { this.change_map_radius(0.001); };
  };

  ReciprocalViewer.prototype.file_drop_callback = function file_drop_callback (file) {
    var self = this;
    var reader = new FileReader();
    if (/\.(map|ccp4)$/.test(file.name)) {
      reader.onloadend = function (evt) {
        if (evt.target.readyState == 2) {
          self.load_map_from_ab(evt.target.result );
        }
      };
      reader.readAsArrayBuffer(file);
    } else {
      reader.onload = function (evt) {
        self.load_from_string(evt.target.result , {});
      };
      reader.readAsText(file);
    }
  };

  ReciprocalViewer.prototype.load_data = function load_data (url, options) {
    if ( options === void 0 ) options = {};

    var self = this;
    this.load_file(url, {binary: false, progress: true}, function (req) {
      var ok = self.load_from_string(req.responseText, options);
      if (ok && options.callback) { options.callback(); }
    });
  };

  ReciprocalViewer.prototype.load_from_string = function load_from_string (text, options) {
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

  ReciprocalViewer.prototype.load_map_from_ab = function load_map_from_ab (buffer) {
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
    });
    this.axes = makeLineSegments(material, vertices, colors);
    this.scene.add(this.axes);
  };

  ReciprocalViewer.prototype.set_points = function set_points (data) {
    if (this.points != null) {
      this.remove_and_dispose(this.points);
      this.points = null;
    }
    if (data == null || data.lattice_ids == null || data.pos == null) { return; }
    var color_arr = new Float32Array(3 * data.lattice_ids.length);
    this.colorize_by_id(color_arr, data.lattice_ids);
    var geometry = new BufferGeometry();
    geometry.setAttribute('position', new BufferAttribute(data.pos, 3));
    geometry.setAttribute('color', new BufferAttribute(color_arr, 3));
    var groups = new Float32Array(data.lattice_ids);
    geometry.setAttribute('group', new BufferAttribute(groups, 1));
    this.points = new Points(geometry, this.point_material);
    this.scene.add(this.points);
    this.request_render();
  };

  ReciprocalViewer.prototype.colorize_by_id = function colorize_by_id (color_arr, group_id) {
    var palette = this.config.colors.lattices;
    for (var i = 0; i < group_id.length; i++) {
      var c = palette[(group_id[i] + 1) % 4];
      color_arr[3*i] = c.r;
      color_arr[3*i+1] = c.g;
      color_arr[3*i+2] = c.b;
    }
  };

  ReciprocalViewer.prototype.mousewheel_action = function mousewheel_action (delta) {
    this.change_zoom_by_factor(1 + 0.0005 * delta);
  };

  ReciprocalViewer.prototype.change_point_size = function change_point_size (delta) {
    var size = this.point_material.uniforms.size;
    size.value = Math.max(size.value + delta, 0.5);
    this.hud('point size: ' + size.value.toFixed(1));
  };

  ReciprocalViewer.prototype.change_dmin = function change_dmin (delta) {
    this.d_min = Math.max(this.d_min + delta, 0.1);
    var dmax = this.d_max_inv > 0 ? 1 / this.d_max_inv : null;
    if (dmax !== null && this.d_min > dmax) { this.d_min = dmax; }
    this.point_material.uniforms.r2_max.value = 1 / (this.d_min * this.d_min);
    var low_res = dmax !== null ? dmax.toFixed(2) : '∞';
    this.hud('res. limit: ' + low_res + ' - ' + this.d_min.toFixed(2) + 'Å');
  };

  ReciprocalViewer.prototype.change_dmax = function change_dmax (delta) {
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
    // here the map is ReciprocalSpaceMap not ElMap
    var a = this.map_bags[0].map.box_size;
    return function (xyz) {
      return [(xyz[0]-0.5) * a[0], (xyz[1]-0.5) * a[1], (xyz[2]-0.5) * a[2]];
    };
  };

  return ReciprocalViewer;
}(Viewer));

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

ReciprocalViewer.prototype.ColorSchemes = ColorSchemes;

function log_timing(t0, text) {
  console.log(text + ': ' + (performance.now() - t0).toFixed(2) + ' ms.');
}

function add_map_from_mtz(viewer, mtz, map_data, is_diff) {
  var map = new ElMap();
  var mc = mtz.cell;
  map.unit_cell = new UnitCell(mc.a, mc.b, mc.c, mc.alpha, mc.beta, mc.gamma);
  map.stats.rms = mtz.rmsd;
  map.grid = new GridArray([mtz.nx, mtz.ny, mtz.nz]);
  map.grid.values.set(map_data);
  viewer.add_map(map, is_diff);
}

function load_maps_from_mtz_buffer(viewer, mtz,
                                   labels) {
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

function load_maps_from_mtz(gemmi, viewer, url,
                            labels, callback) {
  viewer.load_file(url, {binary: true, progress: true}, function (req) {
    var t0 = performance.now();
    try {
      var mtz = gemmi.readMtz(req.response);
      //console.log("[after readMTZ] wasm mem:", gemmi.HEAPU8.length / 1024, "kb");
      load_maps_from_mtz_buffer(viewer, mtz, labels);
    } catch (e) {
      viewer.hud(e.message, 'ERR');
      return;
    }
    log_timing(t0, 'load_maps_from_mtz');
    //console.log("wasm mem:", gemmi.HEAPU8.length / 1024, "kb");
    if (callback) { callback(); }
  });
}

function set_pdb_and_mtz_dropzone(gemmi, viewer,
                                  zone) {
  viewer.set_dropzone(zone, function (file) {
    if (/\.mtz$/.test(file.name)) {
      var reader = new FileReader();
      reader.onloadend = function (evt) {
        if (evt.target.readyState == 2) {
          var t0 = performance.now();
          try {
            var mtz = gemmi.readMtz(evt.target.result);
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
exports.Color = Color;
exports.Controls = Controls;
exports.ElMap = ElMap;
exports.Fog = Fog;
exports.GridArray = GridArray;
exports.Label = Label;
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
exports.UnitCell = UnitCell;
exports.Vector3 = Vector3;
exports.Viewer = Viewer;
exports.WebGLRenderer = WebGLRenderer;
exports.addXyzCross = addXyzCross;
exports.fog_end_fragment = fog_end_fragment;
exports.fog_pars_fragment = fog_pars_fragment;
exports.load_maps_from_mtz = load_maps_from_mtz;
exports.load_maps_from_mtz_buffer = load_maps_from_mtz_buffer;
exports.makeBalls = makeBalls;
exports.makeChickenWire = makeChickenWire;
exports.makeCube = makeCube;
exports.makeGrid = makeGrid;
exports.makeLineMaterial = makeLineMaterial;
exports.makeLineSegments = makeLineSegments;
exports.makeLines = makeLines;
exports.makeMultiColorLines = makeMultiColorLines;
exports.makeRgbBox = makeRgbBox;
exports.makeRibbon = makeRibbon;
exports.makeSticks = makeSticks;
exports.makeUniforms = makeUniforms;
exports.makeWheels = makeWheels;
exports.modelsFromGemmi = modelsFromGemmi;
exports.modelsFromPDB = modelsFromPDB;
exports.set_pdb_and_mtz_dropzone = set_pdb_and_mtz_dropzone;

}));
