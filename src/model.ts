import { UnitCell } from './unitcell';

type Num3 = [number, number, number];

const AMINO_ACIDS = [
  'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
  'LYS', 'MET', 'MSE', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'UNK',
];
const NUCLEIC_ACIDS = [
  'DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U', 'rA', 'rC', 'rG', 'rU',
  'Ar', 'Cr', 'Gr', 'Ur',
];

const NOT_LIGANDS = ['HOH'].concat(AMINO_ACIDS, NUCLEIC_ACIDS);

export function modelsFromPDB(pdb_string: string) {
  const models = [new Model()];
  let pdb_tail = models[0].from_pdb(pdb_string.split('\n'));
  while (pdb_tail != null) {
    const model = new Model();
    pdb_tail = model.from_pdb(pdb_tail);
    if (model.atoms.length > 0) models.push(model);
  }
  return models;
}

export class Model {
  atoms: Atom[];
  unit_cell: UnitCell | null;
  space_group: null;
  has_hydrogens: boolean;
  lower_bound: Num3;
  upper_bound: Num3;
  residue_map: {[id:string]: Atom[]} | null;
  cubes: Cubicles | null;

  constructor() {
    this.atoms = [];
    this.unit_cell = null;
    this.space_group = null;
    this.has_hydrogens = false;
    this.lower_bound = [0, 0, 0];
    this.upper_bound = [0, 0, 0];
    this.residue_map = null;
    this.cubes = null;
  }

  from_pdb(pdb_lines: string[]): string[] | null {
    let chain_index = 0;  // will be ++'ed for the first atom
    let last_chain = null;
    let atom_i_seq = 0;
    let continuation = null;
    for (let i = 0; i < pdb_lines.length; i++) {
      const line = pdb_lines[i];
      const rec_type = line.substring(0, 6).toUpperCase();
      if (rec_type === 'ATOM  ' || rec_type === 'HETATM') {
        const new_atom = new Atom();
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
      } else if (rec_type === 'ANISOU') {
        // we don't use anisotropic B-factors
      } else if (rec_type === 'CRYST1') {
        const a = parseFloat(line.substring(6, 15));
        const b = parseFloat(line.substring(15, 24));
        const c = parseFloat(line.substring(24, 33));
        const alpha = parseFloat(line.substring(33, 40));
        const beta = parseFloat(line.substring(40, 47));
        const gamma = parseFloat(line.substring(47, 54));
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
    if (this.atoms.length === 0) throw Error('No atom records found.');
    this.calculate_bounds();
    this.calculate_connectivity();
    return continuation;
  }

  calculate_bounds() {
    const lower = this.lower_bound = [Infinity, Infinity, Infinity];
    const upper = this.upper_bound = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < this.atoms.length; i++) {
      const atom = this.atoms[i];
      for (let j = 0; j < 3; j++) {
        const v = atom.xyz[j];
        if (v < lower[j]) lower[j] = v;
        if (v > upper[j]) upper[j] = v;
      }
    }
    // with a margin
    for (let k = 0; k < 3; ++k) {
      lower[k] -= 0.001;
      upper[k] += 0.001;
    }
  }

  next_residue(atom: Atom | null, backward: boolean) {
    const len = this.atoms.length;
    const start = (atom ? atom.i_seq : 0) + len;  // +len to avoid idx<0 below
    for (let i = (atom ? 1 : 0); i < len; i++) {
      const idx = (start + (backward ? -i : i)) % len;
      const a = this.atoms[idx];
      if (!a.is_main_conformer()) continue;
      if ((a.name === 'CA' && a.element === 'C') || a.name === 'P') {
        return a;
      }
    }
  }

  extract_trace() {
    const segments = [];
    let current_segment = [];
    let last_atom = null;
    for (let i = 0; i < this.atoms.length; i++) {
      const atom = this.atoms[i];
      if (atom.altloc !== '' && atom.altloc !== 'A') continue;
      if ((atom.name === 'CA' && atom.element === 'C') || atom.name === 'P') {
        let start_new = true;
        if (last_atom !== null && last_atom.chain_index === atom.chain_index) {
          const dxyz2 = atom.distance_sq(last_atom);
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
  }

  get_residues() {
    if (this.residue_map != null) return this.residue_map;
    const residues = {};
    for (let i = 0; i < this.atoms.length; i++) {
      const atom = this.atoms[i];
      const resid = atom.resid();
      const reslist = residues[resid];
      if (reslist === undefined) {
        residues[resid] = [atom];
      } else {
        reslist.push(atom);
      }
    }
    this.residue_map = residues;
    return residues;
  }

  // tangent vector to the ribbon representation
  calculate_tangent_vector(residue: Atom[]): Num3 {
    let a1 = null;
    let a2 = null;
    // it may be too simplistic
    const peptide = (residue[0].resname.length === 3);
    const name1 = peptide ? 'C' : 'C2\'';
    const name2 = peptide ? 'O' : 'O4\'';
    for (let i = 0; i < residue.length; i++) {
      const atom = residue[i];
      if (!atom.is_main_conformer()) continue;
      if (atom.name === name1) {
        a1 = atom.xyz;
      } else if (atom.name === name2) {
        a2 = atom.xyz;
      }
    }
    if (a1 === null || a2 === null) return [0, 0, 1]; // arbitrary value
    const d = [a1[0]-a2[0], a1[1]-a2[1], a1[2]-a2[2]];
    const len = Math.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    return [d[0]/len, d[1]/len, d[2]/len];
  }

  get_center(): Num3 {
    let xsum = 0, ysum = 0, zsum = 0;  // eslint-disable-line
    const n_atoms = this.atoms.length;
    for (let i = 0; i < n_atoms; i++) {
      const xyz = this.atoms[i].xyz;
      xsum += xyz[0];
      ysum += xyz[1];
      zsum += xyz[2];
    }
    return [xsum / n_atoms, ysum / n_atoms, zsum / n_atoms];
  }

  calculate_connectivity() {
    const atoms = this.atoms;
    const cubes = new Cubicles(atoms, 3.0, this.lower_bound, this.upper_bound);
    //let cnt = 0;
    for (let i = 0; i < cubes.boxes.length; i++) {
      const box = cubes.boxes[i];
      if (box.length === 0) continue;
      const nearby_atoms = cubes.get_nearby_atoms(i);
      for (let a = 0; a < box.length; a++) {
        const atom_id = box[a];
        for (let k = 0; k < nearby_atoms.length; k++) {
          const j = nearby_atoms[k];
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
  }

  get_nearest_atom(x: number, y: number, z: number, atom_name?: string) {
    const cubes = this.cubes;
    if (cubes == null) throw Error('Missing Cubicles');
    const box_id = cubes.find_box_id(x, y, z);
    const indices = cubes.get_nearby_atoms(box_id);
    let nearest = null;
    let min_d2 = Infinity;
    for (let i = 0; i < indices.length; i++) {
      const atom = this.atoms[indices[i]];
      if (atom_name != null && atom_name !== atom.name) continue;
      const dx = atom.xyz[0] - x;
      const dy = atom.xyz[1] - y;
      const dz = atom.xyz[2] - z;
      const d2 = dx*dx + dy*dy + dz*dz;
      if (d2 < min_d2) {
        nearest = atom;
        min_d2 = d2;
      }
    }
    return nearest;
  }
}

// Single atom and associated labels
class Atom {
  hetero: boolean;
  name: string;
  altloc: string;
  resname: string;
  chain: string;
  chain_index: number;
  resseq: number;
  icode: string | null;
  xyz: Num3;
  occ: number;
  b: number;
  element: string;
  i_seq: number;
  is_ligand: boolean | null;
  bonds: number[];

  constructor() {
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
  }

  // http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
  from_pdb_line(pdb_line: string) {
    if (pdb_line.length < 66) {
      throw Error('ATOM or HETATM record is too short: ' + pdb_line);
    }
    const rec_type = pdb_line.substring(0, 6);
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
    const x = parseFloat(pdb_line.substring(30, 38));
    const y = parseFloat(pdb_line.substring(38, 46));
    const z = parseFloat(pdb_line.substring(46, 54));
    this.xyz = [x, y, z];
    this.occ = parseFloat(pdb_line.substring(54, 60));
    this.b = parseFloat(pdb_line.substring(60, 66));
    if (pdb_line.length >= 78) {
      this.element = pdb_line.substring(76, 78).trim().toUpperCase();
    }
    //if (pdb_line.length >= 80) {
    //  this.charge = pdb_line.substring(78, 80).trim();
    //}
    this.is_ligand = (NOT_LIGANDS.indexOf(this.resname) === -1);
  }

  b_as_u() {
    // B = 8 * pi^2 * u^2
    return Math.sqrt(this.b / (8 * 3.14159 * 3.14159));
  }

  distance_sq(other: Atom) {
    const dx = this.xyz[0] - other.xyz[0];
    const dy = this.xyz[1] - other.xyz[1];
    const dz = this.xyz[2] - other.xyz[2];
    return dx*dx + dy*dy + dz*dz;
  }

  distance(other: Atom) {
    return Math.sqrt(this.distance_sq(other));
  }

  midpoint(other: Atom) {
    return [(this.xyz[0] + other.xyz[0]) / 2,
            (this.xyz[1] + other.xyz[1]) / 2,
            (this.xyz[2] + other.xyz[2]) / 2];
  }

  is_hydrogen() {
    return this.element === 'H' || this.element === 'D';
  }

  is_ion() {
    return this.element === this.resname;
  }

  is_water() {
    return this.resname === 'HOH';
  }

  is_same_conformer(other: Atom) {
    return this.altloc === '' || other.altloc === '' ||
           this.altloc === other.altloc;
  }

  is_main_conformer() {
    return this.altloc === '' || this.altloc === 'A';
  }

  bond_radius() { // rather crude
    if (this.element === 'H') return 1.3;
    if (this.element === 'S' || this.element === 'P') return 2.43;
    return 1.99;
  }

  is_bonded_to(other: Atom) {
    const MAX_DIST = 2.2 * 2.2;
    if (!this.is_same_conformer(other)) return false;
    const dxyz2 = this.distance_sq(other);
    if (dxyz2 > MAX_DIST) return false;
    if (this.element === 'H' && other.element === 'H') return false;
    return dxyz2 <= this.bond_radius() * other.bond_radius();
  }

  resid() {
    return this.resseq + '/' + this.chain;
  }

  long_label() {
    const a = this;  // eslint-disable-line @typescript-eslint/no-this-alias
    return a.name + ' /' + a.resseq + ' ' + a.resname + '/' + a.chain +
           ' - occ: ' + a.occ.toFixed(2) + ' bf: ' + a.b.toFixed(2) +
           ' ele: ' + a.element + ' pos: (' + a.xyz[0].toFixed(2) + ',' +
           a.xyz[1].toFixed(2) + ',' + a.xyz[2].toFixed(2) + ')';
  }

  short_label() {
    const a = this;  // eslint-disable-line @typescript-eslint/no-this-alias
    return a.name + ' /' + a.resseq + ' ' + a.resname + '/' + a.chain;
  }
}


// Partition atoms into boxes for quick neighbor searching.
class Cubicles {
  boxes: number[][];
  box_length: number;
  lower_bound: Num3;
  upper_bound: Num3;
  xdim: number;
  ydim: number;
  zdim: number;

  constructor(atoms: Atom[], box_length: number,
              lower_bound: Num3, upper_bound: Num3) {
    this.boxes = [];
    this.box_length = box_length;
    this.lower_bound = lower_bound;
    this.upper_bound = upper_bound;
    this.xdim = Math.ceil((upper_bound[0] - lower_bound[0]) / box_length);
    this.ydim = Math.ceil((upper_bound[1] - lower_bound[1]) / box_length);
    this.zdim = Math.ceil((upper_bound[2] - lower_bound[2]) / box_length);
    //console.log("Cubicles: " + this.xdim + "x" + this.ydim + "x" + this.zdim);
    const nxyz = this.xdim * this.ydim * this.zdim;
    for (let j = 0; j < nxyz; j++) {
      this.boxes.push([]);
    }
    for (let i = 0; i < atoms.length; i++) {
      const xyz = atoms[i].xyz;
      const box_id = this.find_box_id(xyz[0], xyz[1], xyz[2]);
      if (box_id === null) {
        throw Error('wrong cubicle');
      }
      this.boxes[box_id].push(i);
    }
  }

  find_box_id(x: number, y: number, z: number) {
    const xstep = Math.floor((x - this.lower_bound[0]) / this.box_length);
    const ystep = Math.floor((y - this.lower_bound[1]) / this.box_length);
    const zstep = Math.floor((z - this.lower_bound[2]) / this.box_length);
    const box_id = (zstep * this.ydim + ystep) * this.xdim + xstep;
    if (box_id < 0 || box_id >= this.boxes.length) throw Error('Ups!');
    return box_id;
  }

  get_nearby_atoms(box_id: number) {
    const indices = [];
    const xydim = this.xdim * this.ydim;
    const uv = Math.max(box_id % xydim, 0);
    const u = Math.max(uv % this.xdim, 0);
    const v = Math.floor(uv / this.xdim);
    const w = Math.floor(box_id / xydim);
    console.assert((w * xydim) + (v * this.xdim) + u === box_id);
    for (let iu = u-1; iu <= u+1; iu++) {
      if (iu < 0 || iu >= this.xdim) continue;
      for (let iv = v-1; iv <= v+1; iv++) {
        if (iv < 0 || iv >= this.ydim) continue;
        for (let iw = w-1; iw <= w+1; iw++) {
          if (iw < 0 || iw >= this.zdim) continue;
          const other_box_id = (iw * xydim) + (iv * this.xdim) + iu;
          if (other_box_id >= this.boxes.length || other_box_id < 0) {
            throw Error('Box out of bounds: ID ' + other_box_id);
          }
          const box = this.boxes[other_box_id];
          for (let i = 0; i < box.length; i++) {
            indices.push(box[i]);
          }
        }
      }
    }
    return indices;
  }
}

export type { Atom };
