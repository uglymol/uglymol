// @flow

import { UnitCell } from './unitcell.js';
import { Block } from './isosurface.js';

function modulo(a, b) {
  const reminder = a % b;
  return reminder >= 0 ? reminder : reminder + b;
}

class GridArray {
  /*::
  dim: number[]
  values: Float32Array
  */
  constructor(dim /*:number[]*/) {
    this.dim = dim; // dimensions of the grid for the entire unit cell
    this.values = new Float32Array(dim[0] * dim[1] * dim[2]);
  }

  grid2index(i/*:number*/, j/*:number*/, k/*:number*/) {
    i = modulo(i, this.dim[0]);
    j = modulo(j, this.dim[1]);
    k = modulo(k, this.dim[2]);
    return this.dim[2] * (this.dim[1] * i + j) + k;
  }

  grid2index_unchecked(i/*:number*/, j/*:number*/, k/*:number*/) {
    return this.dim[2] * (this.dim[1] * i + j) + k;
  }

  grid2frac(i/*:number*/, j/*:number*/, k/*:number*/) {
    return [i / this.dim[0], j / this.dim[1], k / this.dim[2]];
  }

  // return grid coordinates (rounded down) for the given fractional coordinates
  frac2grid(xyz/*:number[]*/) {
    // at one point "| 0" here made extract_block() 40% faster on V8 3.14,
    // but I don't see any effect now
    return [Math.floor(xyz[0] * this.dim[0]) | 0,
            Math.floor(xyz[1] * this.dim[1]) | 0,
            Math.floor(xyz[2] * this.dim[2]) | 0];
  }

  set_grid_value(i/*:number*/, j/*:number*/, k/*:number*/, value/*:number*/) {
    const idx = this.grid2index(i, j, k);
    this.values[idx] = value;
  }

  get_grid_value(i/*:number*/, j/*:number*/, k/*:number*/) {
    const idx = this.grid2index(i, j, k);
    return this.values[idx];
  }
}

function calculate_stddev(a, offset) {
  let sum = 0;
  let sq_sum = 0;
  const alen = a.length;
  for (let i = offset; i < alen; i++) {
    sum += a[i];
    sq_sum += a[i] * a[i];
  }
  const mean = sum / (alen - offset);
  const variance = sq_sum / (alen - offset) - mean * mean;
  return {mean: mean, rms: Math.sqrt(variance)};
}

export class ElMap {
  /*::
  unit_cell: ?UnitCell
  grid: ?GridArray
  stats: { mean: number, rms: number }
  block: Block
  unit: string
  */
  constructor() {
    this.unit_cell = null;
    this.grid = null;
    this.stats = { mean: 0.0, rms: 1.0 };
    this.block = new Block();
  }

  abs_level(sigma /*:number*/) {
    return sigma * this.stats.rms + this.stats.mean;
  }

  // http://www.ccp4.ac.uk/html/maplib.html#description
  // eslint-disable-next-line complexity
  from_ccp4(buf /*:ArrayBuffer*/, expand_symmetry /*:?boolean*/) {
    if (expand_symmetry === undefined) expand_symmetry = true;
    if (buf.byteLength < 1024) throw Error('File shorter than 1024 bytes.');
    //console.log('buf type: ' + Object.prototype.toString.call(buf));
    // for now we assume both file and host are little endian
    const iview = new Int32Array(buf, 0, 256);
    // word 53 - character string 'MAP ' to identify file type
    if (iview[52] !== 0x2050414d) throw Error('not a CCP4 map');
    // map has 3 dimensions referred to as columns (fastest changing), rows
    // and sections (c-r-s)
    const n_crs = [iview[0], iview[1], iview[2]];
    const mode = iview[3];
    let nb;
    if (mode === 2) nb = 4;
    else if (mode === 0) nb = 1;
    else throw Error('Only Mode 2 and Mode 0 of CCP4 map is supported.');
    const start = [iview[4], iview[5], iview[6]];
    const n_grid = [iview[7], iview[8], iview[9]];
    const nsymbt = iview[23]; // size of extended header in bytes
    if (1024 + nsymbt + nb*n_crs[0]*n_crs[1]*n_crs[2] !== buf.byteLength) {
      throw Error('ccp4 file too short or too long');
    }
    const fview = new Float32Array(buf, 0, buf.byteLength / 4);
    this.unit_cell = new UnitCell(fview[10], fview[11], fview[12],
                                  fview[13], fview[14], fview[15]);
    // MAPC, MAPR, MAPS - axis corresp to cols, rows, sections (1,2,3 for X,Y,Z)
    const map_crs = [iview[16], iview[17], iview[18]];
    const ax = map_crs.indexOf(1);
    const ay = map_crs.indexOf(2);
    const az = map_crs.indexOf(3);

    const min = fview[19];
    const max = fview[20];
    //const sg_number = iview[22];
    //const lskflg = iview[24];
    const grid = new GridArray(n_grid);
    if (nsymbt % 4 !== 0) {
      throw Error('CCP4 map with NSYMBT not divisible by 4 is not supported.');
    }
    let data_view;
    if (mode === 2) data_view = fview;
    else /* mode === 0 */ data_view = new Int8Array(buf);
    let idx = (1024 + nsymbt) / nb | 0;

    // We assume that if DMEAN and RMS from the header are not clearly wrong
    // they are what the user wants. Because the map can cover a small part
    // of the asu and its rmsd may be different than the total rmsd.
    this.stats.mean = fview[21];
    this.stats.rms = fview[54];
    if (this.stats.mean < min || this.stats.mean > max || this.stats.rms <= 0) {
      this.stats = calculate_stddev(data_view, idx);
    }
    let b1 = 1;
    let b0 = 0;
    // if the file was converted by mapmode2to0 - scale the data
    if (mode === 0 && iview[39] === -128 && iview[40] === 127) {
      // scaling f(x)=b1*x+b0 such that f(-128)=min and f(127)=max
      b1 = (max - min) / 255.0;
      b0 = 0.5 * (min + max + b1);
    }

    const end = [start[0] + n_crs[0], start[1] + n_crs[1], start[2] + n_crs[2]];
    let it = [0, 0, 0];
    for (it[2] = start[2]; it[2] < end[2]; it[2]++) { // sections
      for (it[1] = start[1]; it[1] < end[1]; it[1]++) { // rows
        for (it[0] = start[0]; it[0] < end[0]; it[0]++) { // cols
          grid.set_grid_value(it[ax], it[ay], it[az], b1 * data_view[idx] + b0);
          idx++;
        }
      }
    }
    if (expand_symmetry && nsymbt > 0) {
      const u8view = new Uint8Array(buf);
      for (let i = 0; i+80 <= nsymbt; i += 80) {
        let j;
        let symop = '';
        for (j = 0; j < 80; ++j) {
          symop += String.fromCharCode(u8view[1024 + i + j]);
        }
        if (/^\s*x\s*,\s*y\s*,\s*z\s*$/i.test(symop)) continue;  // skip x,y,z
        //console.log('sym ops', symop.trim());
        let mat = parse_symop(symop);
        // Note: we apply here symops to grid points instead of coordinates.
        // In the cases we came across it is equivalent, but in general not.
        for (j = 0; j < 3; ++j) {
          mat[j][3] = Math.round(mat[j][3] * n_grid[j]) | 0;
        }
        idx = (1024 + nsymbt) / nb | 0;
        let xyz = [0, 0, 0];
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
  }

  // DSN6 MAP FORMAT
  // http://www.uoxray.uoregon.edu/tnt/manual/node104.html
  // Density values are stored as bytes.
  from_dsn6(buf /*: ArrayBuffer*/) {
    //console.log('buf type: ' + Object.prototype.toString.call(buf));
    const u8data = new Uint8Array(buf);
    const iview = new Int16Array(u8data.buffer);
    if (iview[18] !== 100) {
      const len = iview.length;  // or only header, 256?
      for (let n = 0; n < len; n++) {
        // swapping bytes with Uint8Array like this:
        // var tmp=u8data[n*2]; u8data[n*2]=u8data[n*2+1]; u8data[n*2+1]=tmp;
        // was slowing down this whole function 5x times (!?) on V8.
        const val = iview[n];
        iview[n] = ((val & 0xff) << 8) | ((val >> 8) & 0xff);
      }
    }
    if (iview[18] !== 100) {
      throw Error('Endian swap failed');
    }
    const origin = [iview[0], iview[1], iview[2]];
    const n_real = [iview[3], iview[4], iview[5]];
    const n_grid = [iview[6], iview[7], iview[8]];
    const cell_mult = 1.0 / iview[17];
    this.unit_cell = new UnitCell(cell_mult * iview[9],
                                  cell_mult * iview[10],
                                  cell_mult * iview[11],
                                  cell_mult * iview[12],
                                  cell_mult * iview[13],
                                  cell_mult * iview[14]);
    const grid = new GridArray(n_grid);
    const prod = iview[15] / 100;
    const plus = iview[16];
    //var data_scale_factor = iview[15] / iview[18] + iview[16];
    // bricks have 512 (8x8x8) values
    let offset = 512;
    const n_blocks = [Math.ceil(n_real[0] / 8),
                      Math.ceil(n_real[1] / 8),
                      Math.ceil(n_real[2] / 8)];
    for (let zz = 0; zz < n_blocks[2]; zz++) {
      for (let yy = 0; yy < n_blocks[1]; yy++) {
        for (let xx = 0; xx < n_blocks[0]; xx++) { // loop over bricks
          for (let k = 0; k < 8; k++) {
            const z = 8 * zz + k;
            for (let j = 0; j < 8; j++) {
              const y = 8 * yy + j;
              for (let i = 0; i < 8; i++) { // loop inside brick
                const x = 8 * xx + i;
                if (x < n_real[0] && y < n_real[1] && z < n_real[2]) {
                  const density = (u8data[offset] - plus) / prod;
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
  }

  show_debug_info() {
    console.log('unit cell:', this.unit_cell && this.unit_cell.parameters);
    console.log('grid:', this.grid && this.grid.dim);
  }

  // Extract a block of density for calculating an isosurface using the
  // separate marching cubes implementation.
  extract_block(radius/*:number*/, center /*:[number,number,number]*/) {
    const grid = this.grid;
    const unit_cell = this.unit_cell;
    if (grid == null || unit_cell == null) return;
    const fc = unit_cell.fractionalize(center);
    const r = [radius / unit_cell.parameters[0],
               radius / unit_cell.parameters[1],
               radius / unit_cell.parameters[2]];
    const grid_min = grid.frac2grid([fc[0] - r[0], fc[1] - r[1], fc[2] - r[2]]);
    const grid_max = grid.frac2grid([fc[0] + r[0], fc[1] + r[1], fc[2] + r[2]]);
    const size = [grid_max[0] - grid_min[0] + 1,
                  grid_max[1] - grid_min[1] + 1,
                  grid_max[2] - grid_min[2] + 1];
    let points = [];
    let values = [];
    for (let i = grid_min[0]; i <= grid_max[0]; i++) {
      for (let j = grid_min[1]; j <= grid_max[1]; j++) {
        for (let k = grid_min[2]; k <= grid_max[2]; k++) {
          const frac = grid.grid2frac(i, j, k);
          const orth = unit_cell.orthogonalize(frac);
          points.push(orth);
          const map_value = grid.get_grid_value(i, j, k);
          values.push(map_value);
        }
      }
    }
    this.block.set(points, values, size);
  }

  isomesh_in_block(sigma/*:number*/, method/*:string*/) {
    const abs_level = this.abs_level(sigma);
    return this.block.isosurface(abs_level, method);
  }
}

ElMap.prototype.unit = 'e/\u212B\u00B3';

// symop -> matrix ([x,y,z] = matrix * [x,y,z,1])
function parse_symop(symop) {
  const ops = symop.toLowerCase().replace(/\s+/g, '').split(',');
  if (ops.length !== 3) throw Error('Unexpected symop: ' + symop);
  let mat = [];
  for (let i = 0; i < 3; i++) {
    const terms = ops[i].split(/(?=[+-])/);
    let row = [0, 0, 0, 0];
    for (let j = 0; j < terms.length; j++) {
      const term = terms[j];
      const sign = (term[0] === '-' ? -1 : 1);
      let m = terms[j].match(/^[+-]?([xyz])$/);
      if (m) {
        const pos = {x: 0, y: 1, z: 2}[m[1]];
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

