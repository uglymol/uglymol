import { UnitCell } from './unitcell';
import { Block } from './isosurface';
import type { Module as GemmiModule, Ccp4Map as WasmCcp4Map } from './gemmi_wasm';

type Num3 = [number, number, number];

function modulo(a: number, b: number) {
  const reminder = a % b;
  return reminder >= 0 ? reminder : reminder + b;
}

export class GridArray {
  dim: Num3;
  values: Float32Array;

  constructor(dim: Num3) {
    this.dim = dim; // dimensions of the grid for the entire unit cell
    this.values = new Float32Array(dim[0] * dim[1] * dim[2]);
  }

  grid2index(i: number, j: number, k: number) {
    i = modulo(i, this.dim[0]);
    j = modulo(j, this.dim[1]);
    k = modulo(k, this.dim[2]);
    return this.dim[2] * (this.dim[1] * i + j) + k;
  }

  grid2index_unchecked(i: number, j: number, k: number) {
    return this.dim[2] * (this.dim[1] * i + j) + k;
  }

  grid2frac(i: number, j: number, k: number): Num3 {
    return [i / this.dim[0], j / this.dim[1], k / this.dim[2]];
  }

  // return grid coordinates (rounded down) for the given fractional coordinates
  frac2grid(xyz: Num3) {
    // at one point "| 0" here made extract_block() 40% faster on V8 3.14,
    // but I don't see any effect now
    return [Math.floor(xyz[0] * this.dim[0]) | 0,
            Math.floor(xyz[1] * this.dim[1]) | 0,
            Math.floor(xyz[2] * this.dim[2]) | 0];
  }

  set_grid_value(i: number, j: number, k: number, value: number) {
    const idx = this.grid2index(i, j, k);
    this.values[idx] = value;
  }

  get_grid_value(i: number, j: number, k: number) {
    const idx = this.grid2index(i, j, k);
    return this.values[idx];
  }
}

function calculate_stddev(a: Float32Array|Int8Array, offset: number) {
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
  unit_cell: UnitCell | null;
  grid: GridArray | null;
  stats: { mean: number, rms: number };
  block: Block;
  declare unit: string;
  box_size?: Num3; // used in ReciprocalSpaceMap

  constructor() {
    this.unit_cell = null;
    this.grid = null;
    this.stats = { mean: 0.0, rms: 1.0 };
    this.block = new Block();
  }

  abs_level(sigma: number) {
    return sigma * this.stats.rms + this.stats.mean;
  }

  from_ccp4(buf: ArrayBuffer, expand_symmetry?: boolean, gemmi?: GemmiModule) {
    if (expand_symmetry === undefined) expand_symmetry = true;
    if (gemmi == null || typeof gemmi.readCcp4Map !== 'function') {
      throw Error('Gemmi is required for CCP4 map loading.');
    }
    const ccp4 = gemmi.readCcp4Map(buf, expand_symmetry);
    try {
      this.set_from_ccp4_map(ccp4);
    } finally {
      ccp4.delete();
    }
  }

  // DSN6 MAP FORMAT
  // http://www.uoxray.uoregon.edu/tnt/manual/node104.html
  // Density values are stored as bytes.
  from_dsn6(buf: ArrayBuffer) {
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
    const origin: Num3 = [iview[0], iview[1], iview[2]];
    const n_real: Num3 = [iview[3], iview[4], iview[5]];
    const n_grid: Num3 = [iview[6], iview[7], iview[8]];
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
  extract_block(radius: number, center: Num3) {
    const grid = this.grid;
    const unit_cell = this.unit_cell;
    if (grid == null || unit_cell == null) return;
    const fc = unit_cell.fractionalize(center);
    const r = [radius / unit_cell.parameters[0],
               radius / unit_cell.parameters[1],
               radius / unit_cell.parameters[2]];
    const grid_min = grid.frac2grid([fc[0] - r[0], fc[1] - r[1], fc[2] - r[2]]);
    const grid_max = grid.frac2grid([fc[0] + r[0], fc[1] + r[1], fc[2] + r[2]]);
    const size: Num3 = [grid_max[0] - grid_min[0] + 1,
                        grid_max[1] - grid_min[1] + 1,
                        grid_max[2] - grid_min[2] + 1];
    const points = [];
    const values = [];
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

  isomesh_in_block(sigma: number, method: string) {
    const abs_level = this.abs_level(sigma);
    return this.block.isosurface(abs_level, method);
  }

  private set_from_ccp4_map(ccp4: WasmCcp4Map) {
    const cell = ccp4.cell;
    this.unit_cell = new UnitCell(cell.a, cell.b, cell.c,
                                  cell.alpha, cell.beta, cell.gamma);
    this.stats.mean = ccp4.mean;
    this.stats.rms = ccp4.rms;
    const dim: Num3 = [ccp4.nx, ccp4.ny, ccp4.nz];
    const grid = new GridArray(dim);
    const values = ccp4.data();
    for (let x = 0; x < dim[0]; ++x) {
      for (let y = 0; y < dim[1]; ++y) {
        for (let z = 0; z < dim[2]; ++z) {
          const src = (z * dim[1] + y) * dim[0] + x;
          grid.values[grid.grid2index_unchecked(x, y, z)] = values[src];
        }
      }
    }
    this.grid = grid;
  }

}

ElMap.prototype.unit = 'e/\u212B\u00B3';
