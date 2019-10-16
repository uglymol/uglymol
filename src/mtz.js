import { UnitCell } from './unitcell.js';
import { ElMap, GridArray } from './elmap.js';

/*::
 import type {Viewer} from './viewer.js'
 */

export
function load_maps_from_mtz_buffer(viewer/*:Viewer*/, mtz_buf/*:ArrayBuffer*/) {
  let t0 = performance.now();
  /* global Module, HEAPF32 */
  let arr = new Uint8Array(mtz_buf);
  let buffer = Module._malloc(arr.length);
  Module.writeArrayToMemory(arr, buffer);
  let mtz = new Module.MtzMap(buffer, arr.length);
  let t1 = performance.now();
  let t2 = [];
  let t3 = [];
  for (let nmap = 0; nmap < 2; ++nmap) {
    let is_diff = (nmap == 1);
    let map_data = mtz.calculate_map(is_diff);
    t2.push(performance.now());
    let map = new ElMap();
    map.unit_cell = new UnitCell(
      mtz.cell_param(0), mtz.cell_param(1), mtz.cell_param(2),
      mtz.cell_param(3), mtz.cell_param(4), mtz.cell_param(5));
    map.stats.rms = mtz.rmsd;
    map.grid = new GridArray([mtz.nz, mtz.ny, mtz.nx]);
    let len = mtz.nx * mtz.ny * mtz.nz;
    console.log('fft size', mtz.nx, mtz.ny, mtz.nz);
    map.grid.values.set(HEAPF32.subarray(map_data/4, map_data/4 + len));
    viewer.add_map(map, is_diff);
    t3.push(performance.now());
  }
  Module._free(buffer);
  mtz.delete();
  let t4 = performance.now();
  console.log('reading mtz: ' + (t1 - t0) + ' ms.');
  console.log('map 1 fft: ' + (t2[0] - t1) + ' ms.');
  console.log('map 1 copy: ' + (t3[0] - t2[0]) + ' ms.');
  console.log('map 2 fft: ' + (t2[1] - t3[0]) + ' ms.');
  console.log('map 2 copy: ' + (t3[1] - t2[1]) + ' ms.');
  console.log('total: ' + (t4 - t0) + ' ms.');
}

export
function load_maps_from_mtz(viewer/*:Viewer*/, url/*:string*/,
                            callback/*:?Function*/) {
  viewer.load_file(url, {binary: true, progress: true}, function (req) {
    load_maps_from_mtz_buffer(viewer, req.response);
    if (callback) callback();
  });
}

