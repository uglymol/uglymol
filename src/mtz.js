// @flow

import { UnitCell } from './unitcell.js';
import { ElMap, GridArray } from './elmap.js';

/*::
 import type {Viewer} from './viewer.js'
 declare var Module: any;
 declare var HEAPF32: Float32Array;
 */

function add_map_from_mtz(viewer, mtz, map_data, is_diff/*:boolean*/) {
  //t2.push(performance.now());
  let map = new ElMap();
  map.unit_cell = new UnitCell(
    mtz.cell_param(0), mtz.cell_param(1), mtz.cell_param(2),
    mtz.cell_param(3), mtz.cell_param(4), mtz.cell_param(5));
  map.stats.rms = mtz.rmsd;
  let len = mtz.nx * mtz.ny * mtz.nz;
  console.log('fft size', mtz.nx, mtz.ny, mtz.nz);
  map.grid = new GridArray([mtz.nz, mtz.ny, mtz.nx]);
  map.grid.values.set(HEAPF32.subarray(map_data/4, map_data/4 + len));
  viewer.add_map(map, is_diff);
  //t3.push(performance.now());
}

export
function load_maps_from_mtz_buffer(viewer/*:Viewer*/, mtz_buf/*:ArrayBuffer*/,
                                   labels/*:?string[]*/) {
  //let t0 = performance.now();
  /* global Module, HEAPF32 */
  let arr = new Uint8Array(mtz_buf);
  let buffer = Module._malloc(arr.length);
  Module.writeArrayToMemory(arr, buffer);
  let mtz = new Module.Mtz(buffer, arr.length);
  //let t1 = performance.now();
  //let t2 = [];
  //let t3 = [];
  if (labels != null) {
    for (let n = 0; n < labels.length; n += 2) {
      if (labels[n] === '') continue;
      let map_data = mtz.calculate_map_from_labels(labels[n], labels[n+1]);
      if (map_data !== 0) {
        let is_diff = (n % 4 == 2);
        add_map_from_mtz(viewer, mtz, map_data, is_diff);
      }
    }
  } else {
    for (let nmap = 0; nmap < 2; ++nmap) {
      let is_diff = (nmap == 1);
      let map_data = mtz.calculate_map(is_diff);
      if (map_data !== 0) {
        add_map_from_mtz(viewer, mtz, map_data, is_diff);
      }
    }
  }
  //Module._free(buffer);
  mtz.delete();
  //let t4 = performance.now();
  //console.log('reading mtz: ' + (t1 - t0) + ' ms.');
  //console.log('map 1 fft: ' + (t2[0] - t1) + ' ms.');
  //console.log('map 1 copy: ' + (t3[0] - t2[0]) + ' ms.');
  //console.log('map 2 fft: ' + (t2[1] - t3[0]) + ' ms.');
  //console.log('map 2 copy: ' + (t3[1] - t2[1]) + ' ms.');
  //console.log('total: ' + (t4 - t0) + ' ms.');
}

export
function load_maps_from_mtz(viewer/*:Viewer*/, url/*:string*/,
                            labels/*:?string[]*/,
                            callback/*:?Function*/) {
  viewer.load_file(url, {binary: true, progress: true}, function (req) {
    load_maps_from_mtz_buffer(viewer, req.response, labels);
    if (callback) callback();
  });
}

export
function set_pdb_and_mtz_dropzone(viewer/*:Viewer*/, zone/*:Object*/) {
  viewer.set_dropzone(zone, function (file) {
    if (/\.mtz$/.test(file.name)) {
      const reader = new FileReader();
      reader.onloadend = function (evt/*:any*/) {
        if (evt.target.readyState == 2) {
          load_maps_from_mtz_buffer(viewer, evt.target.result);
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
