// @flow

import { UnitCell } from './unitcell.js';
import { ElMap, GridArray } from './elmap.js';

/*::
 import type {Viewer} from './viewer.js'
 */

function log_timing(t0/*:number*/, text/*:string*/) {
  console.log(text + ': ' + (performance.now() - t0).toFixed(2) + ' ms.');
}

function add_map_from_mtz(viewer, mtz, map_data, is_diff/*:boolean*/) {
  let map = new ElMap();
  let mc = mtz.cell;
  map.unit_cell = new UnitCell(mc.a, mc.b, mc.c, mc.alpha, mc.beta, mc.gamma);
  map.stats.rms = mtz.rmsd;
  map.grid = new GridArray([mtz.nx, mtz.ny, mtz.nz]);
  map.grid.values.set(map_data);
  viewer.add_map(map, is_diff);
}

export
function load_maps_from_mtz_buffer(viewer/*:Viewer*/, mtz/*:Object*/,
                                   labels/*:?string[]*/) {
  if (labels != null) {
    for (let n = 0; n < labels.length; n += 2) {
      if (labels[n] === '') continue;
      let t0 = performance.now();
      let map_data = mtz.calculate_map_from_labels(labels[n], labels[n+1]);
      log_timing(t0, 'map ' + mtz.nx + 'x' + mtz.ny + 'x' + mtz.nz +
                     ' calculated in');
      if (map_data == null) {
        viewer.hud(mtz.last_error, 'ERR');
        continue;
      }
      let is_diff = (n % 4 == 2);
      add_map_from_mtz(viewer, mtz, map_data, is_diff);
    }
  } else {  // use default labels
    for (let nmap = 0; nmap < 2; ++nmap) {
      let is_diff = (nmap == 1);
      let t0 = performance.now();
      let map_data = mtz.calculate_map(is_diff);
      log_timing(t0, 'map ' + mtz.nx + 'x' + mtz.ny + 'x' + mtz.nz +
                     ' calculated in');
      if (map_data != null) {
        add_map_from_mtz(viewer, mtz, map_data, is_diff);
      }
    }
  }
  mtz.delete();
}

export
function load_maps_from_mtz(Gemmi/*:Object*/, viewer/*:Viewer*/, url/*:string*/,
                            labels/*:?string[]*/, callback/*:?Function*/) {
  viewer.load_file(url, {binary: true, progress: true}, function (req) {
    let t0 = performance.now();
    let mtz = Gemmi.readMtz(req.response);
    load_maps_from_mtz_buffer(viewer, mtz, labels);
    log_timing(t0, 'load_maps_from_mtz');
    if (callback) callback();
  });
}

export
function set_pdb_and_mtz_dropzone(Gemmi/*:Object*/, viewer/*:Viewer*/,
                                  zone/*:Object*/) {
  viewer.set_dropzone(zone, function (file) {
    if (/\.mtz$/.test(file.name)) {
      const reader = new FileReader();
      reader.onloadend = function (evt/*:any*/) {
        if (evt.target.readyState == 2) {
          let t0 = performance.now();
          let mtz = Gemmi.readMtz(evt.target.result);
          load_maps_from_mtz_buffer(viewer, mtz);
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
