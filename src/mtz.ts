import { UnitCell } from './unitcell';
import { ElMap, GridArray } from './elmap';
import type { Viewer } from './viewer';
import type { Module as MtzModule, Mtz as WasmMtz } from './wasm/mtz.d.ts';

function log_timing(t0: number, text: string) {
  console.log(text + ': ' + (performance.now() - t0).toFixed(2) + ' ms.');
}

function add_map_from_mtz(viewer, mtz, map_data, is_diff: boolean) {
  const map = new ElMap();
  const mc = mtz.cell;
  map.unit_cell = new UnitCell(mc.a, mc.b, mc.c, mc.alpha, mc.beta, mc.gamma);
  map.stats.rms = mtz.rmsd;
  map.grid = new GridArray([mtz.nx, mtz.ny, mtz.nz]);
  map.grid.values.set(map_data);
  viewer.add_map(map, is_diff);
}

export
function load_maps_from_mtz_buffer(viewer: Viewer, mtz: WasmMtz,
                                   labels?: string[]) {
  if (labels != null) {
    for (let n = 0; n < labels.length; n += 2) {
      if (labels[n] === '') continue;
      const t0 = performance.now();
      const map_data = mtz.calculate_map_from_labels(labels[n], labels[n+1]);
      log_timing(t0, 'map ' + mtz.nx + 'x' + mtz.ny + 'x' + mtz.nz +
                     ' calculated in');
      if (map_data == null) {
        viewer.hud(mtz.last_error, 'ERR');
        continue;
      }
      const is_diff = (n % 4 == 2);
      add_map_from_mtz(viewer, mtz, map_data, is_diff);
    }
  } else {  // use default labels
    for (let nmap = 0; nmap < 2; ++nmap) {
      const is_diff = (nmap == 1);
      const t0 = performance.now();
      const map_data = mtz.calculate_map(is_diff);
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
function load_maps_from_mtz(gemmi: MtzModule, viewer: Viewer, url: string,
                            labels?: string[], callback?: () => void) {
  viewer.load_file(url, {binary: true, progress: true}, function (req) {
    const t0 = performance.now();
    try {
      const mtz = gemmi.readMtz(req.response);
      //console.log("[after readMTZ] wasm mem:", gemmi.HEAPU8.length / 1024, "kb");
      load_maps_from_mtz_buffer(viewer, mtz, labels);
    } catch (e) {
      viewer.hud(e.message, 'ERR');
      return;
    }
    log_timing(t0, 'load_maps_from_mtz');
    //console.log("wasm mem:", gemmi.HEAPU8.length / 1024, "kb");
    if (callback) callback();
  });
}

export
function set_pdb_and_mtz_dropzone(gemmi: MtzModule, viewer: Viewer,
                                  zone: HTMLElement) {
  viewer.set_dropzone(zone, function (file) {
    if (/\.mtz$/.test(file.name)) {
      const reader = new FileReader();
      reader.onloadend = function (evt) {
        if (evt.target.readyState == 2) {
          const t0 = performance.now();
          try {
            const mtz = gemmi.readMtz(evt.target.result);
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
