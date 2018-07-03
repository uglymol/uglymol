// @flow
import { ElMap } from './elmap.js';
import { Viewer } from './viewer.js';
import { addXyzCross, makeLineMaterial, makeLineSegments,
         makeUniforms } from './lines.js';
import * as THREE from 'three';


// options handled by Viewer#select_next()
const SPOT_SEL = ['all', 'unindexed', '#1']; //extended when needed
const SHOW_AXES = ['two', 'three', 'none'];
const SPOT_SHAPES = ['wheel', 'square'];

// Modified ElMap for handling output of dials.rs_mapper.
// rs_mapper outputs map in ccp4 format, but we need to rescale it,
// shift it so the box is centered at 0,0,0,
// and the translational symmetry doesn't apply.
export class ReciprocalSpaceMap extends ElMap {
  /*::
    box_size: [number, number, number]
   */
  constructor(buf /*:ArrayBuffer*/) {
    super();
    this.box_size = [1, 1, 1];
    this.from_ccp4(buf, false);
    if (this.unit_cell == null) return;
    // unit of the map from dials.rs_mapper is (100A)^-1, we scale it to A^-1
    // We assume the "unit cell" is cubic -- as it is in rs_mapper.
    const par = this.unit_cell.parameters;
    this.box_size = [par[0]/ 100, par[1] / 100, par[2] / 100];
    this.unit_cell = null;
  }

  extract_block(radius/*:number*/, center/*:[number,number,number]*/) {
    const grid = this.grid;
    if (grid == null) return;
    const b = this.box_size;
    const lo_bounds = [];
    const hi_bounds = [];
    for (let n = 0; n < 3; n++) {
      let lo = Math.floor(grid.dim[n] * ((center[n] - radius) / b[n] + 0.5));
      let hi = Math.floor(grid.dim[n] * ((center[n] + radius) / b[n] + 0.5));
      lo = Math.min(Math.max(0, lo), grid.dim[n] - 1);
      hi = Math.min(Math.max(0, hi), grid.dim[n] - 1);
      if (lo === hi) return;
      lo_bounds.push(lo);
      hi_bounds.push(hi);
    }

    const points = [];
    const values = [];
    for (let i = lo_bounds[0]; i <= hi_bounds[0]; i++) {
      for (let j = lo_bounds[1]; j <= hi_bounds[1]; j++) {
        for (let k = lo_bounds[2]; k <= hi_bounds[2]; k++) {
          points.push([(i / grid.dim[0] - 0.5) * b[0],
                       (j / grid.dim[1] - 0.5) * b[1],
                       (k / grid.dim[2] - 0.5) * b[2]]);
          const index = grid.grid2index_unchecked(i, j, k);
          values.push(grid.values[index]);
        }
      }
    }

    const size = [hi_bounds[0] - lo_bounds[0] + 1,
                  hi_bounds[1] - lo_bounds[1] + 1,
                  hi_bounds[2] - lo_bounds[2] + 1];
    this.block.set(points, values, size);
  }
}

ReciprocalSpaceMap.prototype.unit = '';

function find_max_dist(pos) {
  let max_sq = 0;
  for (let i = 0; i < pos.length; i += 3) {
    const n = 3 * i;
    const sq = pos[n]*pos[n] + pos[n+1]*pos[n+1] + pos[n+2]*pos[n+2];
    if (sq > max_sq) max_sq = sq;
  }
  return Math.sqrt(max_sq);
}

function max_val(arr) {
  let max = -Infinity;
  for (let i = 0; i < arr.length; i++) {
    if (arr[i] > max) max = arr[i];
  }
  return max;
}

function parse_csv(text) {
  const lines = text.split('\n').filter(function (line) {
    return line.length > 0 && line[0] !== '#';
  });
  let pos = new Float32Array(lines.length * 3);
  let lattice_ids = [];
  for (let i = 0; i < lines.length; i++) {
    const nums = lines[i].split(',').map(Number);
    for (let j = 0; j < 3; j++) {
      pos[3*i+j] = nums[j];
    }
    lattice_ids.push(nums[3]);
  }
  return { pos, lattice_ids };
}

function minus_ones(n) {
  const a = [];
  for (let i = 0; i < n; i++) a.push(-1);
  return a;
}

function parse_json(text) {
  const d = JSON.parse(text);
  const n = d.rlp.length;
  let pos;
  if (n > 0 && d.rlp[0] instanceof Array) { // deprecated format
    pos = new Float32Array(3*n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < 3; j++) {
        pos[3*i+j] = d.rlp[i][j];
      }
    }
  } else { // flat array - new format
    pos = new Float32Array(d.rlp);
  }
  const lattice_ids = d.experiment_id || minus_ones(n);
  return { pos, lattice_ids };
}

const point_vert = `
attribute float group;
uniform float show_only;
uniform float r2_max;
uniform float r2_min;
uniform float size;
varying vec3 vcolor;
void main() {
  vcolor = color;
  float r2 = dot(position, position);
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  if (r2 < r2_min || r2 >= r2_max || (show_only != -2.0 && show_only != group))
    gl_Position.x = 2.0;
  gl_PointSize = size;
}`;

const point_frag = `
#include <fog_pars_fragment>
varying vec3 vcolor;
void main() {
  float alpha = 1.0;
#ifdef ROUND
  // not sure how reliable is such rounding of points
  vec2 diff = gl_PointCoord - vec2(0.5, 0.5);
  float dist_sq = 4.0 * dot(diff, diff);
  if (dist_sq >= 1.0) discard;
  alpha = 1.0 - dist_sq * dist_sq * dist_sq;
#endif
  gl_FragColor = vec4(vcolor, alpha);
#include <fog_fragment>
}`;


export class ReciprocalViewer extends Viewer {
  /*::
  axes: Object | null
  points: THREE.Points | null
  max_dist: number
  d_min: number
  d_max_inv: number
  data: Object
  point_material: THREE.ShaderMaterial
  */
  constructor(options/*:{[key:string]: any}*/) {
    super(options);
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
    this.set_dropzone();
    this.point_material = new THREE.ShaderMaterial({
      uniforms: makeUniforms({
        size: 3,
        show_only: -2,
        r2_max: 100,
        r2_min: 0,
      }),
      defines: {
        ROUND: true,
      },
      vertexShader: point_vert,
      fragmentShader: point_frag,
      vertexColors: THREE.VertexColors,
      fog: true,
      transparent: true,
    });
  }

  set_reciprocal_key_bindings() {
    let kb = this.key_bindings;
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
      this.point_material.defines.ROUND = (this.config.spot_shape === 'wheel');
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
      const idx = SPOT_SEL.indexOf(this.config.show_only);
      this.point_material.uniforms.show_only.value = idx - 2;
    };
    // x
    kb[88] = function (evt) {
      evt.ctrlKey ? this.change_map_line(0.1) : this.change_point_size(0.5);
    };
    // z
    kb[90] = function (evt) {
      evt.ctrlKey ? this.change_map_line(-0.1) : this.change_point_size(-0.5);
    };
    // comma
    kb[188] = function (evt) { if (evt.shiftKey) this.shift_clip(0.1); };
    // period
    kb[190] = function (evt) { if (evt.shiftKey) this.shift_clip(-0.1); };
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
  }

  set_dropzone() {
    if (typeof document === 'undefined') return;  // for testing on node
    const zone = this.renderer.domElement;
    const self = this;
    zone.addEventListener('dragover', function (e) {
      e.stopPropagation();
      e.preventDefault();
      e.dataTransfer.dropEffect = 'copy';
      self.hud('ready for drop...');
    });
    zone.addEventListener('drop', function (e) {
      e.stopPropagation();
      e.preventDefault();
      let names = [];
      for (const file of e.dataTransfer.files) {
        const reader = new FileReader();
        if (/\.(map|ccp4)$/.test(file.name)) {
          reader.onloadend = function (evt) {
            if (evt.target.readyState == 2) {
              self.load_map_from_ab(evt.target.result);
            }
          };
          reader.readAsArrayBuffer(file);
        } else {
          reader.onload = function (evt) {
            self.load_from_string(evt.target.result, {});
          };
          reader.readAsText(file);
        }
        names.push(file.name);
      }
      self.hud('loading ' + names.join(', '));
    });
  }

  load_data(url/*:string*/, options/*:Object*/ = {}) {
    let self = this;
    this.load_file(url, {binary: false, progress: true}, function (req) {
      let ok = self.load_from_string(req.responseText, options);
      if (ok && options.callback) options.callback();
    });
  }

  load_from_string(text/*:string*/, options/*:Object*/) {
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
    const last_group = max_val(this.data.lattice_ids);
    SPOT_SEL.splice(3);
    for (let i = 1; i <= last_group; i++) {
      SPOT_SEL.push('#' + (i + 1));
    }
    this.set_axes();
    this.set_points(this.data);
    this.camera.zoom = 0.5 * (this.camera.top - this.camera.bottom);
    // default scale is set to 100 - same as default_camera_pos
    const d = 1.01 * this.max_dist;
    this.controls.slab_width = [d, d, 100];
    this.set_view(options);
    this.hud('Loaded ' + this.data.pos.length + ' spots.');
    return true;
  }

  load_map_from_ab(buffer/*:ArrayBuffer*/) {
    if (this.map_bags.length > 0) {
      this.clear_el_objects(this.map_bags.pop());
    }
    let map = new ReciprocalSpaceMap(buffer);
    if (map == null) return;
    const map_range = map.box_size[0] / 2;
    this.config.map_radius = Math.round(map_range / 2 * 100) / 100;
    this.config.max_map_radius = Math.round(1.5 * map_range * 100) / 100;
    this.config.default_isolevel = 0.3;
    this.add_map(map, false);
    const map_dmin = 1 / map_range;
    let msg = 'Loaded density map (' + map_dmin.toFixed(2) + 'Å).\n';
    if (this.points !== null && map_dmin > this.d_min) {
      msg += 'Adjusted spot clipping. ';
      this.change_dmin(map_dmin - this.d_min);
    }
    this.hud(msg + 'Use +/- to change the isolevel.');
  }

  set_axes() {
    if (this.axes != null) {
      this.remove_and_dispose(this.axes);
      this.axes = null;
    }
    if (this.config.show_axes === 'none') return;
    const axis_length = 1.2 * this.max_dist;
    let vertices = [];
    addXyzCross(vertices, [0, 0, 0], axis_length);
    const ca = this.config.colors.axes;
    const colors = [ca[0], ca[0], ca[1], ca[1], ca[2], ca[2]];
    if (this.config.show_axes === 'two') {
      vertices.splice(4);
      colors.splice(4);
    }
    const material = makeLineMaterial({
      win_size: this.window_size,
      linewidth: 3,
      segments: true,
    });
    this.axes = makeLineSegments(material, vertices, colors);
    this.scene.add(this.axes);
  }

  set_points(data/*:Object*/) {
    if (this.points != null) {
      this.remove_and_dispose(this.points);
      this.points = null;
    }
    if (data == null || data.lattice_ids == null || data.pos == null) return;
    let color_arr = new Float32Array(3 * data.lattice_ids.length);
    this.colorize_by_id(color_arr, data.lattice_ids);
    let geometry = new THREE.BufferGeometry();
    geometry.addAttribute('position', new THREE.BufferAttribute(data.pos, 3));
    geometry.addAttribute('color', new THREE.BufferAttribute(color_arr, 3));
    const groups = new Float32Array(data.lattice_ids);
    geometry.addAttribute('group', new THREE.BufferAttribute(groups, 1));
    this.points = new THREE.Points(geometry, this.point_material);
    this.scene.add(this.points);
    this.request_render();
  }

  colorize_by_id(color_arr/*:Float32Array*/, group_id/*:number[]*/) {
    const palette = this.config.colors.lattices;
    for (let i = 0; i < group_id.length; i++) {
      const c = palette[(group_id[i] + 1) % 4];
      color_arr[3*i] = c.r;
      color_arr[3*i+1] = c.g;
      color_arr[3*i+2] = c.b;
    }
  }

  mousewheel_action(delta/*:number*/, evt/*:Event*/) {
    this.change_zoom_by_factor(1 + 0.0005 * delta);
  }

  change_point_size(delta/*:number*/) {
    let size = this.point_material.uniforms.size;
    size.value = Math.max(size.value + delta, 0.5);
    this.hud('point size: ' + size.value.toFixed(1));
  }

  change_dmin(delta/*:number*/) {
    this.d_min = Math.max(this.d_min + delta, 0.1);
    const dmax = this.d_max_inv > 0 ? 1 / this.d_max_inv : null;
    if (dmax !== null && this.d_min > dmax) this.d_min = dmax;
    this.point_material.uniforms.r2_max.value = 1 / (this.d_min * this.d_min);
    const low_res = dmax !== null ? dmax.toFixed(2) : '∞';
    this.hud('res. limit: ' + low_res + ' - ' + this.d_min.toFixed(2) + 'Å');
  }

  change_dmax(delta/*:number*/) {
    let v = Math.min(this.d_max_inv + delta, 1 / this.d_min);
    if (v < 1e-6) v = 0;
    this.d_max_inv = v;
    this.point_material.uniforms.r2_min.value = v * v;
    const low_res = v > 0 ? (1 / v).toFixed(2) : '∞';
    this.hud('res. limit: ' + low_res + ' - ' + this.d_min.toFixed(2) + 'Å');
  }

  redraw_models() {
    this.set_points(this.data);
  }

  get_cell_box_func() {
    if (this.map_bags.length === 0) return null;
    // $FlowFixMe: here the map is ReciprocalSpaceMap not ElMap
    const a = this.map_bags[0].map.box_size;
    return function (xyz/*:[number,number,number]*/) {
      return [(xyz[0]-0.5) * a[0], (xyz[1]-0.5) * a[1], (xyz[2]-0.5) * a[2]];
    };
  }
}

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
  },
];

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
  '+/- = map level',
].join('\n');

ReciprocalViewer.prototype.MOUSE_HELP =
    Viewer.prototype.MOUSE_HELP.split('\n').slice(0, -2).join('\n');

