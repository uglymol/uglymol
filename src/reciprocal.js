// @flow
import { UnitCell } from './unitcell.js';
import { ElMap } from './elmap.js';
import { Viewer } from './viewer.js';
import { addXyzCross, makeLineMaterial, makeLineSegments } from './lines.js';
import * as THREE from 'three';


// options handled by Viewer#select_next()
const SPOT_SEL = ['all', 'unindexed', '#1']; //extended when needed
const SHOW_AXES = ['three', 'two', 'none'];

// modified ElMap for handling output of dials.rs_mapper
class ReciprocalSpaceMap extends ElMap {
  constructor(buf /*:ArrayBuffer*/) {
    super();
    this.from_ccp4(buf, false);
    if (this.unit_cell == null) return;
    // unit of the map from dials.rs_mapper is (100A)^-1, we scale it to A^-1
    const par = this.unit_cell.parameters;
    const scale = 0.01;
    this.unit_cell = new UnitCell(scale*par[0], scale*par[1], scale*par[2],
                                  par[3], par[4], par[5]); // always ==90
    // the map needs to be shifted, so it's centered at 0,0,0
    // at last, avoid translational symmetry
  }
}

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
  let pos = new Float32Array(3*n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < 3; j++) {
      pos[3*i+j] = d.rlp[i][j];
    }
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
varying vec3 vcolor;
void main() {
  // not sure how reliable is such rounding of points
  vec2 diff = gl_PointCoord - vec2(0.5, 0.5);
  float dist_sq = 4.0 * dot(diff, diff);
  if (dist_sq >= 1.0) discard;
  gl_FragColor = vec4(vcolor, 1.0 - dist_sq * dist_sq * dist_sq);
}`;


export class ReciprocalViewer extends Viewer {
  /*::
  axes: Object | null
  points: THREE.Points | null
  max_dist: number
  d_min: ?number
  d_max_inv: number
  data: Object
  */
  constructor(options /*: {[key: string]: any}*/) {
    super(options);
    this.default_camera_pos = [100, 0, 0];
    this.axes = null;
    this.points = null;
    this.max_dist = -1;
    this.d_min = null;
    this.d_max_inv = 0;
    this.data = {};
    this.config.show_only = SPOT_SEL[0];
    this.config.show_axes = SHOW_AXES[0];
    this.set_reciprocal_key_bindings();
    this.set_dropzone();
  }

  set_reciprocal_key_bindings() {
    let kb = this.key_bindings;
    // a
    kb[65] = function (evt) {
      this.select_next('axes', 'show_axes', SHOW_AXES, evt.shiftKey);
      this.set_axes();
    };
    // p
    kb[80] = function (evt) { this.permalink(); };
    // v
    kb[86] = function (evt) {
      this.select_next('show', 'show_only', SPOT_SEL, evt.shiftKey);
      const idx = SPOT_SEL.indexOf(this.config.show_only);
      this.points.material.uniforms.show_only.value = idx - 2;
    };
    // x
    kb[88] = function (evt) {
      evt.ctrlKey ? this.change_map_line(0.1) : this.change_point_size(0.5);
    };
    // z
    kb[90] = function (evt) {
      evt.ctrlKey ? this.change_map_line(-0.1) : this.change_point_size(-0.5);
    };
    // 3, numpad 3
    kb[51] = kb[99] = function () { this.shift_clip(0.1); };
    // numpad period (Linux), decimal point (Mac)
    kb[108] = kb[110] = function () { this.shift_clip(-0.1); };
    // <-
    kb[37] = function () { this.change_dmin(0.05); };
    // ->
    kb[39] = function () { this.change_dmin(-0.05); };
    // up arrow
    kb[38] = function () { this.change_dmax(0.025); };
    // down arrow
    kb[40] = function () { this.change_dmax(-0.025); };
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

  load_data(url/*:string*/, options /*:Object*/ = {}) {
    let self = this;
    this.load_file(url, {binary: false, progress: true}, function (req) {
      self.load_from_string(req.responseText, options);
      if (options.callback) options.callback();
    });
  }

  load_from_string(text/*:string*/, options /*:Object*/) {
    if (text[0] === '{') {
      this.data = parse_json(text);
    } else {
      this.data = parse_csv(text);
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
  }

  load_map_from_ab(buffer/*:ArrayBuffer*/) {
    if (this.map_bags.length > 0) {
      this.clear_el_objects(this.map_bags.pop());
    }
    let map = new ReciprocalSpaceMap(buffer);
    if (map == null || map.unit_cell == null) return;
    this.config.map_radius = map.unit_cell.parameters[0] / 2.;
    this.add_map(map, false);
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
    if (data.lattice_ids == null || data.pos == null) return;
    let color_arr = new Float32Array(3 * data.lattice_ids.length);
    this.colorize_by_id(color_arr, data.lattice_ids);
    let geometry = new THREE.BufferGeometry();
    geometry.addAttribute('position', new THREE.BufferAttribute(data.pos, 3));
    geometry.addAttribute('color', new THREE.BufferAttribute(color_arr, 3));
    const groups = new Float32Array(data.lattice_ids);
    geometry.addAttribute('group', new THREE.BufferAttribute(groups, 1));
    let material = new THREE.ShaderMaterial({
      uniforms: {
        size: { value: 3 },
        show_only: { value: -2 },
        r2_max: { value: 100 },
        r2_min: { value: 0 },
      },
      vertexShader: point_vert,
      fragmentShader: point_frag,
      vertexColors: THREE.VertexColors,
    });
    material.transparent = true;
    this.points = new THREE.Points(geometry, material);
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

  redraw_center() {}

  mousewheel_action(delta/*:number*/, evt/*:Event*/) {
    this.change_zoom_by_factor(1 + 0.0005 * delta);
  }

  change_point_size(delta/*:number*/) {
    if (this.points === null) return;
    let size = this.points.material.uniforms.size;
    size.value = Math.max(size.value + delta, 0.5);
    this.hud('point size: ' + size.value.toFixed(1));
  }

  change_dmin(delta/*:number*/) {
    if (this.d_min == null) return;
    this.d_min = Math.max(this.d_min + delta, 0.1);
    const dmax = this.d_max_inv > 0 ? 1 / this.d_max_inv : null;
    if (dmax !== null && this.d_min > dmax) this.d_min = dmax;
    if (this.points === null) return;
    this.points.material.uniforms.r2_max.value = 1 / (this.d_min * this.d_min);
    const low_res = dmax !== null ? dmax.toFixed(2) : '∞';
    this.hud('res. limit: ' + low_res + ' - ' + this.d_min.toFixed(2) + 'Å');
  }

  change_dmax(delta/*:number*/) {
    if (this.d_min == null) return;
    let v = Math.min(this.d_max_inv + delta, 1 / this.d_min);
    if (v < 1e-6) v = 0;
    this.d_max_inv = v;
    if (this.points === null) return;
    this.points.material.uniforms.r2_min.value = v * v;
    const low_res = v > 0 ? (1 / v).toFixed(2) : '∞';
    this.hud('res. limit: ' + low_res + ' - ' + this.d_min.toFixed(2) + 'Å');
  }

  redraw_models() {
    this.set_points(this.data);
  }
}

ReciprocalViewer.prototype.ColorSchemes = [
  {
    name: 'solarized dark',
    bg: 0x002b36,
    fg: 0xfdf6e3,
    map_den: 0xeee8d5,
    lattices: [0xdc322f, 0x2aa198, 0x268bd2, 0x859900,
               0xd33682, 0xb58900, 0x6c71c4, 0xcb4b16],
    axes: [0xffaaaa, 0xaaffaa, 0xaaaaff],
  },
  {
    name: 'solarized light',
    bg: 0xfdf6e3,
    fg: 0x002b36,
    map_den: 0x073642,
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
  'B = bg color',
  'M/N = zoom',
  'D/F = clip width',
  'R = center view',
  'Z/X = point size',
  'Shift+P = permalink',
  'Shift+F = full screen',
  '←/→ = max resol.',
  '↑/↓ = min resol.',
].join('\n');

ReciprocalViewer.prototype.MOUSE_HELP = Viewer.prototype.MOUSE_HELP
                                        .split('\n').slice(0, -2).join('\n');

