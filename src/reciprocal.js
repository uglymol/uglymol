// @flow
import { Viewer } from './viewer.js';
import { addXyzCross, makeLineMaterial, makeLineSegments } from './lines.js';
import * as THREE from 'three';


var SPOT_SEL = ['all', 'indexed', 'not indexed'];
var SHOW_AXES = ['three', 'two', 'none'];

export function ReciprocalViewer(options /*: {[key: string]: any}*/) {
  Viewer.call(this, options);
  this.default_camera_pos = [100, 0, 0];
  this.axes = null;
  this.points = null;
  this.max_dist = null;
  this.d_min = null;
  this.d_max_inv = 0;
  this.data = null;
  this.config.show_only = SPOT_SEL[0];
  this.config.show_axes = SHOW_AXES[0];
  this.set_reciprocal_key_bindings();
}

ReciprocalViewer.prototype = Object.create(Viewer.prototype);
ReciprocalViewer.prototype.constructor = ReciprocalViewer;

ReciprocalViewer.prototype.KEYBOARD_HELP = [
  '<b>keyboard:</b>',
  'H = toggle help',
  'V = show (un)indexed',
  'A = toggle axes',
  'B = bg color',
  'M/N = zoom',
  'D/F = clip width',
  'R = center view',
  'Home/End = point size',
  'Shift+P = permalink',
  'Shift+F = full screen',
  '←/→ = max resol.',
  '↑/↓ = min resol.',
].join('\n');

ReciprocalViewer.prototype.MOUSE_HELP = Viewer.prototype.MOUSE_HELP
                                        .split('\n').slice(0, -2).join('\n');

ReciprocalViewer.prototype.set_reciprocal_key_bindings = function () {
  var kb = this.key_bindings;
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
    var show_only = this.points.material.uniforms.show_only;
    var sel_map = { 'all': -2, 'indexed': 0, 'not indexed': -1 };
    show_only.value = sel_map[this.config.show_only];
  };
  // Home
  kb[36] = function (evt) {
    evt.ctrlKey ? this.change_map_line(0.1) : this.change_point_size(0.5);
  };
  // End
  kb[35] = function (evt) {
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
};

ReciprocalViewer.prototype.load_data = function (url, options) {
  options = options || {};
  var self = this;
  this.load_file(url, {binary: false, progress: true}, function (req) {
    self.parse_data(req.responseText);
    self.set_axes();
    self.set_points();
    self.camera.zoom = 0.5 * (self.camera.top - self.camera.bottom);
    // default scale is set to 100 - same as default_camera_pos
    var d = 1.01 * self.max_dist;
    self.controls.slab_width = [d, d, 100];
    self.set_view(options);
    if (options.callback) options.callback();
  });
};

ReciprocalViewer.prototype.parse_data = function (text) {
  var lines = text.split('\n').filter(function (line) {
    return line.length > 0 && line[0] !== '#';
  });
  var pos = new Float32Array(lines.length * 3);
  var lattice_ids = [];
  var max_sq = 0;
  for (var i = 0; i < lines.length; i++) {
    var nums = lines[i].split(',').map(Number);
    var sq = nums[0]*nums[0] + nums[1]*nums[1] + nums[2]*nums[2];
    if (sq > max_sq) max_sq = sq;
    for (var j = 0; j < 3; j++) {
      pos[3*i+j] = nums[j];
    }
    lattice_ids.push(nums[3]);
  }
  this.max_dist = Math.sqrt(max_sq);
  this.d_min = 1 / this.max_dist;
  this.data = { pos: pos, lattice_ids: lattice_ids };
};

ReciprocalViewer.prototype.set_axes = function () {
  if (this.axes != null) {
    this.remove_and_dispose(this.axes);
    this.axes = null;
  }
  if (this.config.show_axes === 'none') return;
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
    segments: true,
  });
  this.axes = makeLineSegments(material, vertices, colors);
  this.scene.add(this.axes);
};

var point_vert = [
  'attribute float group;',
  'uniform float show_only;',
  'uniform float r2_max;',
  'uniform float r2_min;',
  'uniform float size;',
  'varying vec3 vcolor;',
  'varying float vsel;',
  'void main() {',
  '  vcolor = color;',
  '  bool sel = show_only == -2.0 || show_only == group;',
  '  float r2 = dot(position, position);',
  '  vsel = sel && r2_min <= r2 && r2 < r2_max ? 1.0 : 0.0;',
  '  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);',
  '  gl_PointSize = size;',
  '}'].join('\n');

var point_frag = [
  'varying vec3 vcolor;',
  'varying float vsel;',
  'void main() {',
  // not sure how reliable is such rounding of points
  '  vec2 diff = gl_PointCoord - vec2(0.5, 0.5);',
  '  float dist_sq = 4.0 * dot(diff, diff);',
  '  if (vsel == 0.0 || dist_sq >= 1.0) discard;',
  '  gl_FragColor = vec4(vcolor, 1.0 - dist_sq * dist_sq * dist_sq);',
  '}'].join('\n');


ReciprocalViewer.prototype.set_points = function () {
  if (this.data == null) return;
  var pos = this.data.pos;
  var lattice_ids = this.data.lattice_ids;
  var color_arr = new Float32Array(3 * lattice_ids.length);
  this.colorize_by_id(color_arr, lattice_ids);
  var geometry = new THREE.BufferGeometry();
  geometry.addAttribute('position', new THREE.BufferAttribute(pos, 3));
  geometry.addAttribute('color', new THREE.BufferAttribute(color_arr, 3));
  var groups = new Float32Array(lattice_ids);
  geometry.addAttribute('group', new THREE.BufferAttribute(groups, 1));
  var material = new THREE.ShaderMaterial({
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
};

ReciprocalViewer.prototype.colorize_by_id = function (color_arr, group_id) {
  var palette = this.config.colors.lattices;
  for (var i = 0; i < group_id.length; i++) {
    var c = palette[(group_id[i] + 1) % 4];
    color_arr[3*i] = c.r;
    color_arr[3*i+1] = c.g;
    color_arr[3*i+2] = c.b;
  }
};

ReciprocalViewer.prototype.redraw_center = function () {};

ReciprocalViewer.prototype.mousewheel_action = function (delta, evt) {
  this.change_zoom_by_factor(1 + 0.0005 * delta);
};

ReciprocalViewer.prototype.change_point_size = function (delta) {
  if (this.points === null) return;
  var size = this.points.material.uniforms.size;
  size.value = Math.max(size.value + delta, 0.5);
  this.hud('point size: ' + size.value.toFixed(1));
};

ReciprocalViewer.prototype.change_dmin = function (delta) {
  if (this.d_min == null) return;
  this.d_min = Math.max(this.d_min + delta, 0.1);
  var dmax = this.d_max_inv > 0 ? 1 / this.d_max_inv : null;
  if (dmax !== null && this.d_min > dmax) this.d_min = dmax;
  this.points.material.uniforms.r2_max.value = 1 / (this.d_min * this.d_min);
  var low_res = dmax !== null ? dmax.toFixed(2) : '∞';
  this.hud('res. limit: ' + low_res + ' - ' + this.d_min.toFixed(2) + 'Å');
};

ReciprocalViewer.prototype.change_dmax = function (delta) {
  if (this.d_min == null) return;
  var v = Math.min(this.d_max_inv + delta, 1 / this.d_min);
  if (v < 1e-6) v = 0;
  this.d_max_inv = v;
  this.points.material.uniforms.r2_min.value = v * v;
  var low_res = v > 0 ? (1 / v).toFixed(2) : '∞';
  this.hud('res. limit: ' + low_res + ' - ' + this.d_min.toFixed(2) + 'Å');
};

ReciprocalViewer.prototype.redraw_models = function () {
  if (this.points) this.remove_and_dispose(this.points);
  this.set_points();
};

ReciprocalViewer.prototype.ColorSchemes = [
  {
    name: 'solarized dark',
    bg: 0x002b36,
    fg: 0xfdf6e3,
    lattices: [0xdc322f, 0x2aa198, 0x268bd2, 0x859900,
               0xd33682, 0xb58900, 0x6c71c4, 0xcb4b16],
    axes: [0xffaaaa, 0xaaffaa, 0xaaaaff],
  },
  {
    name: 'solarized light',
    bg: 0xfdf6e3,
    fg: 0x002b36,
    lattices: [0xdc322f, 0x2aa198, 0x268bd2, 0x859900,
               0xd33682, 0xb58900, 0x6c71c4, 0xcb4b16],
    axes: [0xffaaaa, 0xaaffaa, 0xaaaaff],
  },
];
