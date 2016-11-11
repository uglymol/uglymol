// @flow
//import { Viewer } from './viewer.js';
//import * as THREE from 'three';
var Viewer = UM.Viewer;

var SPOT_SEL = ['all', 'indexed', 'not indexed'];
var sel_map = { all: -2, indexed: 0, 'not indexed': -1 };

//export
function ReciprocalViewer(options /*: {[key: string]: any}*/) {
  Viewer.call(this, options);
  this.points = null;
  this.config.show_only = SPOT_SEL[0];
  var self = this;
  this.custom_keydown = {
    86/*v*/: function (evt) {
      self.select_next('show', 'show_only', SPOT_SEL, evt.shiftKey);
      var show_only = self.points.material.uniforms.show_only;
      show_only.value = sel_map[self.config.show_only];
    }
  };
}

ReciprocalViewer.prototype = Object.create(Viewer.prototype);
ReciprocalViewer.prototype.constructor = ReciprocalViewer;

ReciprocalViewer.prototype.load_data = function (url, options) {
  options = options || {};

  var self = this;
  this.load_file(url, false, function (req) {
    var lines = req.responseText.split('\n').filter(function (line) {
      return line.length > 0 && line[0] !== '#';
    });
    var pos = new Float32Array(lines.length * 3);
    var experiment_ids = [];
    var n_col = 5;
    var bounds = [];
    var i;
    for (i = 0; i < n_col; i++) {
      bounds.push([Infinity, -Infinity]);
    }
    for (i = 0; i < lines.length; i++) {
      var line = lines[i];
      var nums = line.split(',').map(Number);
      var j;
      for (j = 0; j < 3; j++) {
        pos[3*i+j] = nums[j];
      }
      for (j = 0; j < n_col; j++) {
        if (nums[j] < bounds[j][0]) bounds[j][0] = nums[j];
        if (nums[j] > bounds[j][1]) bounds[j][1] = nums[j];
      }
      experiment_ids.push(nums[3]);
    }
    self.add_points(pos, experiment_ids);
    self.camera.zoom = 0.5 * (self.camera.top - self.camera.bottom);

    self.set_view(options);
    if (options.callback) options.callback();
    //self.redraw_center();
  });
};


var point_vert = [
  'uniform float size;',
  'varying vec3 vcolor;',
  'void main() {',
  '  vcolor = color;',
  '  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);',
  '  gl_PointSize = size;',
  '}'].join('\n');

var point_frag = [
  'uniform int show_only;',
  'varying vec3 vcolor;',
  'void main() {',
  '  if (show_only == -1 && vcolor.r != 1.0 || show_only == 0 && vcolor.g != 1.0) discard;',
  // not sure how portable it is
  '  vec2 diff = gl_PointCoord - vec2(0.5, 0.5);',
  '  float dist_sq = 4.0 * dot(diff, diff);',
  '  if (dist_sq >= 1.0) discard;',
  '  gl_FragColor = vec4(vcolor, 1.0 - dist_sq * dist_sq * dist_sq);',
  '}'].join('\n');


ReciprocalViewer.prototype.add_points = function (pos, experiment_ids) {
  var colors = new Float32Array(3 * experiment_ids.length);
  for (var i = 0; i < experiment_ids.length; i++) {
    var col = experiment_ids[i] === -1 ? [1, 0, 0] : [0, 1, 0];
    colors[3*i] = col[0];
    colors[3*i+1] = col[1];
    colors[3*i+2] = col[2];
  }

  var geometry = new THREE.BufferGeometry();
  geometry.addAttribute('position', new THREE.BufferAttribute(pos, 3));
  geometry.addAttribute('color', new THREE.BufferAttribute(colors, 3));
  var uniforms = {
    size: { value: 3 },
    show_only: { value: -2 },
  };
  var material = new THREE.ShaderMaterial({
    uniforms: uniforms,
    vertexShader: point_vert,
    fragmentShader: point_frag,
    vertexColors: THREE.VertexColors,
  });
  material.transparent = true;
  this.points = new THREE.Points(geometry, material);
  this.scene.add(this.points);
  this.request_render();
};

// temporary hack
ReciprocalViewer.prototype.change_isolevel_by = function (map_idx, delta) {
  this.change_zoom_by_factor(1 + delta);
};

ReciprocalViewer.prototype.redraw_center = function () {};

ReciprocalViewer.prototype.change_bond_line = function (delta) {
  if (this.points === null) return;
  var size = this.points.material.uniforms.size;
  size.value = Math.max(size.value + (delta > 0 ? 0.5 : -0.5), 0.5);
  this.hud('point size: ' + size.value.toFixed(1));
};
