// @flow
//import { Viewer } from './viewer.js';
//import * as THREE from 'three';
var Viewer = UM.Viewer;

//export
function ReciprocalViewer(options /*: {[key: string]: any}*/) {
  Viewer.call(this, options);
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

    var frag = {};//parse_url_fragment();
    if (frag.zoom) self.camera.zoom = frag.zoom;
    self.recenter(frag.xyz || [0, 0, 0], frag.eye, 1);
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
  //'  gl_PointSize = size;',
  '}'].join('\n');

var point_frag = [
  'varying vec3 vcolor;',
  'void main() {',
  // not sure how portable it is
  //'  vec2 diff = gl_PointCoord - vec2(0.5, 0.5);',
  //'  if (dot(diff, diff) >= 0.25) discard;',
  '  gl_FragColor = vec4(vcolor, 1.0);',
  '}'].join('\n');


ReciprocalViewer.prototype.add_points = function (pos, experiment_ids) {
  var colors = new Float32Array(3 * experiment_ids.length);
  for (var i = 0; i < experiment_ids.length; i++) {
    var col = experiment_ids[i] === -1 ? [0.9, 0.1, 0.1] : [0.2, 0.8, 0.2];
    colors[3*i] = col[0];
    colors[3*i+1] = col[1];
    colors[3*i+2] = col[2];
  }

  var geometry = new THREE.BufferGeometry();
  geometry.addAttribute('position', new THREE.BufferAttribute(pos, 3));
  geometry.addAttribute('color', new THREE.BufferAttribute(colors, 3));
  var uniforms = {
    size: { value: 2 },
  };
  var material = new THREE.ShaderMaterial({
    uniforms: uniforms,
    vertexShader: point_vert,
    fragmentShader: point_frag,
    vertexColors: THREE.VertexColors,
  });
  var obj = new THREE.Points(geometry, material);
  this.scene.add(obj);
  this.request_render();
};

// temporary hack
ReciprocalViewer.prototype.change_isolevel_by = function (map_idx, delta) {
  this.change_zoom_by_factor(1 + delta);
};

ReciprocalViewer.prototype.redraw_center = function () {};
