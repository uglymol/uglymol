var THREE = THREE || require('three'); // eslint-disable-line

var LineFactory = (function () {
'use strict';

// two array must be of the same length
function wide_line_geometry(vertex_arr, color_arr) {
  var len = vertex_arr.length;
  var i;
  var pos = [];
  for (i = 0; i < len; i++) {
    var v = vertex_arr[i];
    pos.push(v.x, v.y, v.z);
    pos.push(v.x, v.y, v.z);
  }
  var position = new Float32Array(pos);
  // could we use three overlapping views of the same buffer?
  var previous = new Float32Array(6*len);
  for (i = 0; i < 6; i++) previous[i] = pos[i];
  for (; i < 6 * len; i++) previous[i] = pos[i-6];
  var next = new Float32Array(6*len);
  for (i = 0; i < 6 * (len-1); i++) next[i] = pos[i+6];
  for (; i < 6 * len; i++) next[i] = pos[i];
  var side = new Float32Array(2*len);
  for (i = 0; i < len; i++) {
    side[2*i] = 1;
    side[2*i+1] = -1;
  }
  var color = new Float32Array(6*len);
  for (i = 0; i < len; i++) {
    var col = color_arr[i];
    color[6*i] = col.r;
    color[6*i+1] = col.g;
    color[6*i+2] = col.b;
    color[6*i+3] = col.r;
    color[6*i+4] = col.g;
    color[6*i+5] = col.b;
  }
  var geometry = new THREE.BufferGeometry();
  geometry.addAttribute('position', new THREE.BufferAttribute(position, 3));
  geometry.addAttribute('previous', new THREE.BufferAttribute(previous, 3));
  geometry.addAttribute('next', new THREE.BufferAttribute(next, 3));
  geometry.addAttribute('side', new THREE.BufferAttribute(side, 1));
  geometry.addAttribute('color', new THREE.BufferAttribute(color, 3));
  return geometry;
}


var wide_line_vert = [
  'precision highp float;',

  'attribute vec3 position;',
  'attribute vec3 previous;',
  'attribute vec3 next;',
  'attribute float side;',
  'attribute vec3 color;',
  'uniform mat4 projectionMatrix;',
  'uniform mat4 modelViewMatrix;',
  'uniform vec2 size;',
  'uniform float linewidth;',
  'varying vec3 vcolor;',

  'void main() {',
  '  vcolor = color;',
  '  mat4 mat = projectionMatrix * modelViewMatrix;',
  '  vec2 dir1 = (mat * vec4(next - position, 0.0)).xy;',
  '  float len = length(dir1);',
  '  if (len > 0.0) dir1 /= len;',
  '  vec2 dir2 = (mat * vec4(position - previous, 0.0)).xy;',
  '  len = length(dir2);',
  '  dir2 = len > 0.0 ? dir2 / len : dir1;',
  '  vec2 tang = normalize(dir1 + dir2);',
  '  vec2 normal = vec2(-tang.y, tang.x);',
  // Now we have more or less a miter join. Bavel join could be more
  // appropriate, but it'd require one more triangle and more complex shader.
  // max() is a trade-off between too-long miters and too-thin lines.
  // The outer vertex should not go too far, the inner is not a problem.
  '  float angle_factor = max(dot(tang, dir2), side > 0.0 ? 0.5 : 0.1);',
  '  gl_Position = mat * vec4(position, 1.0);',
  '  gl_Position.xy += side * linewidth / angle_factor * normal / size;',
  '}'].join('\n');

var wide_line_frag = [
  'precision mediump float;',
  'uniform vec3 fogColor;',
  'uniform float fogNear;',
  'uniform float fogFar;',
  'varying vec3 vcolor;',
  'void main() {',
  '  gl_FragColor = vec4(vcolor, 1.0);',
  '  float depth = gl_FragCoord.z / gl_FragCoord.w;',
  '  float fogFactor = smoothstep(fogNear, fogFar, depth);',
  '  gl_FragColor.rgb = mix(gl_FragColor.rgb, fogColor, fogFactor);',
  '}'].join('\n');


function interpolate_vertices(segment, smooth) {
  var vertices = [];
  for (var i = 0; i < segment.length; i++) {
    var xyz = segment[i].xyz;
    vertices.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]));
  }
  if (!smooth || smooth < 2) return vertices;
  var curve = new THREE.CatmullRomCurve3(vertices);
  return curve.getPoints((segment.length - 1) * smooth);
}

function interpolate_colors(colors, smooth) {
  if (!smooth || smooth < 2) return colors;
  var ret = [];
  for (var i = 0; i < colors.length - 1; i++) {
    for (var j = 0; j < smooth; ++j) {
      // currently we don't interpolate them, just adjust length
      ret.push(colors[i]);
    }
  }
  ret.push(colors[colors.length - 1]);
  return ret;
}

function LineFactory(use_gl_lines, param) {
  this.use_gl_lines = use_gl_lines;
  if (use_gl_lines) {
    if (param.color === undefined) {
      param.vertexColors = THREE.VertexColors;
    }
    delete param.size; // only needed for RawShaderMaterial
    this.material = new THREE.LineBasicMaterial(param);
  } else {
    var uniforms = {
      fogNear: { value: null },  // will be updated in setProgram()
      fogFar: { value: null },
      fogColor: { value: null }
    };
    for (var p in param) {
      uniforms[p] = { value: param[p] };
    }
    this.material = new THREE.RawShaderMaterial({
      uniforms: uniforms,
      vertexShader: wide_line_vert,
      fragmentShader: wide_line_frag
    });
    this.material.fog = true;
  }
}

function rgb_to_buf(colors) {
  var arr = new Float32Array(colors.length * 3);
  for (var i = 0; i < colors.length; i++) {
    var c = colors[i];
    arr[3*i] = c.r;
    arr[3*i+1] = c.g;
    arr[3*i+2] = c.b;
  }
  return arr;
}

function xyz_to_buf(vectors) {
  var arr = new Float32Array(vectors.length * 3);
  for (var i = 0; i < vectors.length; i++) {
    var v = vectors[i];
    arr[3*i] = v.x;
    arr[3*i+1] = v.y;
    arr[3*i+2] = v.z;
  }
  return arr;
}

LineFactory.prototype.produce = function (vertices, colors, smoothness) {
  var vertex_arr = interpolate_vertices(vertices, smoothness);
  var color_arr = interpolate_colors(colors, smoothness);
  if (this.use_gl_lines) {
    var geometry = new THREE.BufferGeometry();
    var pos = xyz_to_buf(vertex_arr);
    geometry.addAttribute('position', new THREE.BufferAttribute(pos, 3));
    var col = rgb_to_buf(color_arr);
    geometry.addAttribute('color', new THREE.BufferAttribute(col, 3));
    return new THREE.Line(geometry, this.material);
  }
  var mesh = new THREE.Mesh(wide_line_geometry(vertex_arr, color_arr),
                            this.material);
  mesh.drawMode = THREE.TriangleStripDrawMode;
  return mesh;
};

LineFactory.make_chickenwire = function (data, parameters) {
  var geom = new THREE.BufferGeometry();
  var position = new Float32Array(data.vertices);
  geom.addAttribute('position', new THREE.BufferAttribute(position, 3));
  /* old version - mesh instead of lines
  geom.setIndex(new THREE.BufferAttribute(new Uint32Array(data.faces), 1));
  var material = new THREE.MeshBasicMaterial({
    color: this.config.colors[mtype],
    wireframe: true,
    wireframeLinewidth: this.config.map_line
  });
  var obj = new THREE.Mesh(geom, material);
  */

  // Although almost all browsers support OES_element_index_uint nowadays,
  // use Uint32 indexes only when needed.
  var arr = (data.vertices.length < 3*65536 ? new Uint16Array(data.segments)
                                            : new Uint32Array(data.segments));
  //console.log('arr len:', data.vertices.length, data.segments.length);
  geom.setIndex(new THREE.BufferAttribute(arr, 1));
  var material = new THREE.LineBasicMaterial(parameters);
  return new THREE.LineSegments(geom, material);
};

return LineFactory;
})();

if (typeof module !== 'undefined') module.exports = LineFactory;
