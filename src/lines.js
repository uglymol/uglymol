var THREE = THREE || require('three'); // eslint-disable-line

var LineFactory = (function () {
'use strict';

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
  var side = new Float32Array(2*len);
  for (i = 0; i < len; i++) {
    side[2*i] = 1;
    side[2*i+1] = -1;
  }
  var previous = new Float32Array(6*len);
  for (i = 0; i < 6; i++) previous[i] = pos[i];
  for (; i < 6 * len; i++) previous[i] = pos[i-6];
  var next = new Float32Array(6*len);
  for (i = 0; i < 6 * (len-1); i++) next[i] = pos[i+6];
  for (; i < 6 * len; i++) next[i] = pos[i];
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
  '  vec4 diff = mat * vec4(next - previous, 0.0);',
  '  vec2 normal = normalize(vec2(-diff.y, diff.x));',
  '  vec4 pos = mat * vec4(position, 1.0);',
  '  pos.xy += side * linewidth * normal / size;',
  '  gl_Position = pos;',
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


function expand_vertices(segment, smooth) {
  var vertices = [];
  for (var i = 0; i < segment.length; i++) {
    var xyz = segment[i].xyz;
    vertices.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]));
  }
  if (!smooth || smooth < 2) return vertices;
  var curve = new THREE.CatmullRomCurve3(vertices);
  return curve.getPoints((segment.length - 1) * smooth);
}

function expand_colors(colors, smooth) {
  if (!smooth || smooth < 2) return colors;
  var ret = [];
  for (var i = 0; i < colors.length - 1; i++) {
    for (var j = 0; j < smooth; ++j) {
      ret.push(colors[i]);
    }
  }
  ret.push(colors[colors.length - 1]);
  return ret;
}

function LineFactory(param) {
  this.use_gl_lines = param.use_gl_lines;
  if (this.use_gl_lines) {
    this.material = new THREE.LineBasicMaterial({
      vertexColors: THREE.VertexColors,
      linewidth: param.linewidth
    });
  } else {
    this.material = new THREE.RawShaderMaterial({
      uniforms: {
        linewidth: { value: param.linewidth },
        color: { value: param.color },
        size: { value: param.size },
        fogNear: { value: null },  // will be updated in setProgram()
        fogFar: { value: null },
        fogColor: { value: null }
      },
      vertexShader: wide_line_vert,
      fragmentShader: wide_line_frag
    });
    this.material.fog = true;
  }
}

LineFactory.prototype.produce = function (vertices, colors, smoothness) {
  var vertex_arr = expand_vertices(vertices, smoothness);
  var color_arr = expand_colors(colors, smoothness);
  if (this.use_gl_lines) {
    var geom = new THREE.Geometry();
    geom.vertices = vertex_arr;
    geom.colors = color_arr;
    return new THREE.Line(geom, this.material);
  }
  var mesh = new THREE.Mesh(wide_line_geometry(vertex_arr, color_arr),
                            this.material);
  mesh.drawMode = THREE.TriangleStripDrawMode;
  return mesh;
};

return LineFactory;
})();

if (typeof module !== 'undefined') module.exports = LineFactory;
