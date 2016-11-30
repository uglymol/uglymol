// @flow

import * as THREE from 'three';

/*:: type Num3 = [number, number, number] */
/*:: type Atom = {xyz: Num3} */
/*:: type Color = {r: number, g: number, b: number} */
/*:: type Vector3 = {x: number, y: number, z: number} */

const CUBE_EDGES = [[0, 0, 0], [1, 0, 0],
                    [0, 0, 0], [0, 1, 0],
                    [0, 0, 0], [0, 0, 1],
                    [1, 0, 0], [1, 1, 0],
                    [1, 0, 0], [1, 0, 1],
                    [0, 1, 0], [1, 1, 0],
                    [0, 1, 0], [0, 1, 1],
                    [0, 0, 1], [1, 0, 1],
                    [0, 0, 1], [0, 1, 1],
                    [1, 0, 1], [1, 1, 1],
                    [1, 1, 0], [1, 1, 1],
                    [0, 1, 1], [1, 1, 1]];

export function makeCube(size /*:number*/,
                         ctr /*:Vector3*/,
                         options /*:{[key:string]: any}*/) {
  const vertices = CUBE_EDGES.map(function (a) {
    return {
      x: ctr.x + size * (a[0] - 0.5),
      y: ctr.y + size * (a[1] - 0.5),
      z: ctr.z + size * (a[2] - 0.5)};
  });
  const material = makeLineMaterial({
    gl_lines: true,
    color: options.color,
    linewidth: options.linewidth,
    win_size: options.win_size,
    segments: true,
  });
  return makeLineSegments(material, vertices);
}

// A cube with 3 edges (for x, y, z axes) colored in red, green and blue.
export function makeRgbBox(transform_func /*:Num3 => Num3*/,
                           options /*:{[key:string]: any}*/) {
  const vertices = CUBE_EDGES.map(function (a) {
    return { xyz: transform_func(a) };
  });
  let colors = [
    new THREE.Color(0xff0000), new THREE.Color(0xffaa00),
    new THREE.Color(0x00ff00), new THREE.Color(0xaaff00),
    new THREE.Color(0x0000ff), new THREE.Color(0x00aaff),
  ];
  for (let j = 6; j < CUBE_EDGES.length; j++) {
    colors.push(options.color);
  }
  const material = makeLineMaterial({
    gl_lines: true,
    linewidth: 1,
    segments: true,
  });
  // flow-ignore-line - the type of vertices confuses flow
  return makeLineSegments(material, vertices, colors);
}

function double_pos(vertex_arr /*:Vector3[] | Atom[]*/) {
  let pos = [];
  let i;
  if (vertex_arr && vertex_arr[0].xyz) {
    for (i = 0; i < vertex_arr.length; i++) {
      // flow-ignore-line - disjoint unions not smart enough
      const xyz /*:Num3*/ = vertex_arr[i].xyz;
      pos.push(xyz[0], xyz[1], xyz[2]);
      pos.push(xyz[0], xyz[1], xyz[2]);
    }
  } else {
    for (i = 0; i < vertex_arr.length; i++) {
      // flow-ignore-line
      const v /*:Vector3*/ = vertex_arr[i];
      pos.push(v.x, v.y, v.z);
      pos.push(v.x, v.y, v.z);
    }
  }
  return pos;
}

function double_color(color_arr /*:Color[]*/) {
  const len = color_arr.length;
  let color = new Float32Array(6*len);
  for (let i = 0; i < len; i++) {
    const col = color_arr[i];
    color[6*i] = col.r;
    color[6*i+1] = col.g;
    color[6*i+2] = col.b;
    color[6*i+3] = col.r;
    color[6*i+4] = col.g;
    color[6*i+5] = col.b;
  }
  return color;
}

// input arrays must be of the same length
function wide_line_geometry(vertex_arr, color_arr) {
  const len = vertex_arr.length;
  const pos = double_pos(vertex_arr);
  const position = new Float32Array(pos);
  // could we use three overlapping views of the same buffer?
  let previous = new Float32Array(6*len);
  let i;
  for (i = 0; i < 6; i++) previous[i] = pos[i];
  for (; i < 6 * len; i++) previous[i] = pos[i-6];
  let next = new Float32Array(6*len);
  for (i = 0; i < 6 * (len-1); i++) next[i] = pos[i+6];
  for (; i < 6 * len; i++) next[i] = pos[i];
  let side = new Float32Array(2*len);
  for (i = 0; i < len; i++) {
    side[2*i] = 1;
    side[2*i+1] = -1;
  }
  const color = double_color(color_arr);
  let geometry = new THREE.BufferGeometry();
  geometry.addAttribute('position', new THREE.BufferAttribute(position, 3));
  geometry.addAttribute('previous', new THREE.BufferAttribute(previous, 3));
  geometry.addAttribute('next', new THREE.BufferAttribute(next, 3));
  geometry.addAttribute('side', new THREE.BufferAttribute(side, 1));
  geometry.addAttribute('color', new THREE.BufferAttribute(color, 3));
  return geometry;
}

// input arrays must be of the same length
function wide_segments_geometry(vertex_arr, color_arr) {
  // n input vertices => 2n output vertices, n triangles, 3n indexes
  const len = vertex_arr.length;
  let i;
  let j;
  const pos = double_pos(vertex_arr);
  const position = new Float32Array(pos);
  let other_vert = new Float32Array(6*len);
  for (i = 0; i < 6 * len; i += 12) {
    for (j = 0; j < 6; j++) other_vert[i+j] = pos[i+j+6];
    for (; j < 12; j++) other_vert[i+j] = pos[i+j-6];
  }
  let side = new Float32Array(2*len);
  for (i = 0; i < len; i++) {
    side[2*i] = -1;
    side[2*i+1] = 1;
  }
  let index = (2*len < 65536 ? new Uint16Array(3*len)
                             : new Uint32Array(3*len));
  const vert_order = [0, 1, 2, 0, 2, 3];
  for (i = 0; i < len / 2; i++) {
    for (j = 0; j < 6; j++) {
      index[6*i+j] = 4*i + vert_order[j];
    }
  }
  let geometry = new THREE.BufferGeometry();
  geometry.addAttribute('position', new THREE.BufferAttribute(position, 3));
  geometry.addAttribute('other', new THREE.BufferAttribute(other_vert, 3));
  geometry.addAttribute('side', new THREE.BufferAttribute(side, 1));
  if (color_arr != null) {
    const color = double_color(color_arr);
    geometry.addAttribute('color', new THREE.BufferAttribute(color, 3));
  }
  geometry.setIndex(new THREE.BufferAttribute(index, 1));
  return geometry;
}


const wide_line_vert = [
  'attribute vec3 previous;',
  'attribute vec3 next;',
  'attribute float side;',
  'uniform vec2 win_size;',
  'uniform float linewidth;',
  'varying vec3 vcolor;',

  'void main() {',
  '  vcolor = color;',
  '  mat4 mat = projectionMatrix * modelViewMatrix;',
  '  vec2 dir1 = (mat * vec4(next - position, 0.0)).xy * win_size;',
  '  float len = length(dir1);',
  '  if (len > 0.0) dir1 /= len;',
  '  vec2 dir2 = (mat * vec4(position - previous, 0.0)).xy * win_size;',
  '  len = length(dir2);',
  '  dir2 = len > 0.0 ? dir2 / len : dir1;',
  '  vec2 tang = normalize(dir1 + dir2);',
  '  vec2 normal = vec2(-tang.y, tang.x);',
  // Now we have more or less a miter join. Bavel join could be more
  // appropriate, but it'd require one more triangle and more complex shader.
  // max() is a trade-off between too-long miters and too-thin lines.
  // The outer vertex should not go too far, the inner is not a problem.
  '  float outer = side * dot(dir2, normal);',
  '  float angle_factor = max(dot(tang, dir2), outer > 0.0 ? 0.5 : 0.1);',
  '  gl_Position = mat * vec4(position, 1.0);',
  '  gl_Position.xy += side * linewidth / angle_factor * normal / win_size;',
  '}'].join('\n');

const wide_segments_vert = [
  'attribute vec3 other;',
  'attribute float side;',
  'uniform vec2 win_size;',
  'uniform float linewidth;',
  'varying vec3 vcolor;',

  'void main() {',
  '  vcolor = color;',
  '  mat4 mat = projectionMatrix * modelViewMatrix;',
  '  vec2 dir = normalize((mat * vec4(position - other, 0.0)).xy);',
  '  vec2 normal = vec2(-dir.y, dir.x);',
  '  gl_Position = mat * vec4(position, 1.0);',
  '  gl_Position.xy += side * linewidth * normal / win_size;',
  '}'].join('\n');

const wide_line_frag = [
  '#include <fog_pars_fragment>',
  'varying vec3 vcolor;',
  'void main() {',
  '  gl_FragColor = vec4(vcolor, 1.0);',
  '#include <fog_fragment>',
  '}'].join('\n');


function interpolate_vertices(segment, smooth) /*:Vector3[]*/{
  let vertices = [];
  for (let i = 0; i < segment.length; i++) {
    const xyz = segment[i].xyz;
    vertices.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]));
  }
  if (!smooth || smooth < 2) return vertices;
  const curve = new THREE.CatmullRomCurve3(vertices);
  return curve.getPoints((segment.length - 1) * smooth);
}

function interpolate_colors(colors, smooth) {
  if (!smooth || smooth < 2) return colors;
  let ret = [];
  for (let i = 0; i < colors.length - 1; i++) {
    for (let j = 0; j < smooth; j++) {
      // currently we don't really interpolate colors
      ret.push(colors[i]);
    }
  }
  ret.push(colors[colors.length - 1]);
  return ret;
}

// a simplistic linear interpolation, no need to SLERP
function interpolate_directions(dirs, smooth) {
  smooth = smooth || 1;
  let ret = [];
  let i;
  for (i = 0; i < dirs.length - 1; i++) {
    const p = dirs[i];
    const n = dirs[i+1];
    for (let j = 0; j < smooth; j++) {
      const an = j / smooth;
      const ap = 1 - an;
      ret.push(ap*p[0] + an*n[0], ap*p[1] + an*n[1], ap*p[2] + an*n[2]);
    }
  }
  ret.push(dirs[i][0], dirs[i][1], dirs[i][2]);
  return ret;
}

function make_uniforms(params) {
  let uniforms = {
    fogNear: { value: null },  // will be updated in setProgram()
    fogFar: { value: null },
    fogColor: { value: null },
  };
  for (let p in params) {  // eslint-disable-line guard-for-in
    uniforms[p] = { value: params[p] };
  }
  return uniforms;
}

const ribbon_vert = [
  //'attribute vec3 normal;' is added by default for ShaderMaterial
  'uniform float shift;',
  'varying vec3 vcolor;',
  'void main() {',
  '  vcolor = color;',
  '  vec3 pos = position + shift * normalize(normal);',
  '  gl_Position = projectionMatrix * modelViewMatrix * vec4(pos, 1.0);',
  '}'].join('\n');

const ribbon_frag = [
  '#include <fog_pars_fragment>',
  'varying vec3 vcolor;',
  'void main() {',
  '  gl_FragColor = vec4(vcolor, 1.0);',
  '#include <fog_fragment>',
  '}'].join('\n');

// 9-line ribbon
export function makeRibbon(vertices /*:Atom[]*/,
                           colors /*:Color[]*/,
                           tangents /*:Num3[]*/,
                           smoothness /*:number*/) {
  const vertex_arr = interpolate_vertices(vertices, smoothness);
  const color_arr = interpolate_colors(colors, smoothness);
  const tang_arr = interpolate_directions(tangents, smoothness);
  let obj = new THREE.Object3D();
  let geometry = makeSimpleGeometry(vertex_arr, color_arr);
  const tan = new Float32Array(tang_arr);
  // it's not 'normal', but it doesn't matter
  geometry.addAttribute('normal', new THREE.BufferAttribute(tan, 3));
  const material0 = new THREE.ShaderMaterial({
    uniforms: make_uniforms({shift: 0}),
    vertexShader: ribbon_vert,
    fragmentShader: ribbon_frag,
    fog: true,
    vertexColors: THREE.VertexColors,
  });
  for (let n = -4; n < 5; n++) {
    let material = n === 0 ? material0 : material0.clone();
    material.uniforms.shift.value = 0.1 * n;
    obj.add(new THREE.Line(geometry, material));
  }
  return obj;
}


export
function makeChickenWire(data /*:{vertices: number[], segments: number[]}*/,
                         parameters /*:{[key: string]: any}*/) {
  let geom = new THREE.BufferGeometry();
  const position = new Float32Array(data.vertices);
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
  const arr = (data.vertices.length < 3*65536 ? new Uint16Array(data.segments)
                                              : new Uint32Array(data.segments));
  //console.log('arr len:', data.vertices.length, data.segments.length);
  geom.setIndex(new THREE.BufferAttribute(arr, 1));
  let material = new THREE.LineBasicMaterial(parameters);
  return new THREE.LineSegments(geom, material);
}


const grid_vert = [
  'uniform vec3 ucolor;',
  'uniform vec3 fogColor;',
  'varying vec4 vcolor;',
  'void main() {',
  '  vec2 scale = vec2(projectionMatrix[0][0], projectionMatrix[1][1]);',
  '  float z = position.z;',
  '  float fogFactor = (z > 0.5 ? 0.2 : 0.7);',
  '  float alpha = 0.8 * smoothstep(z > 1.5 ? -10.0 : 0.01, 0.1, scale.y);',
  '  vcolor = vec4(mix(ucolor, fogColor, fogFactor), alpha);',
  '  gl_Position = vec4(position.xy * scale, -0.99, 1.0);',
  '}'].join('\n');

const grid_frag = [
  'varying vec4 vcolor;',
  'void main() {',
  '  gl_FragColor = vcolor;',
  '}'].join('\n');

export function makeGrid() {
  const N = 50;
  let pos = [];
  for (let i = -N; i <= N; i++) {
    let z = 0; // z only marks major/minor axes
    if (i % 5 === 0) z = i % 2 === 0 ? 2 : 1;
    pos.push(-N, i, z, N, i, z);  // horizontal line
    pos.push(i, -N, z, i, N, z);  // vertical line
  }
  let geom = new THREE.BufferGeometry();
  geom.addAttribute('position',
                    new THREE.BufferAttribute(new Float32Array(pos), 3));
  let material = new THREE.ShaderMaterial({
    uniforms: make_uniforms({ucolor: new THREE.Color(0x888888)}),
    //linewidth: 3,
    vertexShader: grid_vert,
    fragmentShader: grid_frag,
    fog: true, // no really, but we use fogColor
  });
  material.transparent = true;
  let obj = new THREE.LineSegments(geom, material);
  obj.frustumCulled = false;  // otherwise the renderer could skip it
  obj.color_value = material.uniforms.ucolor.value; // shortcut
  return obj;
}


function makeSimpleLineMaterial(options) {
  let mparams = {};
  mparams.linewidth = options.linewidth;
  if (options.color === undefined) {
    mparams.vertexColors = THREE.VertexColors;
  } else {
    mparams.color = options.color;
  }
  return new THREE.LineBasicMaterial(mparams);
}

function makeThickLineMaterial(options) {
  let uniforms = make_uniforms({
    linewidth: options.linewidth,
    win_size: options.win_size,
  });
  return new THREE.ShaderMaterial({
    uniforms: uniforms,
    vertexShader: options.segments ? wide_segments_vert : wide_line_vert,
    fragmentShader: wide_line_frag,
    fog: true,
    vertexColors: THREE.VertexColors,
  });
}

export function makeLineMaterial(options /*:{[key: string]: mixed}*/) {
  return options.gl_lines ? makeSimpleLineMaterial(options)
                          : makeThickLineMaterial(options);
}

function makeSimpleGeometry(vertices /*:Vector3[] | Atom[]*/,
                            colors /*:?Color[]*/) {
  let geometry = new THREE.BufferGeometry();
  const pos = new Float32Array(vertices.length * 3);
  let i;
  if (vertices && vertices[0].xyz) {
    for (i = 0; i < vertices.length; i++) {
      // flow-ignore-line - disjoint unions not smart enough
      const xyz /*:Num3*/ = vertices[i].xyz;
      pos[3*i] = xyz[0];
      pos[3*i+1] = xyz[1];
      pos[3*i+2] = xyz[2];
    }
  } else {
    for (i = 0; i < vertices.length; i++) {
      // flow-ignore-line
      const v /*:Vector3*/ = vertices[i];
      pos[3*i] = v.x;
      pos[3*i+1] = v.y;
      pos[3*i+2] = v.z;
    }
  }
  geometry.addAttribute('position', new THREE.BufferAttribute(pos, 3));
  if (colors != null) {
    const col = new Float32Array(colors.length * 3);
    for (i = 0; i < colors.length; i++) {
      const c = colors[i];
      col[3*i] = c.r;
      col[3*i+1] = c.g;
      col[3*i+2] = c.b;
    }
    geometry.addAttribute('color', new THREE.BufferAttribute(col, 3));
  }
  return geometry;
}

function makeThickLine(material, vertices, colors) {
  let mesh = new THREE.Mesh(wide_line_geometry(vertices, colors), material);
  mesh.drawMode = THREE.TriangleStripDrawMode;
  mesh.raycast = line_raycast;
  return mesh;
}

export function makeLine(material /*:THREE.Material*/,
                         vertices /*:Vector3[]*/,
                         colors /*:Color[]*/) {
  if (material.isShaderMaterial) {
    return makeThickLine(material, vertices, colors);
  } else {
    return new THREE.Line(makeSimpleGeometry(vertices, colors), material);
  }
}

function makeThickLineSegments(material, vertices, colors) {
  let mesh = new THREE.Mesh(wide_segments_geometry(vertices, colors), material);
  mesh.raycast = line_raycast;
  return mesh;
}

export function makeLineSegments(material /*:THREE.Material*/,
                                 vertices /*:Vector3[] | Atom[]*/,
                                 colors /*:?Color[]*/) {
  if (material.isShaderMaterial) {
    return makeThickLineSegments(material, vertices, colors);
  } else {
    return new THREE.LineSegments(makeSimpleGeometry(vertices, colors),
                                  material);
  }
}

const wheel_vert = [
  'uniform float size;',
  'varying vec3 vcolor;',
  'void main() {',
  '  vcolor = color;',
  '  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);',
  '  gl_PointSize = size;',
  '}'].join('\n');

// not sure how portable it is
const wheel_frag = [
  '#include <fog_pars_fragment>',
  'varying vec3 vcolor;',
  'void main() {',
  '  vec2 diff = gl_PointCoord - vec2(0.5, 0.5);',
  '  if (dot(diff, diff) >= 0.25) discard;',
  '  gl_FragColor = vec4(vcolor, 1.0);',
  '#include <fog_fragment>',
  '}'].join('\n');

export function makeWheels(atom_arr /*:Atom[]*/,
                           color_arr /*:Color[]*/,
                           size /*:number*/) {
  let geometry = makeSimpleGeometry(atom_arr, color_arr);
  let material = new THREE.ShaderMaterial({
    uniforms: make_uniforms({size: size}),
    vertexShader: wheel_vert,
    fragmentShader: wheel_frag,
    fog: true,
    vertexColors: THREE.VertexColors,
  });
  let obj = new THREE.Points(geometry, material);
  // currently we use only lines for picking
  obj.raycast = function () {};
  return obj;
}


// based on THREE.Line.prototype.raycast(), but skipping duplicated points
let inverseMatrix = new THREE.Matrix4();
let ray = new THREE.Ray();
// this function will be put on prototype
/* eslint-disable no-invalid-this */
function line_raycast(raycaster, intersects) {
  const precisionSq = raycaster.linePrecision * raycaster.linePrecision;
  inverseMatrix.getInverse(this.matrixWorld);
  ray.copy(raycaster.ray).applyMatrix4(inverseMatrix);
  let vStart = new THREE.Vector3();
  let vEnd = new THREE.Vector3();
  let interSegment = new THREE.Vector3();
  let interRay = new THREE.Vector3();
  const step = this.drawMode === THREE.TriangleStripDrawMode ? 1 : 2;
  const positions = this.geometry.attributes.position.array;
  for (let i = 0, l = positions.length / 6 - 1; i < l; i += step) {
    vStart.fromArray(positions, 6 * i);
    vEnd.fromArray(positions, 6 * i + 6);
    let distSq = ray.distanceSqToSegment(vStart, vEnd, interRay, interSegment);
    if (distSq > precisionSq) continue;
    interRay.applyMatrix4(this.matrixWorld);
    const distance = raycaster.ray.origin.distanceTo(interRay);
    if (distance < raycaster.near || distance > raycaster.far) continue;
    intersects.push({
      distance: distance,
      point: interSegment.clone().applyMatrix4(this.matrixWorld),
      index: i,
      object: this,
      line_dist: Math.sqrt(distSq), // extra property, not in Three.js
    });
  }
}

function makeCanvasWithText(text, options) {
  if (typeof document === 'undefined') return;  // for testing on node
  let canvas = document.createElement('canvas');
  // Canvas size should be 2^N.
  canvas.width = 256;  // arbitrary limit, to keep it simple
  canvas.height = 16;  // font size
  let context = canvas.getContext('2d');
  if (!context) return null;
  context.font = (options.font || 'bold 14px') + ' sans-serif';
  //context.fillStyle = 'green';
  //context.fillRect(0, 0, canvas.width, canvas.height);
  context.textBaseline = 'bottom';
  if (options.color) context.fillStyle = options.color;
  context.fillText(text, 0, canvas.height);
  return canvas;
}

const label_vert = [
  'uniform vec2 canvas_size;',
  'uniform vec2 win_size;',
  'varying vec2 vUv;',
  'void main() {',
  '  vUv = uv;',
  '  vec2 rel_offset = vec2(0.02, -0.3);',
  '  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);',
  '  gl_Position.xy += (uv + rel_offset) * 2.0 * canvas_size / win_size;',
  '  gl_Position.z += 1.0 * projectionMatrix[2][2];',
  '}'].join('\n');

const label_frag = [
  '#include <fog_pars_fragment>',
  'varying vec2 vUv;',
  'uniform sampler2D map;',
  'void main() {',
  '  gl_FragColor = texture2D(map, vUv);',
  '#include <fog_fragment>',
  '}'].join('\n');


export function makeLabel(text /*:string*/, options /*:{[key:string]: any}*/) {
  const canvas = makeCanvasWithText(text, options);
  if (!canvas) return;
  let texture = new THREE.Texture(canvas);
  texture.needsUpdate = true;

  // Rectangle geometry.
  let geometry = new THREE.BufferGeometry();
  const pos = options.pos;
  const position = new Float32Array([].concat(pos, pos, pos, pos));
  const uvs = new Float32Array([0, 1, 1, 1, 0, 0, 1, 0]);
  const indices = new Uint16Array([0, 2, 1, 2, 3, 1]);
  geometry.setIndex(new THREE.BufferAttribute(indices, 1));
  geometry.addAttribute('position', new THREE.BufferAttribute(position, 3));
  geometry.addAttribute('uv', new THREE.BufferAttribute(uvs, 2));

  let material = new THREE.ShaderMaterial({
    uniforms: make_uniforms({map: texture,
                             canvas_size: [canvas.width, canvas.height],
                             win_size: options.win_size}),
    vertexShader: label_vert,
    fragmentShader: label_frag,
    fog: true,
  });
  material.transparent = true;
  let mesh = new THREE.Mesh(geometry, material);
  mesh.remake = function (text, options) {
    texture.image = makeCanvasWithText(text, options);
    texture.needsUpdate = true;
  };
  return mesh;
}

// Add vertices of a 3d cross (representation of an unbonded atom)
export
function addXyzCross(vertices /*:Vector3[]*/, xyz /*:Num3*/, r /*:number*/) {
  vertices.push(new THREE.Vector3(xyz[0]-r, xyz[1], xyz[2]));
  vertices.push(new THREE.Vector3(xyz[0]+r, xyz[1], xyz[2]));
  vertices.push(new THREE.Vector3(xyz[0], xyz[1]-r, xyz[2]));
  vertices.push(new THREE.Vector3(xyz[0], xyz[1]+r, xyz[2]));
  vertices.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]-r));
  vertices.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]+r));
}
