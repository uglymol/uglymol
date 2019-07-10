// @flow

import { BufferAttribute, BufferGeometry, ShaderMaterial,
         Object3D, Mesh, Line, LineSegments, Points,
         Color, Vector3, Matrix4, Ray, CatmullRomCurve3,
         TriangleStripDrawMode, Texture } from './fromthree.js';

/*:: type Num3 = [number, number, number] */
/*:: import type {AtomT} from './model.js' */

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

function makeColorAttribute(colors /*:Color[]*/) {
  const col = new Float32Array(colors.length * 3);
  for (let i = 0; i < colors.length; i++) {
    col[3*i+0] = colors[i].r;
    col[3*i+1] = colors[i].g;
    col[3*i+2] = colors[i].b;
  }
  return new BufferAttribute(col, 3);
}

const light_dir = new Vector3(-0.2, 0.3, 1.0); // length affects brightness

export const fog_pars_fragment =
`#ifdef USE_FOG
uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;
#endif`;

export const fog_end_fragment =
`#ifdef USE_FOG
  float depth = gl_FragCoord.z / gl_FragCoord.w;
  float fogFactor = smoothstep(fogNear, fogFar, depth);
  gl_FragColor.rgb = mix(gl_FragColor.rgb, fogColor, fogFactor);
#endif`;


const varcolor_vert = `
attribute vec3 color;
varying vec3 vcolor;
void main() {
  vcolor = color;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
`;

const unicolor_vert = `
void main() {
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
`;

const unicolor_frag = `
${fog_pars_fragment}
uniform vec3 vcolor;
void main() {
  gl_FragColor = vec4(vcolor, 1.0);
${fog_end_fragment}
}`;

const varcolor_frag = `
${fog_pars_fragment}
varying vec3 vcolor;
void main() {
  gl_FragColor = vec4(vcolor, 1.0);
${fog_end_fragment}
}`;

export function makeLines(pos /*:Float32Array*/, color /*:Color*/,
                          linewidth /*:number*/) {
  const material = new ShaderMaterial({
    uniforms: makeUniforms({vcolor: color}),
    vertexShader: unicolor_vert,
    fragmentShader: unicolor_frag,
    fog: true,
    linewidth: linewidth,
    type: 'um_lines',
  });
  let geometry = new BufferGeometry();
  geometry.addAttribute('position', new BufferAttribute(pos, 3));
  return new LineSegments(geometry, material);
}

export function makeCube(size /*:number*/,
                         ctr /*:Vector3*/,
                         options /*:{[key:string]: any}*/) {
  const pos = new Float32Array(CUBE_EDGES.length * 3);
  for (let i = 0; i < CUBE_EDGES.length; i++) {
    let coor = CUBE_EDGES[i];
    pos[3*i+0] = ctr.x + size * (coor[0] - 0.5);
    pos[3*i+1] = ctr.y + size * (coor[1] - 0.5);
    pos[3*i+2] = ctr.z + size * (coor[2] - 0.5);
  }
  return makeLines(pos, options.color, options.linewidth);
}

export function makeMultiColorLines(pos /*:Float32Array*/,
                                    colors /*:Color[]*/,
                                    linewidth /*:number*/) {
  const material = new ShaderMaterial({
    uniforms: makeUniforms({}),
    vertexShader: varcolor_vert,
    fragmentShader: varcolor_frag,
    fog: true,
    linewidth: linewidth,
    type: 'um_multicolor_lines',
  });
  let geometry = new BufferGeometry();
  geometry.addAttribute('position', new BufferAttribute(pos, 3));
  geometry.addAttribute('color', makeColorAttribute(colors));
  return new LineSegments(geometry, material);
}

// A cube with 3 edges (for x, y, z axes) colored in red, green and blue.
export function makeRgbBox(transform_func /*:Num3 => Num3*/, color /*:Color*/) {
  const pos = new Float32Array(CUBE_EDGES.length * 3);
  for (let i = 0; i < CUBE_EDGES.length; i++) {
    let coor = transform_func(CUBE_EDGES[i]);
    pos[3*i+0] = coor[0];
    pos[3*i+1] = coor[1];
    pos[3*i+2] = coor[2];
  }
  let colors = [
    new Color(0xff0000), new Color(0xffaa00),
    new Color(0x00ff00), new Color(0xaaff00),
    new Color(0x0000ff), new Color(0x00aaff),
  ];
  for (let j = 6; j < CUBE_EDGES.length; j++) {
    colors.push(color);
  }
  return makeMultiColorLines(pos, colors, 1);
}

function double_pos(pos /*:Num3[]*/) {
  let double_pos = [];
  for (let i = 0; i < pos.length; i++) {
    const v = pos[i];
    double_pos.push(v[0], v[1], v[2]);
    double_pos.push(v[0], v[1], v[2]);
  }
  return double_pos;
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

// draw quads as 2 triangles: 4 attributes / quad, 6 indices / quad
function make_quad_index_buffer(len) {
  let index = (4*len < 65536 ? new Uint16Array(6*len)
                             : new Uint32Array(6*len));
  const vert_order = [0, 1, 2, 0, 2, 3];
  for (let i = 0; i < len; i++) {
    for (let j = 0; j < 6; j++) {
      index[6*i+j] = 4*i + vert_order[j];
    }
  }
  return new BufferAttribute(index, 1);
}


const wide_line_vert = [
  'attribute vec3 color;',
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

const wide_segments_vert = `
attribute vec3 color;
attribute vec3 other;
attribute float side;
uniform vec2 win_size;
uniform float linewidth;
varying vec3 vcolor;

void main() {
  vcolor = color;
  mat4 mat = projectionMatrix * modelViewMatrix;
  vec2 dir = normalize((mat * vec4(position - other, 0.0)).xy);
  vec2 normal = vec2(-dir.y, dir.x);
  gl_Position = mat * vec4(position, 1.0);
  gl_Position.xy += side * linewidth * normal / win_size;
}`;

function interpolate_vertices(segment, smooth) /*:Vector3[]*/{
  let vertices = [];
  for (let i = 0; i < segment.length; i++) {
    const xyz = segment[i].xyz;
    vertices.push(new Vector3(xyz[0], xyz[1], xyz[2]));
  }
  if (!smooth || smooth < 2) return vertices;
  const curve = new CatmullRomCurve3(vertices);
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

export function makeUniforms(params/*:{[id:string]:mixed}*/) {
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

const ribbon_vert = `
attribute vec3 color;
attribute vec3 tan;
uniform float shift;
varying vec3 vcolor;
void main() {
  vcolor = color;
  vec3 pos = position + shift * normalize(tan);
  gl_Position = projectionMatrix * modelViewMatrix * vec4(pos, 1.0);
}`;

// 9-line ribbon
export function makeRibbon(vertices /*:AtomT[]*/,
                           colors /*:Color[]*/,
                           tangents /*:Num3[]*/,
                           smoothness /*:number*/) {
  const vertex_arr = interpolate_vertices(vertices, smoothness);
  const color_arr = interpolate_colors(colors, smoothness);
  const tang_arr = interpolate_directions(tangents, smoothness);
  let obj = new Object3D();
  let geometry = new BufferGeometry();
  const pos = new Float32Array(vertex_arr.length * 3);
  for (let i = 0; i < vertex_arr.length; i++) {
    const v = vertex_arr[i];
    pos[3*i+0] = v.x;
    pos[3*i+1] = v.y;
    pos[3*i+2] = v.z;
  }
  geometry.addAttribute('position', new BufferAttribute(pos, 3));
  geometry.addAttribute('color', makeColorAttribute(color_arr));
  const tan = new Float32Array(tang_arr);
  geometry.addAttribute('tan', new BufferAttribute(tan, 3));
  for (let n = -4; n < 5; n++) {
    const material = new ShaderMaterial({
      uniforms: makeUniforms({shift: 0.1 * n}),
      vertexShader: ribbon_vert,
      fragmentShader: varcolor_frag,
      fog: true,
      type: 'um_ribbon',
    });
    obj.add(new Line(geometry, material));
  }
  return obj;
}


export
function makeChickenWire(data /*:{vertices: number[], segments: number[]}*/,
                         options /*:{[key: string]: any}*/) {
  let geom = new BufferGeometry();
  const position = new Float32Array(data.vertices);
  geom.addAttribute('position', new BufferAttribute(position, 3));

  // Although almost all browsers support OES_element_index_uint nowadays,
  // use Uint32 indexes only when needed.
  const arr = (data.vertices.length < 3*65536 ? new Uint16Array(data.segments)
                                              : new Uint32Array(data.segments));
  //console.log('arr len:', data.vertices.length, data.segments.length);
  geom.setIndex(new BufferAttribute(arr, 1));
  const material = new ShaderMaterial({
    uniforms: makeUniforms({vcolor: options.color}),
    vertexShader: unicolor_vert,
    fragmentShader: unicolor_frag,
    fog: true,
    linewidth: options.linewidth,
    type: 'um_line_chickenwire',
  });
  return new LineSegments(geom, material);
}


const grid_vert = `
uniform vec3 ucolor;
uniform vec3 fogColor;
varying vec4 vcolor;
void main() {
  vec2 scale = vec2(projectionMatrix[0][0], projectionMatrix[1][1]);
  float z = position.z;
  float fogFactor = (z > 0.5 ? 0.2 : 0.7);
  float alpha = 0.8 * smoothstep(z > 1.5 ? -10.0 : 0.01, 0.1, scale.y);
  vcolor = vec4(mix(ucolor, fogColor, fogFactor), alpha);
  gl_Position = vec4(position.xy * scale, -0.99, 1.0);
}`;

const grid_frag = `
varying vec4 vcolor;
void main() {
  gl_FragColor = vcolor;
}`;

export function makeGrid() {
  const N = 50;
  let pos = [];
  for (let i = -N; i <= N; i++) {
    let z = 0; // z only marks major/minor axes
    if (i % 5 === 0) z = i % 2 === 0 ? 2 : 1;
    pos.push(-N, i, z, N, i, z);  // horizontal line
    pos.push(i, -N, z, i, N, z);  // vertical line
  }
  let geom = new BufferGeometry();
  const pos_arr = new Float32Array(pos);
  geom.addAttribute('position', new BufferAttribute(pos_arr, 3));
  let material = new ShaderMaterial({
    uniforms: makeUniforms({ucolor: new Color(0x888888)}),
    //linewidth: 3,
    vertexShader: grid_vert,
    fragmentShader: grid_frag,
    fog: true, // no really, but we use fogColor
    type: 'um_grid',
  });
  material.transparent = true;
  let obj = new LineSegments(geom, material);
  obj.frustumCulled = false;  // otherwise the renderer could skip it
  obj.color_value = material.uniforms.ucolor.value; // shortcut
  return obj;
}


export function makeLineMaterial(options /*:{[key: string]: mixed}*/) {
  let uniforms = makeUniforms({
    linewidth: options.linewidth,
    win_size: options.win_size,
  });
  return new ShaderMaterial({
    uniforms: uniforms,
    vertexShader: options.segments ? wide_segments_vert : wide_line_vert,
    fragmentShader: varcolor_frag,
    fog: true,
    type: 'um_line',
  });
}

// vertex_arr and color_arr must be of the same length
export function makeLine(material /*:ShaderMaterial*/,
                         vertex_arr /*:Num3[]*/,
                         color_arr /*:Color[]*/) {
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
  let geometry = new BufferGeometry();
  geometry.addAttribute('position', new BufferAttribute(position, 3));
  geometry.addAttribute('previous', new BufferAttribute(previous, 3));
  geometry.addAttribute('next', new BufferAttribute(next, 3));
  geometry.addAttribute('side', new BufferAttribute(side, 1));
  geometry.addAttribute('color', new BufferAttribute(color, 3));

  let mesh = new Mesh(geometry, material);
  mesh.drawMode = TriangleStripDrawMode;
  mesh.raycast = line_raycast;
  return mesh;
}

// vertex_arr and color_arr must be of the same length
export function makeLineSegments(material /*:ShaderMaterial*/,
                                 vertex_arr /*:Num3[]*/,
                                 color_arr /*:?Color[]*/) {
  // n input vertices => 2n output vertices, n triangles, 3n indexes
  const len = vertex_arr.length;
  const pos = double_pos(vertex_arr);
  const position = new Float32Array(pos);
  let other_vert = new Float32Array(6*len);
  for (let i = 0; i < 6 * len; i += 12) {
    let j = 0;
    for (; j < 6; j++) other_vert[i+j] = pos[i+j+6];
    for (; j < 12; j++) other_vert[i+j] = pos[i+j-6];
  }
  let side = new Float32Array(2*len);
  for (let k = 0; k < len; k++) {
    side[2*k] = -1;
    side[2*k+1] = 1;
  }
  let geometry = new BufferGeometry();
  geometry.addAttribute('position', new BufferAttribute(position, 3));
  geometry.addAttribute('other', new BufferAttribute(other_vert, 3));
  geometry.addAttribute('side', new BufferAttribute(side, 1));
  if (color_arr != null) {
    const color = double_color(color_arr);
    geometry.addAttribute('color', new BufferAttribute(color, 3));
  }
  geometry.setIndex(make_quad_index_buffer(len/2));

  let mesh = new Mesh(geometry, material);
  mesh.raycast = line_raycast;
  return mesh;
}

const wheel_vert = `
attribute vec3 color;
uniform float size;
varying vec3 vcolor;
void main() {
  vcolor = color;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  gl_PointSize = size;
}`;

// not sure how portable it is
const wheel_frag = `
${fog_pars_fragment}
varying vec3 vcolor;
void main() {
  vec2 diff = gl_PointCoord - vec2(0.5, 0.5);
  if (dot(diff, diff) >= 0.25) discard;
  gl_FragColor = vec4(vcolor, 1.0);
${fog_end_fragment}
}`;

export function makeWheels(atom_arr /*:AtomT[]*/,
                           color_arr /*:Color[]*/,
                           size /*:number*/) {
  let geometry = new BufferGeometry();
  const pos = new Float32Array(atom_arr.length * 3);
  for (let i = 0; i < atom_arr.length; i++) {
    const xyz = atom_arr[i].xyz;
    pos[3*i+0] = xyz[0];
    pos[3*i+1] = xyz[1];
    pos[3*i+2] = xyz[2];
  }
  geometry.addAttribute('position', new BufferAttribute(pos, 3));
  geometry.addAttribute('color', makeColorAttribute(color_arr));
  let material = new ShaderMaterial({
    uniforms: makeUniforms({size: size}),
    vertexShader: wheel_vert,
    fragmentShader: wheel_frag,
    fog: true,
    type: 'um_wheel',
  });
  let obj = new Points(geometry, material);
  // currently we use only lines for picking
  obj.raycast = function () {};
  return obj;
}

// For the ball-and-stick rendering we use so-called imposters.
// This technique was described in:
// http://doi.ieeecomputersociety.org/10.1109/TVCG.2006.115
// free copy here:
// http://vcg.isti.cnr.it/Publications/2006/TCM06/Tarini_FinalVersionElec.pdf
// and was nicely summarized in:
// http://www.sunsetlakesoftware.com/2011/05/08/enhancing-molecules-using-opengl-es-20

const sphere_vert = `
attribute vec3 color;
attribute vec2 corner;
uniform float radius;
varying vec3 vcolor;
varying vec2 vcorner;
varying vec3 vpos;

void main() {
  vcolor = color;
  vcorner = corner;
  vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
  vpos = mvPosition.xyz;
  mvPosition.xy += corner * radius;
  gl_Position = projectionMatrix * mvPosition;
}
`;

// based on 3Dmol imposter shaders
const sphere_frag = `
${fog_pars_fragment}
uniform mat4 projectionMatrix;
uniform vec3 lightDir;
uniform float radius;
varying vec3 vcolor;
varying vec2 vcorner;
varying vec3 vpos;

void main() {
  float sq = dot(vcorner, vcorner);
  if (sq > 1.0) discard;
  float z = sqrt(1.0-sq);
  vec3 xyz = vec3(vcorner.x, vcorner.y, z);
  vec4 projPos = projectionMatrix * vec4(vpos + radius * xyz, 1.0);
  gl_FragDepthEXT = 0.5 * ((gl_DepthRange.diff * (projPos.z / projPos.w)) +
                           gl_DepthRange.near + gl_DepthRange.far);
  float weight = clamp(dot(xyz, lightDir), 0.0, 1.0) * 0.8 + 0.2;
  gl_FragColor = vec4(weight * vcolor, 1.0);
  ${fog_end_fragment}
}
`;

const stick_vert = `
attribute vec3 color;
attribute vec3 axis;
attribute vec2 corner;
uniform float radius;
varying vec3 vcolor;
varying vec2 vcorner;
varying vec3 vpos;
varying vec3 vaxis;

void main() {
  vcolor = color;
  vcorner = corner;
  vaxis = normalize((modelViewMatrix * vec4(axis, 0.0)).xyz);
  vec2 normal = normalize(vec2(-vaxis.y, vaxis.x));
  vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
  vpos = mvPosition.xyz;
  mvPosition.xy += corner[1] * radius * normal;
  gl_Position = projectionMatrix * mvPosition;
}`;

const stick_frag = `
${fog_pars_fragment}
uniform mat4 projectionMatrix;
uniform vec3 lightDir;
uniform float radius;
varying vec3 vcolor;
varying vec2 vcorner;
varying vec3 vpos;
varying vec3 vaxis;
void main() {
  float central = 1.0 - vcorner[1] * vcorner[1];
  vec4 pos = vec4(vpos, 1.0);
  pos.z += radius * vaxis.z * central;
  vec4 projPos = projectionMatrix * pos;
  gl_FragDepthEXT = 0.5 * ((gl_DepthRange.diff * (projPos.z / projPos.w)) +
                           gl_DepthRange.near + gl_DepthRange.far);
  float weight = length(cross(vaxis, lightDir)) * central * 0.8 + 0.2;
  gl_FragColor = vec4(min(weight, 1.0) * vcolor, 1.0);
${fog_end_fragment}
}`;

export function makeSticks(vertex_arr /*:Num3[]*/,
                           color_arr /*:Color[]*/,
                           radius /*:number*/) {
  let uniforms = makeUniforms({
    radius: radius,
    lightDir: light_dir,
  });
  let material = new ShaderMaterial({
    uniforms: uniforms,
    vertexShader: stick_vert,
    fragmentShader: stick_frag,
    fog: true,
    type: 'um_stick',
  });
  material.extensions.fragDepth = true;

  const len = vertex_arr.length;
  const pos = double_pos(vertex_arr);
  const position = new Float32Array(pos);
  let axis = new Float32Array(6*len);
  for (let i = 0; i < 6 * len; i += 12) {
    for (let j = 0; j < 6; j++) axis[i+j] = pos[i+j+6] - pos[i+j];
    for (let j = 0; j < 6; j++) axis[i+j+6] = axis[i+j];
  }
  let geometry = new BufferGeometry();
  geometry.addAttribute('position', new BufferAttribute(position, 3));
  let corner = new Float32Array(4*len);
  for (let i = 0; 2 * i < len; i++) {
    corner[8*i + 0] = -1;  // 0
    corner[8*i + 1] = -1;  // 0
    corner[8*i + 2] = -1;  // 1
    corner[8*i + 3] = +1;  // 1
    corner[8*i + 4] = +1;  // 2
    corner[8*i + 5] = +1;  // 2
    corner[8*i + 6] = +1;  // 3
    corner[8*i + 7] = -1;  // 3
  }
  geometry.addAttribute('axis', new BufferAttribute(axis, 3));
  geometry.addAttribute('corner', new BufferAttribute(corner, 2));
  const color = double_color(color_arr);
  geometry.addAttribute('color', new BufferAttribute(color, 3));
  geometry.setIndex(make_quad_index_buffer(len/2));

  let mesh = new Mesh(geometry, material);
  mesh.raycast = line_raycast;
  return mesh;
}

export function makeBalls(atom_arr /*:AtomT[]*/,
                          color_arr /*:Color[]*/,
                          radius /*:number*/) {
  const N = atom_arr.length;
  let geometry = new BufferGeometry();

  let pos = new Float32Array(N * 4 * 3);
  for (let i = 0; i < N; i++) {
    const xyz = atom_arr[i].xyz;
    for (let j = 0; j < 4; j++) {
      for (let k = 0; k < 3; k++) {
        pos[3 * (4*i + j) + k] = xyz[k];
      }
    }
  }
  geometry.addAttribute('position', new BufferAttribute(pos, 3));

  let corner = new Float32Array(N * 4 * 2);
  for (let i = 0; i < N; i++) {
    corner[8*i + 0] = -1;  // 0
    corner[8*i + 1] = -1;  // 0
    corner[8*i + 2] = -1;  // 1
    corner[8*i + 3] = +1;  // 1
    corner[8*i + 4] = +1;  // 2
    corner[8*i + 5] = +1;  // 2
    corner[8*i + 6] = +1;  // 3
    corner[8*i + 7] = -1;  // 3
  }
  geometry.addAttribute('corner', new BufferAttribute(corner, 2));

  let colors = new Float32Array(N * 4 * 3);
  for (let i = 0; i < N; i++) {
    const col = color_arr[i];
    for (let j = 0; j < 4; j++) {
      colors[3 * (4*i + j) + 0] = col.r;
      colors[3 * (4*i + j) + 1] = col.g;
      colors[3 * (4*i + j) + 2] = col.b;
    }
  }
  geometry.addAttribute('color', new BufferAttribute(colors, 3));

  geometry.setIndex(make_quad_index_buffer(N));

  let material = new ShaderMaterial({
    uniforms: makeUniforms({
      radius: radius,
      lightDir: light_dir,
    }),
    vertexShader: sphere_vert,
    fragmentShader: sphere_frag,
    fog: true,
    type: 'um_sphere',
  });
  material.extensions.fragDepth = true;
  let obj = new Mesh(geometry, material);
  // currently we use only lines for picking
  obj.raycast = function () {};
  return obj;
}

// based on Line.prototype.raycast(), but skipping duplicated points
let inverseMatrix = new Matrix4();
let ray = new Ray();
// this function will be put on prototype
/* eslint-disable no-invalid-this */
function line_raycast(raycaster, intersects) {
  const precisionSq = raycaster.linePrecision * raycaster.linePrecision;
  inverseMatrix.getInverse(this.matrixWorld);
  ray.copy(raycaster.ray).applyMatrix4(inverseMatrix);
  let vStart = new Vector3();
  let vEnd = new Vector3();
  let interSegment = new Vector3();
  let interRay = new Vector3();
  const step = this.drawMode === TriangleStripDrawMode ? 1 : 2;
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

const label_vert = `
attribute vec2 uvs;
uniform vec2 canvas_size;
uniform vec2 win_size;
uniform float z_shift;
varying vec2 vUv;
void main() {
  vUv = uvs;
  vec2 rel_offset = vec2(0.02, -0.3);
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  gl_Position.xy += (uvs + rel_offset) * 2.0 * canvas_size / win_size;
  gl_Position.z += z_shift * projectionMatrix[2][2];
}`;

const label_frag = `
${fog_pars_fragment}
varying vec2 vUv;
uniform sampler2D map;
void main() {
  gl_FragColor = texture2D(map, vUv);
${fog_end_fragment}
}`;


export function makeLabel(text /*:string*/, options /*:{[key:string]: any}*/) {
  const canvas = makeCanvasWithText(text, options);
  if (!canvas) return;
  let texture = new Texture(canvas);
  texture.needsUpdate = true;

  // Rectangle geometry.
  let geometry = new BufferGeometry();
  const pos = options.pos;
  const position = new Float32Array([].concat(pos, pos, pos, pos));
  const uvs = new Float32Array([0, 1, 1, 1, 0, 0, 1, 0]);
  const indices = new Uint16Array([0, 2, 1, 2, 3, 1]);
  geometry.setIndex(new BufferAttribute(indices, 1));
  geometry.addAttribute('position', new BufferAttribute(position, 3));
  geometry.addAttribute('uvs', new BufferAttribute(uvs, 2));

  let material = new ShaderMaterial({
    uniforms: makeUniforms({map: texture,
                            canvas_size: [canvas.width, canvas.height],
                            win_size: options.win_size,
                            z_shift: options.z_shift}),
    vertexShader: label_vert,
    fragmentShader: label_frag,
    fog: true,
    type: 'um_label',
  });
  material.transparent = true;
  let mesh = new Mesh(geometry, material);
  mesh.remake = function (text, options) {
    texture.image = makeCanvasWithText(text, options);
    texture.needsUpdate = true;
  };
  return mesh;
}

// Add vertices of a 3d cross (representation of an unbonded atom)
export
function addXyzCross(vertices /*:Num3[]*/, xyz /*:Num3*/, r /*:number*/) {
  vertices.push([xyz[0]-r, xyz[1], xyz[2]], [xyz[0]+r, xyz[1], xyz[2]]);
  vertices.push([xyz[0], xyz[1]-r, xyz[2]], [xyz[0], xyz[1]+r, xyz[2]]);
  vertices.push([xyz[0], xyz[1], xyz[2]-r], [xyz[0], xyz[1], xyz[2]+r]);
}
