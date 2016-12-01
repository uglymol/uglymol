// @flow

import * as THREE from 'three';
import { makeLineMaterial, makeLineSegments, makeLine, makeRibbon,
         makeChickenWire, makeGrid, makeWheels, makeCube,
         makeRgbBox, makeLabel, addXyzCross } from './lines.js';
import { ElMap } from './elmap.js';
import { Model } from './model.js';


const ColorSchemes = [ // Viewer.prototype.ColorSchemes
  { // generally mimicks Coot
    name: 'coot dark',
    bg: 0x000000,
    fg: 0xFFFFFF,
    map_den: 0x3362B2,
    map_pos: 0x298029,
    map_neg: 0x8B2E2E,
    center: 0xC997B0,
    // atoms
    H: 0x858585, // H is normally invisible
    // C, N and O are taken approximately (by color-picker) from coot
    C: 0xb3b300,
    N: 0x7EAAFB,
    O: 0xF24984,
    S: 0x40ff40, // S in coot is too similar to C, here it is greener
    // Coot doesn't define other colors (?)
    MG: 0xc0c0c0,
    P: 0xffc040,
    CL: 0xa0ff60,
    CA: 0xffffff,
    MN: 0xff90c0,
    FE: 0xa03000,
    NI: 0x00ff80,
    def: 0xa0a0a0, // default atom color
  },
  // scheme made of "solarized" colors (http://ethanschoonover.com/solarized):
  // base03  base02  base01  base00  base0   base1   base2   base3
  // #002b36 #073642 #586e75 #657b83 #839496 #93a1a1 #eee8d5 #fdf6e3
  // yellow  orange  red     magenta violet  blue    cyan    green
  // #b58900 #cb4b16 #dc322f #d33682 #6c71c4 #268bd2 #2aa198 #859900
  {
    name: 'solarized dark',
    bg: 0x002b36,
    fg: 0xfdf6e3,
    map_den: 0x268bd2,
    map_pos: 0x859900,
    map_neg: 0xd33682,
    center: 0xfdf6e3,
    H: 0x586e75,
    C: 0x93a1a1,
    N: 0x6c71c4,
    O: 0xcb4b16,
    S: 0xb58900,
    def: 0xeee8d5,
  },
  {
    name: 'solarized light',
    bg: 0xfdf6e3,
    fg: 0x002b36,
    map_den: 0x268bd2,
    map_pos: 0x859900,
    map_neg: 0xd33682,
    center: 0x002b36,
    H: 0x93a1a1,
    C: 0x586e75,
    N: 0x6c71c4,
    O: 0xcb4b16,
    S: 0xb58900,
    def: 0x073642,
  },
  { // like in Coot after Edit > Background Color > White
    name: 'coot light',
    bg: 0xFFFFFF,
    fg: 0x000000,
    map_den: 0x3362B2,
    map_pos: 0x298029,
    map_neg: 0x8B2E2E,
    center: 0xC7C769,
    H: 0x999999,
    C: 0xA96464,
    N: 0x1C51B3,
    O: 0xC33869,
    S: 0x9E7B3D,
    def: 0x808080,
  },
];

const auto_speed = 1.0;

// map 2d position to sphere with radius 1.
function project_on_ball(x, y) {
  let z = 0;
  const length_sq = x * x + y * y;
  if (length_sq < 1) {  // in ellipse
    z = Math.sqrt(1.0 - length_sq);
  } else {  // in a corner
    const length = Math.sqrt(length_sq);
    x /= length;
    y /= length;
  }
  return [x, y, z];  // guaranteed to be normalized
}

function scale_by_height(value, size) { // for scaling bond_line
  return value * size[1] / 700;
}

const STATE = {NONE: -1, ROTATE: 0, PAN: 1, ZOOM: 2, PAN_ZOOM: 3, SLAB: 4,
               ROLL: 5, AUTO_ROTATE: 6, GO: 7};


// based on three.js/examples/js/controls/OrthographicTrackballControls.js
function Controls(camera, target) {
  let _state = STATE.NONE;
  let _rotate_start = new THREE.Vector3();
  let _rotate_end = new THREE.Vector3();
  let _zoom_start = new THREE.Vector2();
  let _zoom_end = new THREE.Vector2();
  let _pinch_start = 0;
  let _pinch_end = 0;
  let _pan_start = new THREE.Vector2();
  let _pan_end = new THREE.Vector2();
  let _panned = true;
  let _rotating = null;
  let _auto_stamp = null;
  let _go_func = null;

  // the far plane is more distant from the target than the near plane (3:1)
  this.slab_width = [2.5, 7.5];

  function rotate_camera(eye) {
    let quat = new THREE.Quaternion();
    quat.setFromUnitVectors(_rotate_end, _rotate_start);
    eye.applyQuaternion(quat);
    camera.up.applyQuaternion(quat);
    _rotate_end.applyQuaternion(quat);
    _rotate_start.copy(_rotate_end);
  }

  function zoom_camera(eye) {
    const dx = _zoom_end.x - _zoom_start.x;
    const dy = _zoom_end.y - _zoom_start.y;
    if (_state === STATE.ZOOM) {
      camera.zoom /= (1 - dx + dy);
    } else if (_state === STATE.SLAB) {
      target.addScaledVector(eye, -5.0 / eye.length() * dy);
    } else if (_state === STATE.ROLL) {
      camera.up.applyAxisAngle(eye, 0.05 * (dx - dy));
    }
    _zoom_start.copy(_zoom_end);
    return _state === STATE.SLAB ? 10*dx : null;
  }

  function pan_camera(eye) {
    let dx = _pan_end.x - _pan_start.x;
    let dy = _pan_end.y - _pan_start.y;
    dx *= 0.5 * (camera.right - camera.left) / camera.zoom;
    dy *= 0.5 * (camera.bottom - camera.top) / camera.zoom;
    let pan = eye.clone().cross(camera.up).setLength(dx);
    pan.addScaledVector(camera.up, dy / camera.up.length());
    camera.position.add(pan);
    target.add(pan);
    _pan_start.copy(_pan_end);
  }

  this.toggle_auto = function (param) {
    if (_state === STATE.AUTO_ROTATE && typeof param === typeof _rotating) {
      _state = STATE.NONE;
    } else {
      _state = STATE.AUTO_ROTATE;
      _auto_stamp = null;
      _rotating = param;
    }
  };

  this.is_going = function () { return _state === STATE.GO; };

  this.is_moving = function () { return _state !== STATE.NONE; };

  function auto_rotate(eye) {
    _rotate_start.copy(eye).normalize();
    const now = Date.now();
    const elapsed = (_auto_stamp !== null ? now - _auto_stamp : 16.7);
    let speed = 1.8e-5 * elapsed * auto_speed;
    _auto_stamp = now;
    if (_rotating === true) {
      speed = -speed;
    } else if (_rotating !== false) {
      _rotating += 0.02;
      speed = 4e-5 * auto_speed * Math.cos(_rotating);
    }
    _rotate_end.crossVectors(camera.up, eye).multiplyScalar(speed)
      .add(_rotate_start);
  }

  this.update = function () {
    let changed = false;
    let eye = camera.position.clone().sub(target);
    if (_state === STATE.AUTO_ROTATE) {
      auto_rotate(eye);
    }
    if (!_rotate_start.equals(_rotate_end)) {
      rotate_camera(eye);
      changed = true;
    }
    if (_pinch_end !== _pinch_start) {
      camera.zoom *= _pinch_end / _pinch_start;
      _pinch_start = _pinch_end;
      changed = true;
    }
    if (!_zoom_end.equals(_zoom_start)) {
      const dslab = zoom_camera(eye);
      if (dslab) {
        this.slab_width[0] = Math.max(this.slab_width[0] + dslab, 0.01);
        this.slab_width[1] = Math.max(this.slab_width[1] + dslab, 0.01);
      }
      changed = true;
    }
    if (!_pan_end.equals(_pan_start)) {
      pan_camera(eye);
      _panned = true;
      changed = true;
    }
    camera.position.addVectors(target, eye);
    if (_state === STATE.GO && _go_func) {
      _go_func();
      changed = true;
    }
    camera.lookAt(target);
    return changed;
  };

  this.start = function (new_state, x, y, dist) {
    if (_state === STATE.NONE || _state === STATE.AUTO_ROTATE) {
      _state = new_state;
    }
    this.move(x, y, dist);
    switch (_state) {
      case STATE.ROTATE:
        _rotate_start.copy(_rotate_end);
        break;
      case STATE.ZOOM:
      case STATE.SLAB:
      case STATE.ROLL:
        _zoom_start.copy(_zoom_end);
        break;
      case STATE.PAN:
        _pan_start.copy(_pan_end);
        _panned = false;
        break;
      case STATE.PAN_ZOOM:
        _pinch_start = _pinch_end;
        _pan_start.copy(_pan_end);
        break;
    }
  };

  this.move = function (x, y, dist) {
    switch (_state) {
      case STATE.ROTATE: {
        const xyz = project_on_ball(x, y);
        //console.log(camera.projectionMatrix);
        //console.log(camera.matrixWorld);
        // TODO maybe use project()/unproject()/applyProjection()
        const eye = camera.position.clone().sub(target);
        _rotate_end.crossVectors(camera.up, eye).setLength(xyz[0]);
        _rotate_end.addScaledVector(camera.up, xyz[1] / camera.up.length());
        _rotate_end.addScaledVector(eye, xyz[2] / eye.length());
        break;
      }
      case STATE.ZOOM:
      case STATE.SLAB:
      case STATE.ROLL:
        _zoom_end.set(x, y);
        break;
      case STATE.PAN:
        _pan_end.set(x, y);
        break;
      case STATE.PAN_ZOOM:
        _pan_end.set(x, y);
        _pinch_end = dist;
        break;
    }
  };

  this.stop = function () {
    let ret = null;
    if (_state === STATE.PAN && !_panned) ret = _pan_start;
    _state = STATE.NONE;
    _rotate_start.copy(_rotate_end);
    _pinch_start = _pinch_end;
    _pan_start.copy(_pan_end);
    return ret;
  };

  this.go_to = function (targ, cam_pos, cam_up, steps) {
    if (targ instanceof Array) {
      targ = new THREE.Vector3(targ[0], targ[1], targ[2]);
    }
    if ((!targ || targ.distanceToSquared(target) < 0.1) &&
        (!cam_pos || cam_pos.distanceToSquared(camera.position) < 0.1) &&
        (!cam_up || cam_up.distanceToSquared(camera.up) < 0.1)) {
      return;
    }
    _state = STATE.GO;
    steps = steps || (60 / auto_speed);
    let alphas = [];
    let prev_pos = 0;
    for (let i = 1; i <= steps; ++i) {
      let pos = i / steps;
      // quadratic easing
      pos = pos < 0.5 ? 2 * pos * pos : -2 * pos * (pos-2) - 1;
      alphas.push((pos - prev_pos) / (1 - prev_pos));
      prev_pos = pos;
    }
    _go_func = function () {
      const a = alphas.shift();
      if (targ) {
        // unspecified cam_pos - camera stays in the same distance to target
        if (!cam_pos) camera.position.sub(target);
        target.lerp(targ, a);
        if (!cam_pos) camera.position.add(target);
      }
      if (cam_pos) camera.position.lerp(cam_pos, a);
      if (cam_up) camera.up.lerp(cam_up, a);
      if (alphas.length === 0) {
        _state = STATE.NONE;
        _go_func = null;
      }
    };
  };
}


// options handled by select_next()

const COLOR_AIMS = ['element', 'B-factor', 'occupancy', 'index', 'chain'];
const RENDER_STYLES = ['lines', 'trace', 'ribbon'/*, 'ball&stick'*/];
const MAP_STYLES = ['marching cubes', 'squarish'/*, 'snapped MC'*/];
const LINE_STYLES = ['normal', 'simplistic'];
const LABEL_FONTS = ['bold 14px', '14px', '16px', 'bold 16px'];

function rainbow_value(v, vmin, vmax) {
  let c = new THREE.Color(0xe0e0e0);
  if (vmin < vmax) {
    const ratio = (v - vmin) / (vmax - vmin);
    const hue = (240 - (240 * ratio)) / 360;
    c.setHSL(hue, 1.0, 0.5);
  }
  return c;
}

function color_by(style, atoms, elem_colors) {
  let color_func;
  let i;
  const last_atom = atoms[atoms.length-1];
  if (style === 'index') {
    color_func = function (atom) {
      return rainbow_value(atom.i_seq, 0, last_atom.i_seq);
    };
  } else if (style === 'B-factor') {
    let vmin = Infinity;
    let vmax = -Infinity;
    for (i = 0; i < atoms.length; i++) {
      const v = atoms[i].b;
      if (v > vmax) vmax = v;
      if (v < vmin) vmin = v;
    }
    //console.log('B-factors in [' + vmin + ', ' + vmax + ']');
    color_func = function (atom) {
      return rainbow_value(atom.b, vmin, vmax);
    };
  } else if (style === 'occupancy') {
    color_func = function (atom) {
      return rainbow_value(atom.occ, 0, 1);
    };
  } else if (style === 'chain') {
    color_func = function (atom) {
      return rainbow_value(atom.chain_index, 0, last_atom.chain_index);
    };
  } else { // element
    color_func = function (atom) {
      return elem_colors[atom.element] || elem_colors.def;
    };
  }
  let colors = [];
  for (i = 0; i < atoms.length; i++) {
    colors.push(color_func(atoms[i]));
  }
  return colors;
}

function MapBag(map, is_diff_map) {
  this.map = map;
  this.name = '';
  this.isolevel = is_diff_map ? 3.0 : 1.5;
  this.visible = true;
  this.types = is_diff_map ? ['map_pos', 'map_neg'] : ['map_den'];
  this.block_ctr = new THREE.Vector3(Infinity, 0, 0);
  this.el_objects = []; // three.js objects
}


function ModelBag(model, config, win_size) {
  this.model = model;
  this.name = '';
  this.visible = true;
  this.conf = config;
  this.win_size = win_size;
  this.atomic_objects = null; // list of three.js objects
}

ModelBag.prototype.get_visible_atoms = function () {
  const atoms = this.model.atoms;
  if (this.conf.hydrogens || !this.model.has_hydrogens) {
    return atoms;
  }
  const non_h = [];
  for (let i = 0; i < atoms.length; i++) {
    if (atoms[i].element !== 'H') {
      non_h.push(atoms[i]);
    }
  }
  return non_h;
};

ModelBag.prototype.add_bonds = function (ligands_only, ball_size) {
  const visible_atoms = this.get_visible_atoms();
  const color_style = ligands_only ? 'element' : this.conf.color_aim;
  const colors = color_by(color_style, visible_atoms, this.conf.colors);
  let vertex_arr /*:THREE.Vector3[]*/ = [];
  let color_arr = [];
  const opt = { hydrogens: this.conf.hydrogens,
                ligands_only: ligands_only,
                balls: this.conf.render_style === 'ball&stick' };
  for (let i = 0; i < visible_atoms.length; i++) {
    const atom = visible_atoms[i];
    const color = colors[i];
    if (ligands_only && !atom.is_ligand) continue;
    if (atom.bonds.length === 0 && !opt.balls) { // nonbonded, draw star
      addXyzCross(vertex_arr, atom.xyz, 0.7);
      for (let n = 0; n < 6; n++) {
        color_arr.push(color);
      }
    } else { // bonded, draw lines
      for (let j = 0; j < atom.bonds.length; j++) {
        const other = this.model.atoms[atom.bonds[j]];
        if (!opt.hydrogens && other.element === 'H') continue;
        // Coot show X-H bonds as thinner lines in a single color.
        // Here we keep it simple and render such bonds like all others.
        if (opt.ligands_only && !other.is_ligand) continue;
        const mid = atom.midpoint(other);
        const vmid = new THREE.Vector3(mid[0], mid[1], mid[2]);
        const vatom = new THREE.Vector3(atom.xyz[0], atom.xyz[1], atom.xyz[2]);
        if (opt.balls) {
          const lerp_factor = vatom.distanceTo(vmid) / ball_size;
          vatom.lerp(vmid, lerp_factor);
        }
        vertex_arr.push(vatom, vmid);
        color_arr.push(color, color);
      }
    }
  }
  //console.log('add_bonds() vertex count: ' + vertex_arr.length);
  const linewidth = scale_by_height(this.conf.bond_line, this.win_size);
  const use_gl_lines = this.conf.line_style === 'simplistic';
  const material = makeLineMaterial({
    gl_lines: use_gl_lines,
    linewidth: linewidth,
    win_size: this.win_size,
    segments: true,
  });
  this.atomic_objects.push(makeLineSegments(material, vertex_arr, color_arr));
  if (opt.balls) {
    this.atomic_objects.push(makeWheels(visible_atoms, colors, ball_size));
  } else if (!use_gl_lines && !ligands_only) {
    // wheels (discs) as simplistic round caps
    this.atomic_objects.push(makeWheels(visible_atoms, colors, linewidth));
  }
};

ModelBag.prototype.add_trace = function () {
  const segments = this.model.extract_trace();
  const visible_atoms = [].concat.apply([], segments);
  const colors = color_by(this.conf.color_aim, visible_atoms, this.conf.colors);
  const material = makeLineMaterial({
    gl_lines: this.conf.line_style === 'simplistic',
    linewidth: scale_by_height(this.conf.bond_line, this.win_size),
    win_size: this.win_size,
  });
  let k = 0;
  for (let i = 0; i < segments.length; i++) {
    const seg = segments[i];
    const color_slice = colors.slice(k, k + seg.length);
    k += seg.length;
    const line = makeLine(material, seg, color_slice);
    this.atomic_objects.push(line);
  }
};

ModelBag.prototype.add_ribbon = function (smoothness) {
  const segments = this.model.extract_trace();
  const res_map = this.model.get_residues();
  const visible_atoms = [].concat.apply([], segments);
  const colors = color_by(this.conf.color_aim, visible_atoms, this.conf.colors);
  let k = 0;
  for (let i = 0; i < segments.length; i++) {
    const seg = segments[i];
    let tangents = [];
    let last = [0, 0, 0];
    for (let j = 0; j < seg.length; j++) {
      const residue = res_map[seg[j].resid()];
      const tang = this.model.calculate_tangent_vector(residue);
      // untwisting (usually applies to beta-strands)
      if (tang[0]*last[0] + tang[1]*last[1] + tang[2]*last[2] < 0) {
        tang[0] = -tang[0];
        tang[1] = -tang[1];
        tang[2] = -tang[2];
      }
      tangents.push(tang);
      last = tang;
    }
    const color_slice = colors.slice(k, k + seg.length);
    k += seg.length;
    const obj = makeRibbon(seg, color_slice, tangents, smoothness);
    this.atomic_objects.push(obj);
  }
};

export function Viewer(options /*: {[key: string]: any}*/) {
  // rendered objects
  this.model_bags = [];
  this.map_bags = [];
  this.decor = {cell_box: null, selection: null, zoom_grid: makeGrid() };
  this.labels = {};
  this.nav = null;

  this.config = {
    bond_line: 4.0, // ~ to height, like in Coot (see scale_by_height())
    map_line: 1.25,  // for any height
    map_radius: 10.0,
    map_style: MAP_STYLES[0],
    render_style: RENDER_STYLES[0],
    color_aim: COLOR_AIMS[0],
    line_style: LINE_STYLES[0],
    label_font: LABEL_FONTS[0],
    colors: this.ColorSchemes[0],
    hydrogens: false,
  };
  this.set_colors();
  this.window_size = [1, 1]; // it will be set in resize()
  this.window_offset = [0, 0];

  this.last_ctr = new THREE.Vector3(Infinity, 0, 0);
  this.selected_atom = null;
  this.active_model_bag = null;
  this.scene = new THREE.Scene();
  this.scene.fog = new THREE.Fog(this.config.colors.bg, 0, 1);
  this.light = new THREE.AmbientLight(0xffffff);
  this.scene.add(this.light);
  if (options.share_view) {
    this.target = options.share_view.target;
    this.camera = options.share_view.camera;
    this.controls = options.share_view.controls;
    this.tied_viewer = options.share_view;
    this.tied_viewer.tied_viewer = this; // not GC friendly
  } else {
    this.target = new THREE.Vector3();
    this.camera = new THREE.OrthographicCamera();
    this.controls = new Controls(this.camera, this.target);
  }
  this.raycaster = new THREE.Raycaster();
  this.default_camera_pos = [0, 0, 100];
  this.set_common_key_bindings();
  if (this.constructor === Viewer) this.set_real_space_key_bindings();
  if (typeof document === 'undefined') return;  // for testing on node

  try {
    this.renderer = new THREE.WebGLRenderer({antialias: true});
  } catch (e) {
    this.hud('no WebGL in your browser?', 'ERR');
    this.renderer = null;
    return;
  }

  function get_elem(name) {
    if (options[name] === null) return null;
    return document.getElementById(options[name] || name);
  }
  this.container = get_elem('viewer');
  this.hud_el = get_elem('hud');
  this.help_el = get_elem('help');
  if (this.hud_el) {
    this.initial_hud_html = this.hud_el.innerHTML;
    this.initial_hud_bg = this.hud_el.style['background-color'];
  }

  if (this.container === null) return; // can be null in headless tests
  this.renderer.setClearColor(this.config.colors.bg, 1);
  this.renderer.setPixelRatio(window.devicePixelRatio);
  this.resize();
  this.camera.zoom = this.camera.right / 35.0;  // arbitrary choice
  this.container.appendChild(this.renderer.domElement);
  if (options.focusable) {
    this.renderer.domElement.tabIndex = 0;
  }
  this.decor.zoom_grid.visible = false;
  this.scene.add(this.decor.zoom_grid);
  if (window.Stats) { // set by including three/examples/js/libs/stats.min.js
    this.stats = new window.Stats();
    this.container.appendChild(this.stats.dom);
  }

  window.addEventListener('resize', this.resize.bind(this));
  let el = this.renderer.domElement;
  let keydown_el = (options.focusable ? el : window);
  keydown_el.addEventListener('keydown', this.keydown.bind(this));
  el.addEventListener('contextmenu', function (e) { e.preventDefault(); });
  el.addEventListener('mousewheel', this.mousewheel.bind(this));
  el.addEventListener('MozMousePixelScroll', this.mousewheel.bind(this));
  el.addEventListener('mousedown', this.mousedown.bind(this));
  el.addEventListener('touchstart', this.touchstart.bind(this));
  el.addEventListener('touchmove', this.touchmove.bind(this));
  el.addEventListener('touchend', this.touchend.bind(this));
  el.addEventListener('touchcancel', this.touchend.bind(this));
  el.addEventListener('dblclick', this.dblclick.bind(this));

  let self = this;

  this.mousemove = function (event) {
    event.preventDefault();
    //event.stopPropagation();
    self.controls.move(self.relX(event), self.relY(event));
  };

  this.mouseup = function (event) {
    event.preventDefault();
    event.stopPropagation();
    document.removeEventListener('mousemove', self.mousemove);
    document.removeEventListener('mouseup', self.mouseup);
    self.decor.zoom_grid.visible = false;
    const not_panned = self.controls.stop();
    // special case - centering on atoms after action 'pan' with no shift
    if (not_panned) {
      const atom = self.pick_atom(not_panned, self.camera);
      if (atom != null) {
        self.select_atom(atom, {steps: 60 / auto_speed});
      }
    }
    self.redraw_maps();
  };

  this.scheduled = false;
  this.request_render();
}

Viewer.prototype.pick_atom = function (coords, camera) {
  const bag = this.active_model_bag;
  if (bag === null) return;
  this.raycaster.setFromCamera(coords, camera);
  this.raycaster.near = camera.near;
  // '0.15' b/c the furthest 15% is hardly visible in the fog
  this.raycaster.far = camera.far - 0.15 * (camera.far - camera.near);
  this.raycaster.linePrecision = 0.3;
  let intersects = this.raycaster.intersectObjects(bag.atomic_objects);
  if (intersects.length < 1) return null;
  intersects.sort(function (x) { return x.line_dist || Infinity; });
  const p = intersects[0].point;
  return bag.model.get_nearest_atom(p.x, p.y, p.z);
};

Viewer.prototype.set_colors = function (scheme) {
  function to_col(x) { return new THREE.Color(x); }
  if (scheme == null) {
    scheme = this.ColorSchemes[0];
  } else if (typeof scheme === 'number') {
    scheme = this.ColorSchemes[scheme % this.ColorSchemes.length];
  } else if (typeof scheme === 'string') {
    for (let i = 0; i !== this.ColorSchemes.length; i++) {
      if (this.ColorSchemes[i].name === scheme) {
        scheme = this.ColorSchemes[i];
        break;
      }
    }
  }
  if (scheme.bg === undefined) return;
  if (typeof scheme.bg === 'number') {
    for (let key in scheme) {
      if (key !== 'name') {
        scheme[key] = scheme[key] instanceof Array ? scheme[key].map(to_col)
                                                   : to_col(scheme[key]);
      }
    }
  }
  this.decor.zoom_grid.color_value.set(scheme.fg);
  this.redraw_all();
};

// relative position on canvas in normalized device coordinates [-1, +1]
Viewer.prototype.relX = function (evt) {
  return 2 * (evt.pageX - this.window_offset[0]) / this.window_size[0] - 1;
};

Viewer.prototype.relY = function (evt) {
  return 1 - 2 * (evt.pageY - this.window_offset[1]) / this.window_size[1];
};

Viewer.prototype.hud = function (text, type) {
  if (typeof document === 'undefined') return;  // for testing on node
  let el = this.hud_el;
  if (el) {
    if (text !== undefined) {
      if (type === 'HTML') {
        el.innerHTML = text;
      } else {
        el.textContent = text;
      }
    } else {
      el.innerHTML = this.initial_hud_html;
    }
    const err = (type === 'ERR');
    el.style['background-color'] = (err ? '#b00' : this.initial_hud_bg);
    if (err) console.log('ERR: ' + text);
  } else {
    console.log('hud: ' + text);
  }
};

Viewer.prototype.redraw_center = function () {
  if (this.target.distanceToSquared(this.last_ctr) > 0.0001) {
    this.last_ctr.copy(this.target);
    if (this.mark) {
      this.scene.remove(this.mark);
    }
    this.mark = makeCube(0.1, this.target, {
      color: this.config.colors.center,
      linewidth: 2,
      win_size: this.window_size,
    });
    this.scene.add(this.mark);
  }
};

Viewer.prototype.redraw_maps = function (force) {
  this.redraw_center();
  for (let i = 0; i < this.map_bags.length; i++) {
    const map_bag = this.map_bags[i];
    if (force || this.target.distanceToSquared(map_bag.block_ctr) > 0.01) {
      this.redraw_map(map_bag);
    }
  }
};

Viewer.prototype.remove_and_dispose = function (obj, only_dispose) {
  if (!only_dispose) this.scene.remove(obj);
  if (obj.geometry) obj.geometry.dispose();
  if (obj.material) {
    if (obj.material.uniforms && obj.material.uniforms.map) {
      obj.material.uniforms.map.value.dispose();
    }
    obj.material.dispose();
  }
  for (let i = 0; i < obj.children.length; i++) {
    this.remove_and_dispose(obj.children[i]);
  }
};

Viewer.prototype.clear_el_objects = function (map_bag) {
  for (let i = 0; i < map_bag.el_objects.length; i++) {
    this.remove_and_dispose(map_bag.el_objects[i]);
  }
  map_bag.el_objects = [];
};

Viewer.prototype.clear_atomic_objects = function (model) {
  if (model.atomic_objects) {
    for (let i = 0; i < model.atomic_objects.length; i++) {
      this.remove_and_dispose(model.atomic_objects[i]);
    }
  }
  model.atomic_objects = null;
};

Viewer.prototype.set_atomic_objects = function (model_bag) {
  model_bag.atomic_objects = [];
  switch (model_bag.conf.render_style) {
    case 'lines':
      model_bag.add_bonds();
      break;
    case 'ball&stick': {
      const h_scale = this.camera.projectionMatrix.elements[5];
      const ball_size = Math.max(1, 200 * h_scale);
      model_bag.add_bonds(false, ball_size);
      break;
    }
    case 'trace':  // + lines for ligands
      model_bag.add_trace();
      model_bag.add_bonds(true);
      break;
    case 'ribbon':
      model_bag.add_ribbon(8);
      model_bag.add_bonds(true);
      break;
  }
  for (let i = 0; i < model_bag.atomic_objects.length; i++) {
    this.scene.add(model_bag.atomic_objects[i]);
  }
};

// Add/remove label if `show` is specified, toggle otherwise.
Viewer.prototype.toggle_label = function (atom, show) {
  if (!atom) return;
  const text = atom.short_label();
  const uid = text; // we assume that the labels are unique - often true
  const is_shown = (uid in this.labels);
  if (show === undefined) show = !is_shown;
  if (show) {
    if (is_shown) return;
    const label = makeLabel(text, {
      pos: atom.xyz,
      font: this.config.label_font,
      color: '#' + this.config.colors.fg.getHexString(),
      win_size: this.window_size,
    });
    if (!label) return;
    this.labels[uid] = label;
    this.scene.add(label);
  } else {
    if (!is_shown) return;
    this.remove_and_dispose(this.labels[uid]);
    delete this.labels[uid];
  }
};

Viewer.prototype.redraw_labels = function () {
  for (let uid in this.labels) { // eslint-disable-line guard-for-in
    const text = uid;
    this.labels[uid].remake(text, {
      font: this.config.label_font,
      color: '#' + this.config.colors.fg.getHexString(),
    });
  }
};


Viewer.prototype.toggle_map_visibility = function (map_bag) {
  if (typeof map_bag === 'number') {
    map_bag = this.map_bags[map_bag];
  }
  map_bag.visible = !map_bag.visible;
  this.redraw_map(map_bag);
  this.request_render();
};

Viewer.prototype.redraw_map = function (map_bag) {
  this.clear_el_objects(map_bag);
  if (map_bag.visible) {
    map_bag.map.block = null;
    this.add_el_objects(map_bag);
  }
};

Viewer.prototype.toggle_model_visibility = function (model_bag) {
  model_bag = model_bag || this.active_model_bag;
  model_bag.visible = !model_bag.visible;
  this.redraw_model(model_bag);
  this.request_render();
};

Viewer.prototype.redraw_model = function (model_bag) {
  this.clear_atomic_objects(model_bag);
  if (model_bag.visible) {
    this.set_atomic_objects(model_bag);
  }
};

Viewer.prototype.redraw_models = function () {
  for (let i = 0; i < this.model_bags.length; i++) {
    this.redraw_model(this.model_bags[i]);
  }
};

Viewer.prototype.add_el_objects = function (map_bag) {
  if (!map_bag.visible || this.config.map_radius <= 0) return;
  if (!map_bag.map.block) {
    map_bag.block_ctr.copy(this.target);
    map_bag.map.extract_block(this.config.map_radius,
                              [this.target.x, this.target.y, this.target.z]);
  }
  for (let i = 0; i < map_bag.types.length; i++) {
    const mtype = map_bag.types[i];
    const isolevel = (mtype === 'map_neg' ? -1 : 1) * map_bag.isolevel;
    const iso = map_bag.map.isomesh_in_block(isolevel, this.config.map_style);

    const obj = makeChickenWire(iso, {
      color: this.config.colors[mtype],
      linewidth: this.config.map_line,
    });
    map_bag.el_objects.push(obj);
    this.scene.add(obj);
  }
};

Viewer.prototype.change_isolevel_by = function (map_idx, delta) {
  if (map_idx >= this.map_bags.length) return;
  const map_bag = this.map_bags[map_idx];
  map_bag.isolevel += delta;
  //TODO: move slow part into update()
  this.clear_el_objects(map_bag);
  this.add_el_objects(map_bag);
  const abs_level = map_bag.map.abs_level(map_bag.isolevel);
  let abs_text = abs_level.toFixed(4);
  let tied = this.tied_viewer;
  if (tied && map_idx < tied.map_bags.length) {
    let tied_bag = tied.map_bags[map_idx];
    // Should we tie by sigma or absolute level? Now it's sigma.
    tied_bag.isolevel = map_bag.isolevel;
    abs_text += ' / ' + tied_bag.map.abs_level(tied_bag.isolevel).toFixed(4);
    tied.clear_el_objects(tied_bag);
    tied.add_el_objects(tied_bag);
  }
  this.hud('map ' + (map_idx+1) + ' level =  ' + abs_text +
           ' e/\u212B\u00B3 (' + map_bag.isolevel.toFixed(2) + ' rmsd)');
};

Viewer.prototype.change_map_radius = function (delta) {
  const RMAX = 40;
  const cf = this.config;
  cf.map_radius = Math.min(Math.max(cf.map_radius + delta, 0), RMAX);
  let info = 'map "radius": ' + cf.map_radius;
  if (cf.map_radius === RMAX) info += ' (max)';
  else if (cf.map_radius === 0) info += ' (hidden maps)';
  this.hud(info);
  this.redraw_maps(true);
};

Viewer.prototype.change_slab_width_by = function (delta) {
  let slab_width = this.controls.slab_width;
  slab_width[0] = Math.max(slab_width[0] + delta, 0.01);
  slab_width[1] = Math.max(slab_width[1] + delta, 0.01);
  this.update_camera();
  this.hud('clip width: ' + (this.camera.far-this.camera.near).toPrecision(3));
};

Viewer.prototype.change_zoom_by_factor = function (mult) {
  this.camera.zoom *= mult;
  this.update_camera();
  this.hud('zoom: ' + this.camera.zoom.toPrecision(3));
};

Viewer.prototype.change_bond_line = function (delta) {
  this.config.bond_line = Math.max(this.config.bond_line + delta, 0.1);
  this.redraw_models();
  this.hud('bond width: ' + scale_by_height(this.config.bond_line,
                                            this.window_size).toFixed(1));
};

Viewer.prototype.change_map_line = function (delta) {
  this.config.map_line = Math.max(this.config.map_line + delta, 0.1);
  this.redraw_maps(true);
  this.hud('wireframe width: ' + this.config.map_line.toFixed(1));
};

Viewer.prototype.toggle_full_screen = function () {
  let d = document;
  if (d.fullscreenElement || d.mozFullScreenElement ||
      d.webkitFullscreenElement || d.msFullscreenElement) {
    let ex = d.exitFullscreen || d.webkitExitFullscreen ||
    // flow-ignore-line property `msExitFullscreen` not found in document
             d.mozCancelFullScreen || d.msExitFullscreen;
    // flow-ignore-line cannot call property `exitFullscreen` of unknown type
    if (ex) ex.call(d);
  } else {
    let el = this.container;
    let req = el.requestFullscreen || el.webkitRequestFullscreen ||
              el.mozRequestFullScreen || el.msRequestFullscreen;
    if (req) req.call(el);
  }
};

Viewer.prototype.toggle_cell_box = function () {
  if (this.decor.cell_box) {
    this.scene.remove(this.decor.cell_box);
    this.decor.cell_box = null;
  } else {
    let uc = null;
    if (this.model_bags.length > 0) {
      uc = this.model_bags[0].model.unit_cell;
    }
    // model may not have unit cell
    if (!uc && this.map_bags.length > 0) {
      uc = this.map_bags[0].map.unit_cell;
    }
    if (uc) {
      this.decor.cell_box = makeRgbBox(uc.orthogonalize.bind(uc), {
        color: this.config.colors.fg,
      });
      this.scene.add(this.decor.cell_box);
    }
  }
};

function vec3_to_fixed(vec, n) {
  return [vec.x.toFixed(n), vec.y.toFixed(n), vec.z.toFixed(n)];
}

Viewer.prototype.shift_clip = function (delta) {
  let eye = this.camera.position.clone().sub(this.target);
  eye.multiplyScalar(delta / eye.length());
  this.target.add(eye);
  this.camera.position.add(eye);
  this.update_camera();
  this.redraw_maps();
  this.hud('clip shifted by [' + vec3_to_fixed(eye, 2).join(' ') + ']');
};

Viewer.prototype.go_to_nearest_Ca = function () {
  const t = this.target;
  if (this.active_model_bag === null) return;
  const a = this.active_model_bag.model.get_nearest_atom(t.x, t.y, t.z, 'CA');
  if (a) {
    this.select_atom(a);
  } else {
    this.hud('no nearby CA');
  }
};

Viewer.prototype.permalink = function () {
  if (typeof window === 'undefined') return;
  window.location.hash = '#xyz=' + vec3_to_fixed(this.target, 1).join(',') +
    '&eye=' + vec3_to_fixed(this.camera.position, 1).join(',') +
    '&zoom=' + this.camera.zoom.toFixed(0);
  this.hud('copy URL from the location bar');
};

Viewer.prototype.redraw_all = function () {
  if (!this.renderer) return;
  this.scene.fog.color = this.config.colors.bg;
  if (this.renderer) this.renderer.setClearColor(this.config.colors.bg, 1);
  this.redraw_models();
  this.redraw_maps(true);
  this.redraw_labels();
};

Viewer.prototype.toggle_help = function () {
  let el = this.help_el;
  if (!el) return;
  el.style.display = el.style.display === 'block' ? 'none' : 'block';
  if (el.innerHTML === '') {
    el.innerHTML = [this.MOUSE_HELP, this.KEYBOARD_HELP,
                    this.ABOUT_HELP].join('\n\n');
  }
};

Viewer.prototype.MOUSE_HELP = [
  '<b>mouse:</b>',
  'Left = rotate',
  'Middle or Ctrl+Left = pan',
  'Right = zoom',
  'Ctrl+Right = clipping',
  'Ctrl+Shift+Right = roll',
  'Wheel = σ level',
  'Shift+Wheel = diff map σ',
].join('\n');

Viewer.prototype.KEYBOARD_HELP = [
  '<b>keyboard:</b>',
  'H = toggle help',
  'T = representation',
  'C = coloring',
  'B = bg color',
  'Q = label font',
  '+/- = sigma level',
  ']/[ = map radius',
  'D/F = clip width',
  'numpad 3/. = move clip',
  'M/N = zoom',
  'U = unitcell box',
  'Y = hydrogens',
  'R = center view',
  'W = wireframe style',
  'I = spin',
  'K = rock',
  'Home/End = bond width',
  '\\ = bond caps',
  'P = nearest Cα',
  'Shift+P = permalink',
  '(Shift+)space = next res.',
  'Shift+F = full screen',
].join('\n');

Viewer.prototype.ABOUT_HELP =
  '<a href="https://uglymol.github.io">about uglymol</a>';

Viewer.prototype.select_next = function (info, key, options, back) {
  const old_idx = options.indexOf(this.config[key]);
  const len = options.length;
  const new_idx = (old_idx + (back ? len - 1 : 1)) % len;
  this.config[key] = options[new_idx];
  let html = info + ':';
  for (let i = 0; i < len; i++) {
    const tag = (i === new_idx ? 'u' : 's');
    const opt_name = options[i].name || options[i];
    html += ' <' + tag + '>' + opt_name + '</' + tag + '>';
  }
  this.hud(html, 'HTML');
};

Viewer.prototype.keydown = function (evt) {
  const action = this.key_bindings[evt.keyCode];
  if (action) {
    (action.bind(this))(evt);
  } else {
    if (action === false) evt.preventDefault();
    if (this.help_el) this.hud('Nothing here. Press H for help.');
  }
  this.request_render();
};

Viewer.prototype.set_common_key_bindings = function () {
  let kb = new Array(256);
  // Home
  kb[36] = function (evt) {
    evt.ctrlKey ? this.change_map_line(0.1) : this.change_bond_line(0.2);
  };
  // End
  kb[35] = function (evt) {
    evt.ctrlKey ? this.change_map_line(-0.1) : this.change_bond_line(-0.2);
  };
  // b
  kb[66] = function (evt) {
    this.select_next('color scheme', 'colors', this.ColorSchemes, evt.shiftKey);
    this.set_colors(this.config.colors);
  };
  // c
  kb[67] = function (evt) {
    this.select_next('coloring by', 'color_aim', COLOR_AIMS, evt.shiftKey);
    this.redraw_models();
  };
  // d
  kb[68] = function () { this.change_slab_width_by(-0.1); };
  // f
  kb[70] = function (evt) {
    evt.shiftKey ? this.toggle_full_screen() : this.change_slab_width_by(0.1);
  };
  // h
  kb[72] = this.toggle_help;
  // i
  kb[73] = function (evt) {
    this.hud('toggled spinning');
    this.controls.toggle_auto(evt.shiftKey);
  };
  // k
  kb[75] = function () {
    this.hud('toggled rocking');
    this.controls.toggle_auto(0.0);
  };
  // m
  kb[77] = function (evt) {
    this.change_zoom_by_factor(evt.shiftKey ? 1.2 : 1.03);
  };
  // n
  kb[78] = function (evt) {
    this.change_zoom_by_factor(1 / (evt.shiftKey ? 1.2 : 1.03));
  };
  // q
  kb[81] = function (evt) {
    this.select_next('label font', 'label_font', LABEL_FONTS, evt.shiftKey);
    this.redraw_labels();
  };
  // r
  kb[82] = function (evt) {
    if (evt.shiftKey) {
      this.hud('redraw!');
      this.redraw_all();
    } else {
      this.hud('recentered');
      this.recenter();
    }
  };
  // u
  kb[85] = function () {
    this.hud('toggled unit cell box');
    this.toggle_cell_box();
  };
  // w
  kb[87] = function (evt) {
    this.select_next('map style', 'map_style', MAP_STYLES, evt.shiftKey);
    this.redraw_maps(true);
  };
  // add, equals/firefox, equal sign
  kb[107] = kb[61] = kb[187] = function (evt) {
    this.change_isolevel_by(evt.shiftKey ? 1 : 0, 0.1);
  };
  // subtract, minus/firefox, dash
  kb[109] = kb[173] = kb[189] = function (evt) {
    this.change_isolevel_by(evt.shiftKey ? 1 : 0, -0.1);
  };
  // [
  kb[219] = function () { this.change_map_radius(-2); };
  // ]
  kb[221] = function () { this.change_map_radius(2); };
  // \ (backslash)
  kb[220] = function (evt) {
    this.select_next('bond lines', 'line_style', LINE_STYLES, evt.shiftKey);
    this.redraw_models();
  };
  // 3, numpad 3
  kb[51] = kb[99] = function () { this.shift_clip(1); };
  // numpad period (Linux), decimal point (Mac)
  kb[108] = kb[110] = function () { this.shift_clip(-1); };
  // shift, ctrl, alt, altgr
  kb[16] = kb[17] = kb[18] = kb[225] = function () {};
  // slash, single quote
  kb[191] = kb[222] = false;  // -> preventDefault()

  this.key_bindings = kb;
};

Viewer.prototype.set_real_space_key_bindings = function () {
  let kb = this.key_bindings;
  // Space
  kb[32] = function (evt) { this.center_next_residue(evt.shiftKey); };
  // p
  kb[80] = function (evt) {
    evt.shiftKey ? this.permalink() : this.go_to_nearest_Ca();
  };
  // t
  kb[84] = function (evt) {
    this.select_next('rendering as', 'render_style', RENDER_STYLES,
                     evt.shiftKey);
    this.redraw_models();
  };
  // y
  kb[89] = function (evt) {
    this.config.hydrogens = !this.config.hydrogens;
    this.hud((this.config.hydrogens ? 'show' : 'hide') +
             ' hydrogens (if any)');
    this.redraw_models();
  };
};

Viewer.prototype.mousedown = function (event) {
  //event.preventDefault(); // default involves setting focus, which we need
  event.stopPropagation();
  document.addEventListener('mouseup', this.mouseup);
  document.addEventListener('mousemove', this.mousemove);
  let state = STATE.NONE;
  if (event.button === 1 || (event.button === 0 && event.ctrlKey)) {
    state = STATE.PAN;
  } else if (event.button === 0) {
    // in Coot shift+Left is labeling atoms like dblclick, + rotation
    if (event.shiftKey) {
      this.dblclick(event);
    }
    state = STATE.ROTATE;
  } else if (event.button === 2) {
    if (event.ctrlKey) {
      state = event.shiftKey ? STATE.ROLL : STATE.SLAB;
    } else {
      this.decor.zoom_grid.visible = true;
      state = STATE.ZOOM;
    }
  }
  this.controls.start(state, this.relX(event), this.relY(event));
  this.request_render();
};

Viewer.prototype.dblclick = function (event) {
  if (event.button !== 0) return;
  if (this.decor.selection) {
    this.remove_and_dispose(this.decor.selection);
    this.decor.selection = null;
  }
  const mouse = new THREE.Vector2(this.relX(event), this.relY(event));
  const atom = this.pick_atom(mouse, this.camera);
  if (atom) {
    this.hud(atom.long_label());
    this.toggle_label(atom);
    const color = this.config.colors[atom.element] || this.config.colors.def;
    const size = 2.5 * scale_by_height(this.config.bond_line, this.window_size);
    this.decor.selection = makeWheels([atom], [color], size);
    this.scene.add(this.decor.selection);
  } else {
    this.hud();
  }
  this.request_render();
};

// for two-finger touch events
function touch_info(evt) {
  const touches = evt.touches;
  const dx = touches[0].pageX - touches[1].pageX;
  const dy = touches[0].pageY - touches[1].pageY;
  return {pageX: (touches[0].pageX + touches[1].pageX) / 2,
          pageY: (touches[0].pageY + touches[1].pageY) / 2,
          dist: Math.sqrt(dx * dx + dy * dy)};
}

Viewer.prototype.touchstart = function (event) {
  const touches = event.touches;
  if (touches.length === 1) {
    this.controls.start(STATE.ROTATE,
                        this.relX(touches[0]), this.relY(touches[0]));
  } else { // for now using only two touches
    const info = touch_info(event);
    this.controls.start(STATE.PAN_ZOOM,
                        this.relX(info), this.relY(info), info.dist);
  }
  this.request_render();
};

Viewer.prototype.touchmove = function (event) {
  event.preventDefault();
  event.stopPropagation();
  const touches = event.touches;
  if (touches.length === 1) {
    this.controls.move(this.relX(touches[0]), this.relY(touches[0]));
  } else { // for now using only two touches
    const info = touch_info(event);
    this.controls.move(this.relX(info), this.relY(info), info.dist);
  }
};

Viewer.prototype.touchend = function (/*event*/) {
  this.controls.stop();
  this.redraw_maps();
};

Viewer.prototype.mousewheel = function (evt) {
  evt.preventDefault();
  evt.stopPropagation();
  // evt.wheelDelta for WebKit, evt.detail for Firefox
  const delta = evt.wheelDelta || -2 * (evt.detail || 0);
  this.mousewheel_action(delta, evt);
  this.request_render();
};

Viewer.prototype.mousewheel_action = function (delta, evt) {
  this.change_isolevel_by(evt.shiftKey ? 1 : 0, 0.0005 * delta);
};

Viewer.prototype.resize = function (/*evt*/) {
  const el = this.container;
  const width = el.clientWidth;
  const height = el.clientHeight;
  this.window_offset[0] = el.offsetLeft;
  this.window_offset[1] = el.offsetTop;
  this.camera.left = -width;
  this.camera.right = width;
  this.camera.top = height;
  this.camera.bottom = -height;
  this.camera.updateProjectionMatrix();
  this.renderer.setSize(width, height);
  if (width !== this.window_size[0] || height !== this.window_size[1]) {
    this.window_size[0] = width;
    this.window_size[1] = height;
    this.redraw_models(); // b/c bond_line is scaled by height
  }
  this.request_render();
};

// makes sense only for full-window viewer
function parse_url_fragment() {
  let ret = {};
  if (typeof window === 'undefined') return ret;
  const params = window.location.hash.substr(1).split('&');
  for (let i = 0; i < params.length; i++) {
    const kv = params[i].split('=');
    let val = kv[1];
    if (kv[0] === 'xyz' || kv[0] === 'eye') {
      val = val.split(',').map(Number);
    } else if (kv[0] === 'zoom') {
      val = Number(val);
    }
    ret[kv[0]] = val;
  }
  return ret;
}

// If xyz set recenter on it looking toward the model center.
// Otherwise recenter on the model center looking along the z axis.
Viewer.prototype.recenter = function (xyz, cam, steps) {
  const bag = this.active_model_bag;
  let new_up;
  if (xyz != null && cam == null && bag !== null) {
    // look from specified point toward the center of the molecule,
    // i.e. shift camera away from the molecule center.
    xyz = new THREE.Vector3(xyz[0], xyz[1], xyz[2]);
    const mc = bag.model.get_center();
    let d = new THREE.Vector3(xyz[0] - mc[0], xyz[1] - mc[1], xyz[2] - mc[2]);
    d.setLength(100);
    new_up = d.y < 90 ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0);
    new_up.projectOnPlane(d).normalize();
    cam = d.add(xyz);
  } else {
    xyz = xyz || (bag ? bag.model.get_center() : [0, 0, 0]);
    if (cam != null) {
      cam = new THREE.Vector3(cam[0], cam[1], cam[2]);
      new_up = null; // preserve the up direction
    } else {
      const dc = this.default_camera_pos;
      cam = new THREE.Vector3(xyz[0] + dc[0], xyz[1] + dc[1], xyz[2] + dc[2]);
      new_up = THREE.Object3D.DefaultUp; // Vector3(0, 1, 0)
    }
  }
  this.controls.go_to(xyz, cam, new_up, steps);
};

Viewer.prototype.center_next_residue = function (back) {
  if (!this.active_model_bag) return;
  const a = this.active_model_bag.model.next_residue(this.selected_atom, back);
  if (a) this.select_atom(a);
};

Viewer.prototype.select_atom = function (atom, options) {
  options = options || {};
  this.hud('-> ' + atom.long_label());
  const steps = options.steps || 30. / auto_speed;
  this.controls.go_to(atom.xyz, null, null, steps);
  this.toggle_label(this.selected_atom);
  this.selected_atom = atom;
  this.toggle_label(atom);
};

Viewer.prototype.update_camera = function () {
  const dxyz = this.camera.position.distanceTo(this.target);
  const w = this.controls.slab_width;
  const scale = w.length === 3 ? w[2] : this.camera.zoom;
  this.camera.near = dxyz * (1 - w[0] / scale);
  this.camera.far = dxyz * (1 + w[1] / scale);
  //this.light.position.copy(this.camera.position);
  this.camera.updateProjectionMatrix();
};

// The main loop. Running when a mouse button is pressed or when the view
// is moving (and run once more after the mouse button is released).
// It is also triggered by keydown events.
Viewer.prototype.render = function () {
  this.scheduled = true;
  if (this.renderer === null) return;
  if (this.controls.update()) {
    this.update_camera();
  }
  const tied = this.tied_viewer;
  if (!this.controls.is_going()) {
    this.redraw_maps();
    if (tied && !tied.scheduled) tied.redraw_maps();
  }
  this.renderer.render(this.scene, this.camera);
  if (tied && !tied.scheduled) tied.renderer.render(tied.scene, tied.camera);
  if (this.nav) {
    this.nav.renderer.render(this.nav.scene, this.camera);
  }
  this.scheduled = false;
  if (this.controls.is_moving()) {
    this.request_render();
  }
  if (this.stats) {
    this.stats.update();
  }
};

Viewer.prototype.request_render = function () {
  if (typeof window !== 'undefined' && !this.scheduled) {
    this.scheduled = true;
    window.requestAnimationFrame(this.render.bind(this));
  }
};

Viewer.prototype.set_model = function (model) {
  const model_bag = new ModelBag(model, this.config, this.window_size);
  this.model_bags.push(model_bag);
  this.set_atomic_objects(model_bag);
  this.active_model_bag = model_bag;
  this.request_render();
};

Viewer.prototype.add_map = function (map, is_diff_map) {
  //map.show_debug_info();
  const map_bag = new MapBag(map, is_diff_map);
  this.map_bags.push(map_bag);
  this.add_el_objects(map_bag);
  this.request_render();
};

Viewer.prototype.load_file = function (url/*:string*/,
                                       options/*:{[id:string]: mixed}*/,
                                       callback/*:Function*/) {
  if (this.renderer === null) return;  // no WebGL detected
  let req = new XMLHttpRequest();
  req.open('GET', url, true);
  if (options.binary) {
    req.responseType = 'arraybuffer';
  } else {
    // http://stackoverflow.com/questions/7374911/
    req.overrideMimeType('text/plain');
  }
  let self = this;
  req.onreadystatechange = function () {
    if (req.readyState === 4) {
      // chrome --allow-file-access-from-files gives status 0
      if (req.status === 200 || (req.status === 0 && req.response !== null)) {
        try {
          callback(req);
        } catch (e) {
          self.hud('Error: ' + e.message + '\nin ' + url, 'ERR');
        }
      } else {
        self.hud('Failed to fetch ' + url, 'ERR');
      }
    }
  };
  if (options.progress) {
    // flow-ignore-line  dom.js in flow is incomplete
    req.addEventListener('progress', function (evt /*:ProgressEvent*/) {
      if (evt.lengthComputable && evt.loaded && evt.total) {
        const fn = url.split('/').pop();
        self.hud('loading ' + fn + ' ... ' +
                 (evt.loaded >> 10) + ' / ' + (evt.total >> 10) + ' kB');
        if (evt.loaded === evt.total) self.hud(); // clear progress message
      }
    });
  }
  try {
    req.send(null);
  } catch (e) {
    self.hud('loading ' + url + ' failed:\n' + e, 'ERR');
  }
};

Viewer.prototype.set_view = function (options) {
  const frag = parse_url_fragment();
  if (frag.zoom) this.camera.zoom = frag.zoom;
  this.recenter(frag.xyz || (options && options.center), frag.eye, 1);
};

// Load molecular model from PDB file and centers the view
Viewer.prototype.load_pdb = function (url, options, callback) {
  let self = this;
  this.load_file(url, {binary: false}, function (req) {
    let model = new Model();
    model.from_pdb(req.responseText);
    self.set_model(model);
    self.set_view(options);
    if (callback) callback();
  });
};

Viewer.prototype.load_map = function (url, options, callback) {
  if (options.format !== 'ccp4' && options.format !== 'dsn6') {
    throw Error('Unknown map format.');
  }
  let self = this;
  this.load_file(url, {binary: true, progress: true}, function (req) {
    let map = new ElMap();
    if (options.format === 'ccp4') map.from_ccp4(req.response, true);
    else /* === 'dsn6'*/ map.from_dsn6(req.response);
    self.add_map(map, options.diff_map);
    if (callback) callback();
  });
};

// Load a normal map and a difference map.
// To show the first map ASAP we do not download both maps in parallel.
Viewer.prototype.load_ccp4_maps = function (url1, url2, callback) {
  let self = this;
  this.load_map(url1, {diff_map: false, format: 'ccp4'}, function () {
    self.load_map(url2, {diff_map: true, format: 'ccp4'}, callback);
  });
};

// Load a model (PDB), normal map and a difference map - in this order.
Viewer.prototype.load_pdb_and_ccp4_maps = function (pdb, map1, map2, callback) {
  let self = this;
  this.load_pdb(pdb, {}, function () {
    self.load_ccp4_maps(map1, map2, callback);
  });
};

// TODO: navigation window like in gimp and mifit
/*
Viewer.prototype.show_nav = function (inset_id) {
  var inset = document.getElementById(inset_id);
  if (!inset) return;
  inset.style.display = 'block';
  var nav = {};
  nav.renderer = new THREE.WebGLRenderer();
  nav.renderer.setClearColor(0x555555, 1);
  nav.renderer.setSize(200, 200);
  inset.appendChild(nav.renderer.domElement);
  //nav.scene = new THREE.Scene();
  nav.scene = this.scene;
  //var light = new THREE.AmbientLight(0xffffff);
  //nav.scene.add(light);
  this.nav = nav;
};
*/

Viewer.prototype.ColorSchemes = ColorSchemes;

