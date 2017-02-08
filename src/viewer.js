// @flow

import * as THREE from 'three';
import { makeLineMaterial, makeLineSegments, makeLine, makeRibbon,
         makeChickenWire, makeGrid, makeWheels, makeCube,
         makeRgbBox, makeLabel, addXyzCross } from './lines.js';
import { STATE, Controls } from './controls.js';
import { ElMap } from './elmap.js';
import { modelsFromPDB } from './model.js';

/*::
 import type {AtomT, Model} from './model.js'

 type ColorScheme = {
   name: string,
   bg: number,
   fg: number,
   [name:string]: number | number[],
 };
 type Num3 = [number, number, number];
 */

const ColorSchemes /*:ColorScheme[]*/ = [ // Viewer.prototype.ColorSchemes
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


// options handled by select_next()

const COLOR_AIMS = ['element', 'B-factor', 'occupancy', 'index', 'chain'];
const RENDER_STYLES = ['lines', 'trace', 'ribbon'/*, 'ball&stick'*/];
const MAP_STYLES = ['marching cubes', 'squarish'/*, 'snapped MC'*/];
const LINE_STYLES = ['normal', 'simplistic'];
const LABEL_FONTS = ['bold 14px', '14px', '16px', 'bold 16px'];

function rainbow_value(v/*:number*/, vmin/*:number*/, vmax/*:number*/) {
  let c = new THREE.Color(0xe0e0e0);
  if (vmin < vmax) {
    const ratio = (v - vmin) / (vmax - vmin);
    const hue = (240 - (240 * ratio)) / 360;
    c.setHSL(hue, 1.0, 0.5);
  }
  return c;
}

function color_by(style, atoms /*:AtomT[]*/, elem_colors, hue_shift) {
  let color_func;
  const last_atom = atoms[atoms.length-1];
  if (style === 'index') {
    color_func = function (atom) {
      return rainbow_value(atom.i_seq, 0, last_atom.i_seq);
    };
  } else if (style === 'B-factor') {
    let vmin = Infinity;
    let vmax = -Infinity;
    for (let i = 0; i < atoms.length; i++) {
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
    if (hue_shift === 0) {
      color_func = function (atom) {
        return elem_colors[atom.element] || elem_colors.def;
      };
    } else {
      const c_col = elem_colors['C'].clone().offsetHSL(hue_shift, 0, 0);
      color_func = function (atom) {
        const el = atom.element;
        return el === 'C' ? c_col : (elem_colors[el] || elem_colors.def);
      };
    }
  }
  return atoms.map(color_func);
}

function scale_by_height(value, size) { // for scaling bond_line
  return value * size[1] / 700;
}

class MapBag {
  /*::
  map: ElMap
  name: string
  isolevel: number
  visible: boolean
  types: string[]
  block_ctr: THREE.Vector3
  el_objects: Object[]
  */
  constructor(map, config, is_diff_map) {
    this.map = map;
    this.name = '';
    this.isolevel = is_diff_map ? 3.0 : config.default_isolevel;
    this.visible = true;
    this.types = is_diff_map ? ['map_pos', 'map_neg'] : ['map_den'];
    this.block_ctr = new THREE.Vector3(Infinity, 0, 0);
    this.el_objects = []; // three.js objects
  }
}


class ModelBag {
  /*::
  model: Model
  label: string
  visible: boolean
  hue_shift: number
  conf: Object
  win_size: [number, number]
  atomic_objects: Object[]
  static ctor_counter: number
  */
  constructor(model, config, win_size) {
    this.model = model;
    this.label = '(model #' + ++ModelBag.ctor_counter + ')';
    this.visible = true;
    this.hue_shift = 0;
    this.conf = config;
    this.win_size = win_size;
    this.atomic_objects = []; // list of three.js objects
  }

  get_visible_atoms() {
    const atoms = this.model.atoms;
    if (this.conf.hydrogens || !this.model.has_hydrogens) {
      return atoms;
    }
    // with filter() it's twice slower (on Node 4.2)
    //return atoms.filter(function(a) { return a.element !== 'H'; });
    const non_h = [];
    for (const atom of atoms) {
      if (atom.element !== 'H') non_h.push(atom);
    }
    return non_h;
  }

  add_bonds(ligands_only, ball_size) {
    const visible_atoms = this.get_visible_atoms();
    const color_style = ligands_only ? 'element' : this.conf.color_aim;
    const colors = color_by(color_style, visible_atoms,
                            this.conf.colors, this.hue_shift);
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
          const vatom = new THREE.Vector3(atom.xyz[0], atom.xyz[1],
                                          atom.xyz[2]);
          if (opt.balls && ball_size != null) {
            const lerp_factor = vatom.distanceTo(vmid) / ball_size;
            vatom.lerp(vmid, lerp_factor);
          }
          vertex_arr.push(vatom, vmid);
          color_arr.push(color, color);
        }
      }
    }
    if (vertex_arr.length === 0) return;
    const linewidth = scale_by_height(this.conf.bond_line, this.win_size);
    const use_gl_lines = this.conf.line_style === 'simplistic';
    const material = makeLineMaterial({
      gl_lines: use_gl_lines,
      linewidth: linewidth,
      win_size: this.win_size,
      segments: true,
    });
    this.atomic_objects.push(makeLineSegments(material, vertex_arr, color_arr));
    if (opt.balls && ball_size != null) {
      this.atomic_objects.push(makeWheels(visible_atoms, colors, ball_size));
    } else if (!use_gl_lines && !ligands_only) {
      // wheels (discs) as simplistic round caps
      this.atomic_objects.push(makeWheels(visible_atoms, colors, linewidth));
    }
  }

  add_trace() {
    const segments = this.model.extract_trace();
    const visible_atoms = [].concat.apply([], segments);
    const colors = color_by(this.conf.color_aim, visible_atoms,
                            this.conf.colors, this.hue_shift);
    const material = makeLineMaterial({
      gl_lines: this.conf.line_style === 'simplistic',
      linewidth: scale_by_height(this.conf.bond_line, this.win_size),
      win_size: this.win_size,
    });
    let k = 0;
    for (const seg of segments) {
      const color_slice = colors.slice(k, k + seg.length);
      k += seg.length;
      const line = makeLine(material, seg, color_slice);
      this.atomic_objects.push(line);
    }
  }

  add_ribbon(smoothness) {
    const segments = this.model.extract_trace();
    const res_map = this.model.get_residues();
    const visible_atoms = [].concat.apply([], segments);
    const colors = color_by(this.conf.color_aim, visible_atoms,
                            this.conf.colors, this.hue_shift);
    let k = 0;
    for (const seg of segments) {
      let tangents = [];
      let last = [0, 0, 0];
      for (const atom of seg) {
        const residue = res_map[atom.resid()];
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
  }
}

ModelBag.ctor_counter = 0;

function vec3_to_fixed(vec, n) {
  return [vec.x.toFixed(n), vec.y.toFixed(n), vec.z.toFixed(n)];
}

// for two-finger touch events
function touch_info(evt/*:TouchEvent*/) {
  const touches = evt.touches;
  const dx = touches[0].pageX - touches[1].pageX;
  const dy = touches[0].pageY - touches[1].pageY;
  return {pageX: (touches[0].pageX + touches[1].pageX) / 2,
          pageY: (touches[0].pageY + touches[1].pageY) / 2,
          dist: Math.sqrt(dx * dx + dy * dy)};
}

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


export class Viewer {
  /*::
  model_bags: ModelBag[]
  map_bags: MapBag[]
  decor: {cell_box: ?Object , selection: ?Object, zoom_grid: Object,
          mark: ?Object}
  labels: {[id:string]: {o: THREE.Mesh, bag: ModelBag}}
  nav: ?Object
  config: Object
  window_size: [number, number]
  window_offset: [number, number]
  last_ctr: THREE.Vector3
  selected: {bag: ?ModelBag, atom: ?AtomT}
  scene: THREE.Scene
  light: THREE.Light
  default_camera_pos: [number, number, number]
  target: THREE.Vector3
  camera: THREE.OrthographicCamera
  controls: Controls
  tied_viewer: ?Viewer
  raycaster: THREE.Raycaster
  renderer: THREE.WebGLRenderer
  container: ?HTMLElement
  hud_el: ?HTMLElement
  help_el: ?HTMLElement
  initial_hud_html: string
  scheduled: boolean
  ColorSchemes: ColorScheme[]
  MOUSE_HELP: string
  KEYBOARD_HELP: string
  ABOUT_HELP: string
  stats: ?Object
  mousemove: (MouseEvent) => void
  mouseup: (MouseEvent) => void
  key_bindings: Array<Function|false>
  */
  constructor(options /*: {[key: string]: any}*/) {
    // rendered objects
    this.model_bags = [];
    this.map_bags = [];
    this.decor = { cell_box: null, selection: null, zoom_grid: makeGrid(),
                   mark: null };
    this.labels = {};
    this.nav = null;

    this.config = {
      bond_line: 4.0, // ~ to height, like in Coot (see scale_by_height())
      map_line: 1.25,  // for any height
      map_radius: 10.0,
      max_map_radius: 40,
      default_isolevel: 1.5,
      center_cube_size: 0.1,
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
    this.selected = {bag: null, atom: null};
    this.scene = new THREE.Scene();
    this.scene.fog = new THREE.Fog(this.config.colors.bg, 0, 1);
    this.light = new THREE.AmbientLight(0xffffff);
    this.scene.add(this.light);
    this.default_camera_pos = [0, 0, 100];
    if (options.share_view) {
      this.target = options.share_view.target;
      this.camera = options.share_view.camera;
      this.controls = options.share_view.controls;
      this.tied_viewer = options.share_view;
      this.tied_viewer.tied_viewer = this; // not GC friendly
    } else {
      this.target = new THREE.Vector3(0, 0, 0);
      this.camera = new THREE.OrthographicCamera();
      this.camera.position.fromArray(this.default_camera_pos);
      this.controls = new Controls(this.camera, this.target);
    }
    this.raycaster = new THREE.Raycaster();
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
    }

    if (this.container == null) return; // can be null in headless tests
    this.renderer.setClearColor(this.config.colors.bg, 1);
    this.renderer.setPixelRatio(window.devicePixelRatio);
    this.resize();
    this.camera.zoom = this.camera.right / 35.0;  // arbitrary choice
    this.update_camera();
    const el = this.renderer.domElement;
    // $FlowFixMe: flow can't figure out that this.container != null
    this.container.appendChild(el);
    if (options.focusable) {
      el.tabIndex = 0;
    }
    this.decor.zoom_grid.visible = false;
    this.scene.add(this.decor.zoom_grid);
    if (window.Stats) { // set by including three/examples/js/libs/stats.min.js
      this.stats = new window.Stats();
      // $FlowFixMe: flow can't figure out that this.container != null
      this.container.appendChild(this.stats.dom);
    }

    window.addEventListener('resize', this.resize.bind(this));
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

    this.mousemove = function (event/*:MouseEvent*/) {
      event.preventDefault();
      //event.stopPropagation();
      self.controls.move(self.relX(event), self.relY(event));
    };

    this.mouseup = function (event/*:MouseEvent*/) {
      event.preventDefault();
      event.stopPropagation();
      document.removeEventListener('mousemove', self.mousemove);
      document.removeEventListener('mouseup', self.mouseup);
      self.decor.zoom_grid.visible = false;
      const not_panned = self.controls.stop();
      // special case - centering on atoms after action 'pan' with no shift
      if (not_panned) {
        const pick = self.pick_atom(not_panned, self.camera);
        if (pick != null) {
          self.select_atom(pick, {steps: 60});
        }
      }
      self.redraw_maps();
    };

    this.scheduled = false;
    this.request_render();
  }

  pick_atom(coords/*:THREE.Vector2*/, camera/*:THREE.OrthographicCamera*/) {
    for (const bag of this.model_bags) {
      if (!bag.visible) continue;
      this.raycaster.setFromCamera(coords, camera);
      this.raycaster.near = camera.near;
      // '0.15' b/c the furthest 15% is hardly visible in the fog
      this.raycaster.far = camera.far - 0.15 * (camera.far - camera.near);
      this.raycaster.linePrecision = 0.3;
      let intersects = this.raycaster.intersectObjects(bag.atomic_objects);
      if (intersects.length > 0) {
        intersects.sort(function (x) { return x.line_dist || Infinity; });
        const p = intersects[0].point;
        const atom = bag.model.get_nearest_atom(p.x, p.y, p.z);
        if (atom != null) {
          return {bag, atom};
        }
      }
    }
  }

  set_colors(scheme/*:?number|string|ColorScheme*/) {
    function to_col(x) { return new THREE.Color(x); }
    if (scheme == null) {
      scheme = this.config.colors;
    } else if (typeof scheme === 'number') {
      scheme = this.ColorSchemes[scheme % this.ColorSchemes.length];
    } else if (typeof scheme === 'string') {
      for (const sc of this.ColorSchemes) {
        if (sc.name === scheme) {
          scheme = sc;
          break;
        }
      }
      throw Error('Unknown color scheme.');
    }
    if (typeof scheme.bg === 'number') {
      for (const key in scheme) {
        if (key !== 'name') {
          scheme[key] = scheme[key] instanceof Array ? scheme[key].map(to_col)
                                                     : to_col(scheme[key]);
        }
      }
    }
    this.decor.zoom_grid.color_value.set(scheme.fg);
    this.config.config = scheme;
    this.redraw_all();
  }

  // relative position on canvas in normalized device coordinates [-1, +1]
  relX(evt/*:{pageX: number}*/) {
    return 2 * (evt.pageX - this.window_offset[0]) / this.window_size[0] - 1;
  }

  relY(evt/*:{pageY: number}*/) {
    return 1 - 2 * (evt.pageY - this.window_offset[1]) / this.window_size[1];
  }

  hud(text/*:?string*/, type/*:?string*/) {
    if (typeof document === 'undefined') return;  // for testing on node
    let el = this.hud_el;
    if (el) {
      if (text != null) {
        if (type === 'HTML') {
          el.innerHTML = text;
        } else {
          el.textContent = text;
        }
      } else {
        el.innerHTML = this.initial_hud_html;
      }
      const err = (type === 'ERR');
      el.style.backgroundColor = (err ? '#b00' : '');
      if (err && text) console.log('ERR: ' + text);
    } else {
      console.log('hud:', text);
    }
  }

  redraw_center(force/*:?boolean*/) {
    const size = this.config.center_cube_size;
    if (force ||
        this.target.distanceToSquared(this.last_ctr) > 0.01 * size * size) {
      this.last_ctr.copy(this.target);
      if (this.decor.mark) {
        this.scene.remove(this.decor.mark);
      }
      this.decor.mark = makeCube(size, this.target, {
        color: this.config.colors.center,
        linewidth: 2,
        win_size: this.window_size,
      });
      this.scene.add(this.decor.mark);
    }
  }

  redraw_maps(force/*:?boolean*/) {
    this.redraw_center(force);
    const r = this.config.map_radius;
    for (const map_bag of this.map_bags) {
      if (force || this.target.distanceToSquared(map_bag.block_ctr) > r/100) {
        this.redraw_map(map_bag);
      }
    }
  }

  remove_and_dispose(obj/*:Object*/) {
    this.scene.remove(obj);
    if (obj.geometry) obj.geometry.dispose();
    if (obj.material) {
      if (obj.material.uniforms && obj.material.uniforms.map) {
        obj.material.uniforms.map.value.dispose();
      }
      obj.material.dispose();
    }
    for (let o of obj.children) {
      this.remove_and_dispose(o);
    }
  }

  clear_el_objects(map_bag/*:MapBag*/) {
    for (let o of map_bag.el_objects) {
      this.remove_and_dispose(o);
    }
    map_bag.el_objects = [];
  }

  clear_atomic_objects(model_bag/*:ModelBag*/) {
    for (let o of model_bag.atomic_objects) {
      this.remove_and_dispose(o);
    }
    model_bag.atomic_objects = [];
  }

  set_atomic_objects(model_bag/*:ModelBag*/) {
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
    for (let o of model_bag.atomic_objects) {
      this.scene.add(o);
    }
  }

  // Add/remove label if `show` is specified, toggle otherwise.
  toggle_label(pick/*:{bag:?ModelBag, atom:?AtomT}*/, show/*:?boolean*/) {
    if (pick.atom == null) return;
    const text = pick.atom.short_label();
    const uid = text; // we assume that the labels inside one model are unique
    const is_shown = (uid in this.labels);
    if (show === undefined) show = !is_shown;
    if (show) {
      if (is_shown) return;
      if (pick.atom == null) return; // silly flow
      const label = makeLabel(text, {
        pos: pick.atom.xyz,
        font: this.config.label_font,
        color: '#' + this.config.colors.fg.getHexString(),
        win_size: this.window_size,
      });
      if (!label) return;
      if (pick.bag == null) return;
      this.labels[uid] = { o: label, bag: pick.bag };
      this.scene.add(label);
    } else {
      if (!is_shown) return;
      this.remove_and_dispose(this.labels[uid].o);
      delete this.labels[uid];
    }
  }

  redraw_labels() {
    for (let uid in this.labels) { // eslint-disable-line guard-for-in
      const text = uid;
      this.labels[uid].o.remake(text, {
        font: this.config.label_font,
        color: '#' + this.config.colors.fg.getHexString(),
      });
    }
  }

  toggle_map_visibility(map_bag/*:MapBag*/) {
    if (typeof map_bag === 'number') {
      map_bag = this.map_bags[map_bag];
    }
    map_bag.visible = !map_bag.visible;
    this.redraw_map(map_bag);
    this.request_render();
  }

  redraw_map(map_bag/*:MapBag*/) {
    this.clear_el_objects(map_bag);
    if (map_bag.visible) {
      map_bag.map.block.clear();
      this.add_el_objects(map_bag);
    }
  }

  toggle_model_visibility(model_bag/*:?ModelBag*/) {
    model_bag = model_bag || this.selected.bag;
    if (model_bag == null) return;
    model_bag.visible = !model_bag.visible;
    this.redraw_model(model_bag);
    this.request_render();
  }

  redraw_model(model_bag/*:ModelBag*/) {
    this.clear_atomic_objects(model_bag);
    if (model_bag.visible) {
      this.set_atomic_objects(model_bag);
    }
  }

  redraw_models() {
    for (const model_bag of this.model_bags) {
      this.redraw_model(model_bag);
    }
  }

  add_el_objects(map_bag/*:MapBag*/) {
    if (!map_bag.visible || this.config.map_radius <= 0) return;
    if (map_bag.map.block.empty()) {
      const t = this.target;
      map_bag.block_ctr.copy(t);
      map_bag.map.extract_block(this.config.map_radius, [t.x, t.y, t.z]);
    }
    for (const mtype of map_bag.types) {
      const isolevel = (mtype === 'map_neg' ? -1 : 1) * map_bag.isolevel;
      const iso = map_bag.map.isomesh_in_block(isolevel, this.config.map_style);
      if (iso == null) continue;
      const obj = makeChickenWire(iso, {
        color: this.config.colors[mtype],
        linewidth: this.config.map_line,
      });
      map_bag.el_objects.push(obj);
      this.scene.add(obj);
    }
  }

  change_isolevel_by(map_idx/*:number*/, delta/*:number*/) {
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
    this.hud('map ' + (map_idx+1) + ' level =  ' + abs_text + ' ' +
             map_bag.map.unit + ' (' + map_bag.isolevel.toFixed(2) + ' rmsd)');
  }

  change_map_radius(delta/*:number*/) {
    const rmax = this.config.max_map_radius;
    const cf = this.config;
    cf.map_radius = Math.min(Math.max(cf.map_radius + delta, 0), rmax);
    cf.map_radius = Math.round(cf.map_radius * 1e9) / 1e9;
    let info = 'map "radius": ' + cf.map_radius;
    if (cf.map_radius === rmax) info += ' (max)';
    else if (cf.map_radius === 0) info += ' (hidden maps)';
    if (this.map_bags.length === 0) info += '\nNB: no map is loaded.';
    this.hud(info);
    this.redraw_maps(true);
  }

  change_slab_width_by(delta/*:number*/) {
    let slab_width = this.controls.slab_width;
    slab_width[0] = Math.max(slab_width[0] + delta, 0.01);
    slab_width[1] = Math.max(slab_width[1] + delta, 0.01);
    this.update_camera();
    const final_width = this.camera.far - this.camera.near;
    this.hud('clip width: ' + final_width.toPrecision(3));
  }

  change_zoom_by_factor(mult/*:number*/) {
    this.camera.zoom *= mult;
    this.update_camera();
    this.hud('zoom: ' + this.camera.zoom.toPrecision(3));
  }

  change_bond_line(delta/*:number*/) {
    this.config.bond_line = Math.max(this.config.bond_line + delta, 0.1);
    this.redraw_models();
    this.hud('bond width: ' + scale_by_height(this.config.bond_line,
                                              this.window_size).toFixed(1));
  }

  change_map_line(delta/*:number*/) {
    this.config.map_line = Math.max(this.config.map_line + delta, 0.1);
    this.redraw_maps(true);
    this.hud('wireframe width: ' + this.config.map_line.toFixed(1));
  }

  toggle_full_screen() {
    let d = document;
    if (d.fullscreenElement || d.mozFullScreenElement ||
        d.webkitFullscreenElement || d.msFullscreenElement) {
      let ex = d.exitFullscreen || d.webkitExitFullscreen ||
      // $FlowFixMe: property `msExitFullscreen` not found in document
               d.mozCancelFullScreen || d.msExitFullscreen;
      // $FlowFixMe: cannot call property `exitFullscreen` of unknown type
      if (ex) ex.call(d);
    } else {
      let el = this.container;
      if (!el) return;
      let req = el.requestFullscreen || el.webkitRequestFullscreen ||
      // $FlowFixMe: property `msRequestFullscreen` not found in HTMLElement
                el.mozRequestFullScreen || el.msRequestFullscreen;
      // $FlowFixMe
      if (req) req.call(el);
    }
  }

  toggle_cell_box() {
    if (this.decor.cell_box) {
      this.scene.remove(this.decor.cell_box);
      this.decor.cell_box = null;
    } else {
      const uc_func = this.get_cell_box_func();
      if (uc_func) {
        this.decor.cell_box = makeRgbBox(uc_func, {
          color: this.config.colors.fg,
        });
        this.scene.add(this.decor.cell_box);
      }
    }
  }

  get_cell_box_func() /*:?Function*/ {
    let uc = null;
    if (this.selected.bag != null) {
      uc = this.selected.bag.model.unit_cell;
    }
    // note: model may not have unit cell
    if (uc == null && this.map_bags.length > 0) {
      uc = this.map_bags[0].map.unit_cell;
    }
    return uc && uc.orthogonalize.bind(uc);
  }

  shift_clip(delta/*:number*/) {
    let eye = this.camera.position.clone().sub(this.target);
    eye.multiplyScalar(delta / eye.length());
    this.target.add(eye);
    this.camera.position.add(eye);
    this.update_camera();
    this.redraw_maps();
    this.hud('clip shifted by [' + vec3_to_fixed(eye, 2).join(' ') + ']');
  }

  go_to_nearest_Ca() {
    const t = this.target;
    const bag = this.selected.bag;
    if (bag == null) return;
    const atom = bag.model.get_nearest_atom(t.x, t.y, t.z, 'CA');
    if (atom != null) {
      this.select_atom({bag, atom}, {steps: 30});
    } else {
      this.hud('no nearby CA');
    }
  }

  permalink() {
    if (typeof window === 'undefined') return;
    const xyz_prec = Math.round(-Math.log10(0.001));
    window.location.hash =
      '#xyz=' + vec3_to_fixed(this.target, xyz_prec).join(',') +
      '&eye=' + vec3_to_fixed(this.camera.position, 1).join(',') +
      '&zoom=' + this.camera.zoom.toFixed(0);
    this.hud('copy URL from the location bar');
  }

  redraw_all() {
    if (!this.renderer) return;
    this.scene.fog.color = this.config.colors.bg;
    if (this.renderer) this.renderer.setClearColor(this.config.colors.bg, 1);
    this.redraw_models();
    this.redraw_maps(true);
    this.redraw_labels();
  }

  toggle_help() {
    let el = this.help_el;
    if (!el) return;
    el.style.display = el.style.display === 'block' ? 'none' : 'block';
    if (el.innerHTML === '') {
      el.innerHTML = [this.MOUSE_HELP, this.KEYBOARD_HELP,
                      this.ABOUT_HELP].join('\n\n');
    }
  }

  select_next(info/*:string*/, key/*:string*/,
              options/*:Array<string|ColorScheme>*/, back/*:boolean*/) {
    const old_idx = options.indexOf(this.config[key]);
    const len = options.length;
    const new_idx = (old_idx + (back ? len - 1 : 1)) % len;
    this.config[key] = options[new_idx];
    let html = info + ':';
    for (let i = 0; i < len; i++) {
      const tag = (i === new_idx ? 'u' : 's');
      const opt_name = typeof options[i] === 'string' ? options[i]
                                                      : options[i].name;
      html += ' <' + tag + '>' + opt_name + '</' + tag + '>';
    }
    this.hud(html, 'HTML');
  }

  keydown(evt/*:KeyboardEvent*/) {
    const action = this.key_bindings[evt.keyCode];
    if (action) {
      (action.bind(this))(evt);
    } else {
      if (action === false) evt.preventDefault();
      if (this.help_el) this.hud('Nothing here. Press H for help.');
    }
    this.request_render();
  }

  set_common_key_bindings() {
    let kb = new Array(256);
    // b
    kb[66] = function (evt) {
      this.select_next('color scheme', 'colors', this.ColorSchemes,
                       evt.shiftKey);
      this.set_colors();
    };
    // c
    kb[67] = function (evt) {
      this.select_next('coloring by', 'color_aim', COLOR_AIMS, evt.shiftKey);
      this.redraw_models();
    };
    // e
    kb[69] = function toggle_fog() {
      const fog = this.scene.fog;
      const has_fog = (fog.far === 1);
      fog.far = (has_fog ? 1e9 : 1);
      this.hud((has_fog ? 'dis': 'en') + 'able fog');
      this.redraw_all();
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
    // shift, ctrl, alt, altgr
    kb[16] = kb[17] = kb[18] = kb[225] = function () {};
    // slash, single quote
    kb[191] = kb[222] = false;  // -> preventDefault()

    this.key_bindings = kb;
  }

  set_real_space_key_bindings() {
    let kb = this.key_bindings;
    // Home
    kb[36] = function (evt) {
      evt.ctrlKey ? this.change_map_line(0.1) : this.change_bond_line(0.2);
    };
    // End
    kb[35] = function (evt) {
      evt.ctrlKey ? this.change_map_line(-0.1) : this.change_bond_line(-0.2);
    };
    // Space
    kb[32] = function (evt) { this.center_next_residue(evt.shiftKey); };
    // d
    kb[68] = function () { this.change_slab_width_by(-0.1); };
    // f
    kb[70] = function (evt) {
      evt.shiftKey ? this.toggle_full_screen() : this.change_slab_width_by(0.1);
    };
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
    // u
    kb[85] = function () {
      this.hud('toggled unit cell box');
      this.toggle_cell_box();
    };
    // y
    kb[89] = function (evt) {
      this.config.hydrogens = !this.config.hydrogens;
      this.hud((this.config.hydrogens ? 'show' : 'hide') +
               ' hydrogens (if any)');
      this.redraw_models();
    };
    // comma
    kb[188] = function (evt) { if (evt.shiftKey) this.shift_clip(1); };
    // period
    kb[190] = function (evt) { if (evt.shiftKey) this.shift_clip(-1); };
  }

  mousedown(event/*:MouseEvent*/) {
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
  }

  dblclick(event/*:MouseEvent*/) {
    if (event.button !== 0) return;
    if (this.decor.selection) {
      this.remove_and_dispose(this.decor.selection);
      this.decor.selection = null;
    }
    const mouse = new THREE.Vector2(this.relX(event), this.relY(event));
    const pick = this.pick_atom(mouse, this.camera);
    if (pick) {
      const atom = pick.atom;
      this.hud(pick.bag.label + ' ' + atom.long_label());
      this.toggle_label(pick);
      const color = this.config.colors[atom.element] || this.config.colors.def;
      const size = 2.5 * scale_by_height(this.config.bond_line,
                                         this.window_size);
      this.decor.selection = makeWheels([atom], [color], size);
      this.scene.add(this.decor.selection);
    } else {
      this.hud();
    }
    this.request_render();
  }

  touchstart(event/*:TouchEvent*/) {
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
  }

  touchmove(event/*:TouchEvent*/) {
    event.preventDefault();
    event.stopPropagation();
    const touches = event.touches;
    if (touches.length === 1) {
      this.controls.move(this.relX(touches[0]), this.relY(touches[0]));
    } else { // for now using only two touches
      const info = touch_info(event);
      this.controls.move(this.relX(info), this.relY(info), info.dist);
    }
  }

  touchend(/*event*/) {
    this.controls.stop();
    this.redraw_maps();
  }

  // $FlowFixMe TODO: wheel()+WheelEvent are more standard
  mousewheel(evt/*:MouseWheelEvent*/) {
    evt.preventDefault();
    evt.stopPropagation();
    // evt.wheelDelta for WebKit, evt.detail for Firefox
    const delta = evt.wheelDelta || -2 * (evt.detail || 0);
    this.mousewheel_action(delta, evt);
    this.request_render();
  }

  mousewheel_action(delta/*:number*/, evt/*:WheelEvent*/) {
    this.change_isolevel_by(evt.shiftKey ? 1 : 0, 0.0005 * delta);
  }

  resize(/*evt*/) {
    const el = this.container;
    if (el == null) return;
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
  }

  // If xyz set recenter on it looking toward the model center.
  // Otherwise recenter on the model center looking along the z axis.
  recenter(xyz/*:?Num3*/, cam/*:?Num3*/, steps/*:number*/) {
    const bag = this.selected.bag;
    let new_up;
    if (xyz != null && cam == null && bag != null) {
      // look from specified point toward the center of the molecule,
      // i.e. shift camera away from the molecule center.
      xyz = new THREE.Vector3(xyz[0], xyz[1], xyz[2]);
      const mc = bag.model.get_center();
      let d = new THREE.Vector3(xyz[0] - mc[0], xyz[1] - mc[1], xyz[2] - mc[2]);
      d.setLength(100);
      new_up = d.y < 90 ? new THREE.Vector3(0, 1, 0)
                        : new THREE.Vector3(1, 0, 0);
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
  }

  center_next_residue(back/*:boolean*/) {
    const bag = this.selected.bag;
    if (bag == null) return;
    const atom = bag.model.next_residue(this.selected.atom, back);
    if (atom != null) {
      this.select_atom({bag, atom}, {steps: 30});
    }
  }

  select_atom(pick/*:{bag:ModelBag, atom:AtomT}*/, options/*:Object*/={}) {
    this.hud('-> ' + pick.bag.label + ' ' + pick.atom.long_label());
    this.controls.go_to(pick.atom.xyz, null, null, options.steps);
    this.toggle_label(this.selected, false);
    this.selected = {bag: pick.bag, atom: pick.atom}; // not ...=pick b/c flow
    this.toggle_label(this.selected, true);
  }

  update_camera() {
    const dxyz = this.camera.position.distanceTo(this.target);
    const w = this.controls.slab_width;
    const scale = w[2] || this.camera.zoom;
    this.camera.near = dxyz * (1 - w[0] / scale);
    this.camera.far = dxyz * (1 + w[1] / scale);
    //this.light.position.copy(this.camera.position);
    this.camera.updateProjectionMatrix();
  }

  // The main loop. Running when a mouse button is pressed or when the view
  // is moving (and run once more after the mouse button is released).
  // It is also triggered by keydown events.
  render() {
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
  }

  request_render() {
    if (typeof window !== 'undefined' && !this.scheduled) {
      this.scheduled = true;
      window.requestAnimationFrame(this.render.bind(this));
    }
  }

  add_model(model/*:Model*/, options/*:Object*/={}) {
    const model_bag = new ModelBag(model, this.config, this.window_size);
    model_bag.hue_shift = options.hue_shift || 0.06 * this.model_bags.length;
    this.model_bags.push(model_bag);
    this.set_atomic_objects(model_bag);
    this.request_render();
  }

  add_map(map/*:ElMap*/, is_diff_map/*:boolean*/) {
    //map.show_debug_info();
    const map_bag = new MapBag(map, this.config, is_diff_map);
    this.map_bags.push(map_bag);
    this.add_el_objects(map_bag);
    this.request_render();
  }

  load_file(url/*:string*/, options/*:{[id:string]: mixed}*/,
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
  }

  set_view(options/*:?Object*/) {
    const frag = parse_url_fragment();
    if (frag.zoom) this.camera.zoom = frag.zoom;
    this.recenter(frag.xyz || (options && options.center), frag.eye, 1);
  }

  // Load molecular model from PDB file and centers the view
  load_pdb(url/*:string*/, options/*:?Object*/, callback/*:?Function*/) {
    let self = this;
    this.load_file(url, {binary: false}, function (req) {
      const len = self.model_bags.length;
      const models = modelsFromPDB(req.responseText);
      for (const model of models) {
        self.add_model(model);
      }
      self.selected.bag = self.model_bags[len];
      self.set_view(options);
      if (callback) callback();
    });
  }

  load_map(url/*:string*/, options/*:Object*/, callback/*:?Function*/) {
    if (options.format !== 'ccp4' && options.format !== 'dsn6') {
      throw Error('Unknown map format.');
    }
    let self = this;
    this.load_file(url, {binary: true, progress: true}, function (req) {
      self.load_map_from_buffer(req.response, options);
      if (callback) callback();
    });
  }

  load_map_from_buffer(buffer/*:ArrayBuffer*/, options/*:Object*/) {
    let map = new ElMap();
    if (options.format === 'dsn6') {
      map.from_dsn6(buffer);
    } else {
      map.from_ccp4(buffer, true);
    }
    this.add_map(map, options.diff_map);
  }

  // Load a normal map and a difference map.
  // To show the first map ASAP we do not download both maps in parallel.
  load_ccp4_maps(url1/*:string*/, url2/*:string*/, callback/*:?Function*/) {
    let self = this;
    this.load_map(url1, {diff_map: false, format: 'ccp4'}, function () {
      self.load_map(url2, {diff_map: true, format: 'ccp4'}, callback);
    });
  }

  // Load a model (PDB), normal map and a difference map - in this order.
  load_pdb_and_ccp4_maps(pdb/*:string*/, map1/*:string*/, map2/*:string*/,
                         callback/*:?Function*/) {
    let self = this;
    this.load_pdb(pdb, {}, function () {
      self.load_ccp4_maps(map1, map2, callback);
    });
  }

  load_from_pdbe() {
    if (typeof window === 'undefined') return;
    const url = window.location.href;
    const match = url.match(/[?&]id=([^&#]+)/);
    if (match == null) return;
    const id = match[1].toLowerCase();
    this.load_pdb('https://www.ebi.ac.uk/pdbe/entry-files/pdb' + id + '.ent');
  }

  // TODO: navigation window like in gimp and mifit
  /*
  show_nav(inset_id) {
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
}

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
  'E = toggle fog',
  'Q = label font',
  '+/- = sigma level',
  ']/[ = map radius',
  'D/F = clip width',
  '&lt;/> = move clip',
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
  '&nbsp; <a href="https://uglymol.github.io">uglymol</a> ' +
  // $FlowFixMe
  (typeof VERSION === 'string' ? VERSION : 'dev'); // eslint-disable-line

Viewer.prototype.ColorSchemes = ColorSchemes;

