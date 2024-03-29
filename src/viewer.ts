import { OrthographicCamera, Scene, AmbientLight, Color, Vector3,
         Ray, WebGLRenderer, Fog } from './uthree/main';
import { makeLineMaterial, makeLineSegments, makeRibbon,
         makeChickenWire, makeGrid, makeSticks, makeBalls, makeWheels, makeCube,
         makeRgbBox, Label, addXyzCross } from './draw';
import { STATE, Controls } from './controls';
import { ElMap } from './elmap';
import { modelsFromPDB } from './model';

import type { Atom, Model } from './model';
import type { LineSegments } from './uthree/main';
import type { OrCameraType } from './controls';

type Num2 = [number, number];
type Num3 = [number, number, number];

type ColorScheme = {
  bg: Color,
  fg: Color,
  map_den?: Color,
  map_pos?: Color,
  map_neg?: Color,
  center?: Color,
  H?: Color,
  C?: Color,
  N?: Color,
  O?: Color,
  S?: Color,
  P?: Color,
  MG?: Color,
  CL?: Color,
  CA?: Color,
  MN?: Color,
  FE?: Color,
  NI?: Color,
  def?: Color,
};

export type ViewerConfig = {
  bond_line: number,
  map_line: number,
  map_radius: number,
  max_map_radius: number,
  default_isolevel: number,
  center_cube_size: number,
  map_style: string,
  render_style: string,
  ligand_style: string,
  water_style: string,
  color_prop: string,
  line_style: string,
  label_font: string,
  color_scheme: string,
  colors?: ColorScheme,
  hydrogens: boolean,
  ball_size: number,
  stay?: boolean;
};

const ColorSchemes: Record<string, ColorScheme> = {
  // the default scheme that generally mimicks Coot
  'coot dark': {
    bg: new Color(0x000000),
    fg: new Color(0xFFFFFF),
    map_den: new Color(0x3362B2),
    map_pos: new Color(0x298029),
    map_neg: new Color(0x8B2E2E),
    center: new Color(0xC997B0),
    // atoms
    H: new Color(0x858585), // H is normally invisible
    // C, N and O are taken approximately (by color-picker) from coot
    C: new Color(0xb3b300),
    N: new Color(0x7EAAFB),
    O: new Color(0xF24984),
    S: new Color(0x40ff40), // S in coot is too similar to C, here it's greener
    // Coot doesn't define other colors (?)
    MG: new Color(0xc0c0c0),
    P:  new Color(0xffc040),
    CL: new Color(0xa0ff60),
    CA: new Color(0xffffff),
    MN: new Color(0xff90c0),
    FE: new Color(0xa03000),
    NI: new Color(0x00ff80),
    def: new Color(0xa0a0a0), // default atom color
  },

  // scheme made of "solarized" colors (http://ethanschoonover.com/solarized):
  // base03  base02  base01  base00  base0   base1   base2   base3
  // #002b36 #073642 #586e75 #657b83 #839496 #93a1a1 #eee8d5 #fdf6e3
  // yellow  orange  red     magenta violet  blue    cyan    green
  // #b58900 #cb4b16 #dc322f #d33682 #6c71c4 #268bd2 #2aa198 #859900
  'solarized dark': {
    bg: new Color(0x002b36),
    fg: new Color(0xfdf6e3),
    map_den: new Color(0x268bd2),
    map_pos: new Color(0x859900),
    map_neg: new Color(0xd33682),
    center: new Color(0xfdf6e3),
    H: new Color(0x586e75),
    C: new Color(0x93a1a1),
    N: new Color(0x6c71c4),
    O: new Color(0xcb4b16),
    S: new Color(0xb58900),
    def: new Color(0xeee8d5),
  },

  'solarized light': {
    bg: new Color(0xfdf6e3),
    fg: new Color(0x002b36),
    map_den: new Color(0x268bd2),
    map_pos: new Color(0x859900),
    map_neg: new Color(0xd33682),
    center: new Color(0x002b36),
    H: new Color(0x93a1a1),
    C: new Color(0x586e75),
    N: new Color(0x6c71c4),
    O: new Color(0xcb4b16),
    S: new Color(0xb58900),
    def: new Color(0x073642),
  },

  // like in Coot after Edit > Background Color > White
  'coot light': {
    bg: new Color(0xFFFFFF),
    fg: new Color(0x000000),
    map_den: new Color(0x3362B2),
    map_pos: new Color(0x298029),
    map_neg: new Color(0x8B2E2E),
    center: new Color(0xC7C769),
    H: new Color(0x999999),
    C: new Color(0xA96464),
    N: new Color(0x1C51B3),
    O: new Color(0xC33869),
    S: new Color(0x9E7B3D),
    def: new Color(0x808080),
  },
};


const INIT_HUD_TEXT = 'This is UglyMol not Coot. ' +
  '<a href="#" onclick="V.toggle_help(); return false;">H shows help.</a>';

// options handled by select_next()

const COLOR_PROPS = ['element', 'B-factor', 'pLDDT', 'occupancy', 'index', 'chain'];
const RENDER_STYLES = ['lines', 'trace', 'ribbon', 'ball&stick'];
const LIGAND_STYLES = ['ball&stick', 'lines'];
const WATER_STYLES = ['cross', 'dot', 'invisible'];
const MAP_STYLES = ['marching cubes', 'squarish'/*, 'snapped MC'*/];
const LINE_STYLES = ['normal', 'simplistic'];
const LABEL_FONTS = ['bold 14px', '14px', '16px', 'bold 16px'];

function rainbow_value(v: number, vmin: number, vmax: number) {
  const c = new Color(0xe0e0e0);
  if (vmin < vmax) {
    const ratio = (v - vmin) / (vmax - vmin);
    const hue = (240 - (240 * ratio)) / 360;
    c.setHSL(hue, 1.0, 0.5);
  }
  return c;
}

function color_by(prop: string, atoms: Atom[], elem_colors: ColorScheme,
                  hue_shift: number): Color[] {
  let color_func;
  const last_atom = atoms[atoms.length-1];
  if (prop === 'index') {
    color_func = function (atom: Atom) {
      return rainbow_value(atom.i_seq, 0, last_atom.i_seq);
    };
  } else if (prop === 'B-factor') {
    let vmin = Infinity;
    let vmax = -Infinity;
    for (let i = 0; i < atoms.length; i++) {
      const v = atoms[i].b;
      if (v > vmax) vmax = v;
      if (v < vmin) vmin = v;
    }
    //console.log('B-factors in [' + vmin + ', ' + vmax + ']');
    color_func = function (atom: Atom) {
      return rainbow_value(atom.b, vmin, vmax);
    };
  } else if (prop === 'pLDDT') {
    const steps = [90, 70, 50];
    const colors = [
      new Color(0x0053d6), // dark blue
      new Color(0x65cbf3), // light blue
      new Color(0xffdb13), // yellow
      new Color(0xff7d45)  // orange
    ];
    color_func = function (atom: Atom) {
      let i = 0;
      while (i < 3 && atom.b < steps[i]) {
        ++i;
      }
      return colors[i];
    };
  } else if (prop === 'occupancy') {
    color_func = function (atom: Atom) {
      return rainbow_value(atom.occ, 0, 1);
    };
  } else if (prop === 'chain') {
    color_func = function (atom: Atom) {
      return rainbow_value(atom.chain_index, 0, last_atom.chain_index);
    };
  } else { // element
    if (hue_shift === 0) {
      color_func = function (atom: Atom) {
        return elem_colors[atom.element] || elem_colors.def;
      };
    } else {
      const c_hsl = { h: 0, s: 0, l: 0 };
      elem_colors['C'].getHSL(c_hsl);
      const c_col = new Color(0, 0, 0);
      c_col.setHSL(c_hsl.h + hue_shift, c_hsl.s, c_hsl.l);
      color_func = function (atom: Atom) {
        const el = atom.element;
        return el === 'C' ? c_col : (elem_colors[el] || elem_colors.def);
      };
    }
  }
  return atoms.map(color_func);
}

function scale_by_height(value: number, size: Num2) { // for scaling bond_line
  return value * size[1] / 700;
}

class MapBag {
  map: ElMap;
  name: string;
  isolevel: number;
  visible: boolean;
  types: string[];
  block_ctr: Vector3;
  el_objects: object[];

  constructor(map: ElMap, config: ViewerConfig, is_diff_map: boolean) {
    this.map = map;
    this.name = '';
    this.isolevel = is_diff_map ? 3.0 : config.default_isolevel;
    this.visible = true;
    this.types = is_diff_map ? ['map_pos', 'map_neg'] : ['map_den'];
    this.block_ctr = new Vector3(Infinity, 0, 0);
    this.el_objects = []; // three.js objects
  }
}

class ModelBag {
  model: Model;
  label: string;
  visible: boolean;
  hue_shift: number;
  conf: ViewerConfig;
  win_size: Num2;
  objects: object[];
  atom_array: Atom[]
  static ctor_counter: number;

  constructor(model: Model, config: ViewerConfig, win_size: Num2) {
    this.model = model;
    this.label = '(model #' + ++ModelBag.ctor_counter + ')';
    this.visible = true;
    this.hue_shift = 0;
    this.conf = config;
    this.win_size = win_size;
    this.objects = []; // list of three.js objects
    this.atom_array = [];
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

  add_bonds(polymers: boolean, ligands: boolean, ball_size?: number) {
    const visible_atoms = this.get_visible_atoms();
    const colors = color_by(this.conf.color_prop, visible_atoms,
                            this.conf.colors, this.hue_shift);
    const vertex_arr: Num3[] = [];
    const color_arr = [];
    const sphere_arr = [];
    const sphere_color_arr = [];
    const hydrogens = this.conf.hydrogens;
    for (let i = 0; i < visible_atoms.length; i++) {
      const atom = visible_atoms[i];
      const color = colors[i];
      if (!(atom.is_ligand ? ligands : polymers)) continue;
      if (atom.is_water() && this.conf.water_style === 'invisible') continue;
      if (atom.bonds.length === 0 && ball_size == null) { // nonbonded - cross
        if (!atom.is_water() || this.conf.water_style === 'cross') {
          addXyzCross(vertex_arr, atom.xyz, 0.7);
          for (let n = 0; n < 6; n++) {
            color_arr.push(color);
          }
        }
      } else { // bonded, draw lines
        for (let j = 0; j < atom.bonds.length; j++) {
          const other = this.model.atoms[atom.bonds[j]];
          if (!hydrogens && other.element === 'H') continue;
          // Coot show X-H bonds as thinner lines in a single color.
          // Here we keep it simple and render such bonds like all others.
          const mid = atom.midpoint(other);
          vertex_arr.push(atom.xyz, mid);
          color_arr.push(color, color);
        }
      }
      sphere_arr.push(atom);
      sphere_color_arr.push(color);
    }

    if (ball_size != null) {
      if (vertex_arr.length !== 0) {
        this.objects.push(makeSticks(vertex_arr, color_arr, ball_size / 2));
      }
      if (sphere_arr.length !== 0) {
        this.objects.push(makeBalls(sphere_arr, sphere_color_arr, ball_size));
      }
    } else if (vertex_arr.length !== 0) {
      const linewidth = scale_by_height(this.conf.bond_line, this.win_size);
      const material = makeLineMaterial({
        linewidth: linewidth,
        win_size: this.win_size,
      });
      this.objects.push(makeLineSegments(material, vertex_arr, color_arr));
      if (this.conf.line_style !== 'simplistic') {
        // wheels (discs) as round caps
        this.objects.push(makeWheels(sphere_arr, sphere_color_arr, linewidth));
      }
    }

    sphere_arr.forEach(function (v) { this.atom_array.push(v); }, this);
  }

  add_trace() {
    const segments = this.model.extract_trace();
    const visible_atoms = [].concat.apply([], segments);
    const colors = color_by(this.conf.color_prop, visible_atoms,
                            this.conf.colors, this.hue_shift);
    const vertex_arr: Num3[] = [];
    const color_arr = [];
    let k = 0;
    for (const seg of segments) {
      for (let i = 1; i < seg.length; ++i) {
        vertex_arr.push(seg[i-1].xyz, seg[i].xyz);
        color_arr.push(colors[k+i-1], colors[k+i]);
      }
      k += seg.length;
    }
    const linewidth = scale_by_height(this.conf.bond_line, this.win_size);
    const material = makeLineMaterial({
      linewidth: linewidth,
      win_size: this.win_size,
    });
    this.objects.push(makeLineSegments(material, vertex_arr, color_arr));
    if (this.conf.line_style !== 'simplistic') {
      // wheels (discs) as round caps
      this.objects.push(makeWheels(visible_atoms, colors, linewidth));
    }
    this.atom_array = visible_atoms;
  }

  add_ribbon(smoothness: number) {
    const segments = this.model.extract_trace();
    const res_map = this.model.get_residues();
    const visible_atoms = [].concat.apply([], segments);
    const colors = color_by(this.conf.color_prop, visible_atoms,
                            this.conf.colors, this.hue_shift);
    let k = 0;
    for (const seg of segments) {
      const tangents = [];
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
      this.objects.push(obj);
    }
  }
}

ModelBag.ctor_counter = 0;

function vec3_to_fixed(vec, n) {
  return [vec.x.toFixed(n), vec.y.toFixed(n), vec.z.toFixed(n)];
}

// for two-finger touch events
function touch_info(evt: TouchEvent) {
  const touches = evt.touches;
  const dx = touches[0].pageX - touches[1].pageX;
  const dy = touches[0].pageY - touches[1].pageY;
  return {pageX: (touches[0].pageX + touches[1].pageX) / 2,
          pageY: (touches[0].pageY + touches[1].pageY) / 2,
          dist: Math.sqrt(dx * dx + dy * dy)};
}

// makes sense only for full-window viewer
function parse_url_fragment() {
  const ret : Record<string, any> = {};
  if (typeof window === 'undefined') return ret;
  const params = window.location.hash.substr(1).split('&');
  for (let i = 0; i < params.length; i++) {
    const kv = params[i].split('=');
    const key = kv[0];
    const val = kv[1];
    if (key === 'xyz' || key === 'eye') {
      ret[key] = val.split(',').map(Number);
    } else if (key === 'zoom') {
      ret[key] = Number(val);
    } else {
      ret[key] = val;
    }
  }
  return ret;
}


export class Viewer {
  model_bags: ModelBag[];
  map_bags: MapBag[];
  decor: {
      cell_box: object | null,
      selection: object | null,
      zoom_grid: LineSegments,
      mark: object | null
  };
  labels: {[index:string]: {o: Label, bag: ModelBag}};
  //nav: object | null;
  xhr_headers: Record<string, string>;
  config: ViewerConfig;
  window_size: Num2;
  window_offset: Num2;
  last_ctr: Vector3;
  selected: {bag: ModelBag | null, atom: Atom | null};
  dbl_click_callback: (arg: object) => void;
  scene: Scene;
  light: AmbientLight;
  default_camera_pos: Num3;
  target: Vector3;
  camera: OrCameraType;
  controls: Controls;
  tied_viewer: Viewer | null;
  renderer: WebGLRenderer;
  container: HTMLElement | null;
  hud_el: HTMLElement | null;
  help_el: HTMLElement | null;
  initial_hud_html: string;
  scheduled: boolean;
  declare MOUSE_HELP: string;
  declare KEYBOARD_HELP: string;
  declare ABOUT_HELP: string;
  mousemove: (arg: MouseEvent) => void;
  mouseup: (arg: MouseEvent) => void;
  key_bindings: Array<((evt: KeyboardEvent) => void) | false | undefined>;
  declare ColorSchemes: typeof ColorSchemes;

  constructor(options: Record<string, any>) {
    // rendered objects
    this.model_bags = [];
    this.map_bags = [];
    this.decor = {
      cell_box: null,
      selection: null,
      zoom_grid: makeGrid(),
      mark: null,
    };
    this.labels = {};
    //this.nav = null;
    this.xhr_headers = {};

    this.config = {
      bond_line: 4.0, // ~ to height, like in Coot (see scale_by_height())
      map_line: 1.25,  // for any height
      map_radius: 10.0,
      max_map_radius: 40,
      default_isolevel: 1.5,
      center_cube_size: 0.1,
      map_style: MAP_STYLES[0],
      render_style: RENDER_STYLES[0],
      ligand_style: LIGAND_STYLES[0],
      water_style: WATER_STYLES[0],
      color_prop: COLOR_PROPS[0],
      line_style: LINE_STYLES[0],
      label_font: LABEL_FONTS[0],
      color_scheme: 'coot dark',
      // `colors` is assigned in set_colors()
      hydrogens: false,
      ball_size: 0.4,
    };

    // options of the constructor overwrite default values of the config
    for (const o of Object.keys(options)) {
      if (o in this.config) {
        this.config[o] = options[o];
      }
    }

    this.set_colors();
    this.window_size = [1, 1]; // it will be set in resize()
    this.window_offset = [0, 0];

    this.last_ctr = new Vector3(Infinity, 0, 0);
    this.selected = {bag: null, atom: null};
    this.dbl_click_callback = this.toggle_label;
    this.scene = new Scene();
    this.scene.fog = new Fog(this.config.colors.bg, 0, 1);
    this.scene.add(new AmbientLight(0xffffff));
    this.default_camera_pos = [0, 0, 100];
    if (options.share_view) {
      this.target = options.share_view.target;
      this.camera = options.share_view.camera;
      this.controls = options.share_view.controls;
      this.tied_viewer = options.share_view;
      this.tied_viewer.tied_viewer = this; // not GC friendly
    } else {
      this.target = new Vector3(0, 0, 0);
      this.camera = new OrthographicCamera() as OrCameraType;
      this.camera.position.fromArray(this.default_camera_pos);
      this.controls = new Controls(this.camera, this.target);
    }
    this.set_common_key_bindings();
    if (this.constructor === Viewer) this.set_real_space_key_bindings();

    function get_elem(name) {
      if (options[name] === null || typeof document === 'undefined') return null;
      return document.getElementById(options[name] || name);
    }
    this.hud_el = get_elem('hud');
    this.container = get_elem('viewer');
    this.help_el = get_elem('help');
    if (this.hud_el) {
      if (this.hud_el.innerHTML === '') this.hud_el.innerHTML = INIT_HUD_TEXT;
      this.initial_hud_html = this.hud_el.innerHTML;
    }

    try {
      this.renderer = new WebGLRenderer({antialias: true});
    } catch (e) {
      this.hud('No WebGL in your browser?', 'ERR');
      this.renderer = null;
      return;
    }

    if (this.container == null) return; // can be null in headless tests
    this.renderer.setClearColor(this.config.colors.bg, 1);
    this.renderer.setPixelRatio(window.devicePixelRatio);
    this.resize();
    this.camera.zoom = this.camera.right / 35.0;  // arbitrary choice
    this.update_camera();
    const el = this.renderer.domElement;
    this.container.appendChild(el);
    if (options.focusable) {
      el.tabIndex = 0;
    }
    this.decor.zoom_grid.visible = false;
    this.scene.add(this.decor.zoom_grid);

    window.addEventListener('resize', this.resize.bind(this));
    const keydown_el = (options.focusable ? el : window);
    keydown_el.addEventListener('keydown', this.keydown.bind(this));
    el.addEventListener('contextmenu', function (e) { e.preventDefault(); });
    el.addEventListener('wheel', this.wheel.bind(this));
    el.addEventListener('mousedown', this.mousedown.bind(this));
    el.addEventListener('touchstart', this.touchstart.bind(this));
    el.addEventListener('touchmove', this.touchmove.bind(this));
    el.addEventListener('touchend', this.touchend.bind(this));
    el.addEventListener('touchcancel', this.touchend.bind(this));
    el.addEventListener('dblclick', this.dblclick.bind(this));

    const self = this;

    this.mousemove = function (event: MouseEvent) {
      event.preventDefault();
      //event.stopPropagation();
      self.controls.move(self.relX(event), self.relY(event));
    };

    this.mouseup = function (event: MouseEvent) {
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

  pick_atom(coords: Num2, camera: OrCameraType) {
    let pick = null;
    for (const bag of this.model_bags) {
      if (!bag.visible) continue;
      const z = (camera.near + camera.far) / (camera.near - camera.far);
      const ray = new Ray();
      ray.origin.set(coords[0], coords[1], z).unproject(camera);
      ray.direction.set(0, 0, -1).transformDirection(camera.matrixWorld);
      const near = camera.near;
      // '0.15' b/c the furthest 15% is hardly visible in the fog
      const far = camera.far - 0.15 * (camera.far - camera.near);
      /*
      // previous version - line-based search
      let intersects = [];
      for (const object of bag.objects) {
        if (object.visible === false) continue;
        if (object.userData.bond_lines) {
          line_raycast(object, {ray, near, far, precision: 0.3}, intersects);
        }
      }
      ...
      if (intersects.length > 0) {
        intersects.sort(function (x) { return x.dist2 || Infinity; });
        const p = intersects[0].point;
        const atom = bag.model.get_nearest_atom(p.x, p.y, p.z);
        if (atom != null) {
          return {bag, atom};
        }
      }
      */
      // search directly atom array ignoring matrixWorld
      const vec = new Vector3();
      // required picking precision: 0.35A at zoom 50, 0.27A @z30, 0.44 @z80
      const precision2 = 0.35 * 0.35 * 0.02 * camera.zoom;
      for (const atom of bag.atom_array) {
        vec.set(atom.xyz[0] - ray.origin.x,
                atom.xyz[1] - ray.origin.y,
                atom.xyz[2] - ray.origin.z);
        const distance = vec.dot(ray.direction);
        if (distance < 0 || distance < near || distance > far) continue;
        const diff2 = vec.addScaledVector(ray.direction, -distance).lengthSq();
        if (diff2 > precision2) continue;
        if (pick == null || distance < pick.distance) {
          pick = {bag, atom, distance};
        }
      }
    }
    return pick;
  }

  set_colors() {
    const scheme = ColorSchemes[this.config.color_scheme];
    if (!scheme) throw Error('Unknown color scheme.');
    this.decor.zoom_grid.material.uniforms.ucolor.value.set(scheme.fg);
    this.config.colors = scheme;
    this.redraw_all();
  }

  // relative position on canvas in normalized device coordinates [-1, +1]
  relX(evt: {pageX: number}) {
    return 2 * (evt.pageX - this.window_offset[0]) / this.window_size[0] - 1;
  }

  relY(evt: {pageY: number}) {
    return 1 - 2 * (evt.pageY - this.window_offset[1]) / this.window_size[1];
  }

  hud(text?: string, type?: string) {
    if (typeof document === 'undefined') return;  // for testing on node
    const el = this.hud_el;
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

  redraw_center(force?: boolean) {
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
      });
      this.scene.add(this.decor.mark);
    }
  }

  redraw_maps(force?: boolean) {
    this.redraw_center(force);
    const r = this.config.map_radius;
    for (const map_bag of this.map_bags) {
      if (force || this.target.distanceToSquared(map_bag.block_ctr) > r/100) {
        this.redraw_map(map_bag);
      }
    }
  }

  remove_and_dispose(obj: any) {
    this.scene.remove(obj);
    if (obj.geometry) obj.geometry.dispose();
    if (obj.material) {
      if (obj.material.uniforms && obj.material.uniforms.map) {
        obj.material.uniforms.map.value.dispose();
      }
      obj.material.dispose();
    }
    for (const o of obj.children) {
      this.remove_and_dispose(o);
    }
  }

  clear_el_objects(map_bag: MapBag) {
    for (const o of map_bag.el_objects) {
      this.remove_and_dispose(o);
    }
    map_bag.el_objects = [];
  }

  clear_model_objects(model_bag: ModelBag) {
    for (const o of model_bag.objects) {
      this.remove_and_dispose(o);
    }
    model_bag.objects = [];
  }

  has_frag_depth() {
    return this.renderer && this.renderer.extensions.get('EXT_frag_depth');
  }

  set_model_objects(model_bag: ModelBag) {
    model_bag.objects = [];
    model_bag.atom_array = [];
    let ligand_balls = null;
    if (model_bag.conf.ligand_style === 'ball&stick' && this.has_frag_depth()) {
      ligand_balls = this.config.ball_size;
    }
    switch (model_bag.conf.render_style) {
      case 'lines':
        if (ligand_balls === null) {
          model_bag.add_bonds(true, true);
        } else {
          model_bag.add_bonds(true, false);
          model_bag.add_bonds(false, true, ligand_balls);
        }
        break;
      case 'ball&stick':
        if (!this.has_frag_depth()) {
          this.hud('Ball-and-stick rendering is not working in this browser' +
                   '\ndue to lack of suppport for EXT_frag_depth', 'ERR');
          return;
        }
        if (ligand_balls === null) {
          model_bag.add_bonds(true, false, this.config.ball_size);
          model_bag.add_bonds(false, true);
        } else {
          model_bag.add_bonds(true, true, this.config.ball_size);
        }
        break;
      case 'trace':
        model_bag.add_trace();
        model_bag.add_bonds(false, true, ligand_balls);
        break;
      case 'ribbon':
        model_bag.add_ribbon(8);
        model_bag.add_bonds(false, true, ligand_balls);
        break;
    }
    for (const o of model_bag.objects) {
      this.scene.add(o);
    }
  }

  // Add/remove label if `show` is specified, toggle otherwise.
  toggle_label(pick: {bag?: ModelBag, atom?: Atom}, show?: boolean) {
    if (pick.atom == null) return;
    const text = pick.atom.short_label();
    const uid = text; // we assume that the labels inside one model are unique
    const is_shown = (uid in this.labels);
    if (show === undefined) show = !is_shown;
    if (show) {
      if (is_shown) return;
      const atom_style = pick.atom.is_ligand ? 'ligand_style' : 'render_style';
      const balls = pick.bag && pick.bag.conf[atom_style] === 'ball&stick';
      const label = new Label(text, {
        pos: pick.atom.xyz,
        font: this.config.label_font,
        color: '#' + this.config.colors.fg.getHexString(),
        win_size: this.window_size,
        z_shift: balls ? this.config.ball_size + 0.1 : 0.2,
      });
      if (pick.bag == null || label.mesh == null) return;
      this.labels[uid] = { o: label, bag: pick.bag };
      this.scene.add(label.mesh);
    } else {
      if (!is_shown) return;
      this.remove_and_dispose(this.labels[uid].o.mesh);
      delete this.labels[uid];
    }
  }

  redraw_labels() {
    for (const uid in this.labels) { // eslint-disable-line guard-for-in
      const text = uid;
      this.labels[uid].o.redraw(text, {
        font: this.config.label_font,
        color: '#' + this.config.colors.fg.getHexString(),
      });
    }
  }

  toggle_map_visibility(map_bag: MapBag) {
    if (typeof map_bag === 'number') {
      map_bag = this.map_bags[map_bag];
    }
    map_bag.visible = !map_bag.visible;
    this.redraw_map(map_bag);
    this.request_render();
  }

  redraw_map(map_bag: MapBag) {
    this.clear_el_objects(map_bag);
    if (map_bag.visible) {
      map_bag.map.block.clear();
      this.add_el_objects(map_bag);
    }
  }

  toggle_model_visibility(model_bag?: ModelBag, visible?: boolean) {
    model_bag = model_bag || this.selected.bag;
    if (model_bag == null) return;
    model_bag.visible = visible == null ? !model_bag.visible : visible;
    this.redraw_model(model_bag);
    this.request_render();
  }

  redraw_model(model_bag: ModelBag) {
    this.clear_model_objects(model_bag);
    if (model_bag.visible) {
      this.set_model_objects(model_bag);
    }
  }

  redraw_models() {
    for (const model_bag of this.model_bags) {
      this.redraw_model(model_bag);
    }
  }

  add_el_objects(map_bag: MapBag) {
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

  change_isolevel_by(map_idx: number, delta: number) {
    if (map_idx >= this.map_bags.length) return;
    const map_bag = this.map_bags[map_idx];
    map_bag.isolevel += delta;
    //TODO: move slow part into update()
    this.clear_el_objects(map_bag);
    this.add_el_objects(map_bag);
    const abs_level = map_bag.map.abs_level(map_bag.isolevel);
    let abs_text = abs_level.toFixed(4);
    const tied = this.tied_viewer;
    if (tied && map_idx < tied.map_bags.length) {
      const tied_bag = tied.map_bags[map_idx];
      // Should we tie by sigma or absolute level? Now it's sigma.
      tied_bag.isolevel = map_bag.isolevel;
      abs_text += ' / ' + tied_bag.map.abs_level(tied_bag.isolevel).toFixed(4);
      tied.clear_el_objects(tied_bag);
      tied.add_el_objects(tied_bag);
    }
    this.hud('map ' + (map_idx+1) + ' level =  ' + abs_text + ' ' +
             map_bag.map.unit + ' (' + map_bag.isolevel.toFixed(2) + ' rmsd)');
  }

  change_map_radius(delta: number) {
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

  change_slab_width_by(delta: number) {
    const slab_width = this.controls.slab_width;
    slab_width[0] = Math.max(slab_width[0] + delta, 0.01);
    slab_width[1] = Math.max(slab_width[1] + delta, 0.01);
    this.update_camera();
    const final_width = this.camera.far - this.camera.near;
    this.hud('clip width: ' + final_width.toPrecision(3));
  }

  change_zoom_by_factor(mult: number) {
    this.camera.zoom *= mult;
    this.update_camera();
    this.hud('zoom: ' + this.camera.zoom.toPrecision(3));
  }

  change_bond_line(delta: number) {
    this.config.bond_line = Math.max(this.config.bond_line + delta, 0.1);
    this.redraw_models();
    this.hud('bond width: ' + scale_by_height(this.config.bond_line,
                                              this.window_size).toFixed(1));
  }

  change_map_line(delta: number) {
    this.config.map_line = Math.max(this.config.map_line + delta, 0.1);
    this.redraw_maps(true);
    this.hud('wireframe width: ' + this.config.map_line.toFixed(1));
  }

  toggle_full_screen() {
    const d = document;
    // @ts-expect-error no mozFullScreenElement
    if (d.fullscreenElement || d.mozFullScreenElement ||
        // @ts-expect-error no msFullscreenElement
        d.webkitFullscreenElement || d.msFullscreenElement) {
      // @ts-expect-error no webkitExitFullscreen
      const ex = d.exitFullscreen || d.webkitExitFullscreen ||
                 // @ts-expect-error no msExitFullscreen
                 d.mozCancelFullScreen || d.msExitFullscreen;
      if (ex) ex.call(d);
    } else {
      const el = this.container;
      if (!el) return;
      // @ts-expect-error no webkitRequestFullscreen
      const req = el.requestFullscreen || el.webkitRequestFullscreen ||
                  // @ts-expect-error no msRequestFullscreen
                  el.mozRequestFullScreen || el.msRequestFullscreen;
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
        this.decor.cell_box = makeRgbBox(uc_func, this.config.colors.fg);
        this.scene.add(this.decor.cell_box);
      }
    }
  }

  get_cell_box_func() {
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

  shift_clip(delta: number) {
    const eye = this.camera.position.clone().sub(this.target);
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

  toggle_inactive_models() {
    const n = this.model_bags.length;
    if (n < 2) {
      this.hud((n == 0 ? 'No' : 'Only one') + ' model is loaded. ' +
               '"V" is for working with multiple models.');
      return;
    }
    const show_all = !this.model_bags.every(function (m) { return m.visible; });
    for (const model_bag of this.model_bags) {
      const show = show_all || model_bag === this.selected.bag;
      this.toggle_model_visibility(model_bag, show);
    }
    this.hud(show_all ? 'All models visible' : 'Inactive models hidden');
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
    const el = this.help_el;
    if (!el) return;
    el.style.display = el.style.display === 'block' ? 'none' : 'block';
    if (el.innerHTML === '') {
      el.innerHTML = [this.MOUSE_HELP, this.KEYBOARD_HELP,
                      this.ABOUT_HELP].join('\n\n');
    }
  }

  select_next(info: string, key: string, options: string[], back: boolean) {
    const old_idx = options.indexOf(this.config[key]);
    const len = options.length;
    const new_idx = (old_idx + (back ? len - 1 : 1)) % len;
    this.config[key] = options[new_idx];
    let html = info + ':';
    for (let i = 0; i < len; i++) {
      const tag = (i === new_idx ? 'u' : 's');
      html += ' <' + tag + '>' + options[i] + '</' + tag + '>';
    }
    this.hud(html, 'HTML');
  }

  keydown(evt: KeyboardEvent) {
    if (evt.ctrlKey) return;
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
    const kb = new Array(256);
    // b
    kb[66] = function (this: Viewer, evt: KeyboardEvent) {
      const schemes = Object.keys(this.ColorSchemes);
      this.select_next('color scheme', 'color_scheme', schemes, evt.shiftKey);
      this.set_colors();
    };
    // c
    kb[67] = function (this: Viewer, evt: KeyboardEvent) {
      this.select_next('coloring by', 'color_prop', COLOR_PROPS, evt.shiftKey);
      this.redraw_models();
    };
    // e
    kb[69] = function (this: Viewer) {
      const fog = this.scene.fog;
      const has_fog = (fog.far === 1);
      fog.far = (has_fog ? 1e9 : 1);
      this.hud((has_fog ? 'dis': 'en') + 'able fog');
      this.redraw_all();
    };
    // h
    kb[72] = this.toggle_help;
    // i
    kb[73] = function (this: Viewer, evt: KeyboardEvent) {
      this.hud('toggled spinning');
      this.controls.toggle_auto(evt.shiftKey);
    };
    // k
    kb[75] = function (this: Viewer) {
      this.hud('toggled rocking');
      this.controls.toggle_auto(0.0);
    };
    // m
    kb[77] = function (this: Viewer, evt: KeyboardEvent) {
      this.change_zoom_by_factor(evt.shiftKey ? 1.2 : 1.03);
    };
    // n
    kb[78] = function (this: Viewer, evt: KeyboardEvent) {
      this.change_zoom_by_factor(1 / (evt.shiftKey ? 1.2 : 1.03));
    };
    // q
    kb[81] = function (this: Viewer, evt: KeyboardEvent) {
      this.select_next('label font', 'label_font', LABEL_FONTS, evt.shiftKey);
      this.redraw_labels();
    };
    // r
    kb[82] = function (this: Viewer, evt: KeyboardEvent) {
      if (evt.shiftKey) {
        this.hud('redraw!');
        this.redraw_all();
      } else {
        this.hud('recentered');
        this.recenter();
      }
    };
    // w
    kb[87] = function (this: Viewer, evt: KeyboardEvent) {
      this.select_next('map style', 'map_style', MAP_STYLES, evt.shiftKey);
      this.redraw_maps(true);
    };
    // add, equals/firefox, equal sign
    kb[107] = kb[61] = kb[187] = function (this: Viewer, evt: KeyboardEvent) {
      this.change_isolevel_by(evt.shiftKey ? 1 : 0, 0.1);
    };
    // subtract, minus/firefox, dash
    kb[109] = kb[173] = kb[189] = function (this: Viewer, evt: KeyboardEvent) {
      this.change_isolevel_by(evt.shiftKey ? 1 : 0, -0.1);
    };
    // [
    kb[219] = function (this: Viewer) { this.change_map_radius(-2); };
    // ]
    kb[221] = function (this: Viewer) { this.change_map_radius(2); };
    // \ (backslash)
    kb[220] = function (this: Viewer, evt: KeyboardEvent) {
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
    const kb = this.key_bindings;
    // Home
    kb[36] = function (this: Viewer, evt: KeyboardEvent) {
      evt.shiftKey ? this.change_map_line(0.1) : this.change_bond_line(0.2);
    };
    // End
    kb[35] = function (this: Viewer, evt: KeyboardEvent) {
      evt.shiftKey ? this.change_map_line(-0.1) : this.change_bond_line(-0.2);
    };
    // Space
    kb[32] = function (this: Viewer, evt: KeyboardEvent) {
      this.center_next_residue(evt.shiftKey);
    };
    // d
    kb[68] = function (this: Viewer) {
      this.change_slab_width_by(-0.1);
    };
    // f
    kb[70] = function (this: Viewer, evt: KeyboardEvent) {
      evt.shiftKey ? this.toggle_full_screen() : this.change_slab_width_by(0.1);
    };
    // l
    kb[76] = function (this: Viewer, evt: KeyboardEvent) {
      this.select_next('ligands as', 'ligand_style', LIGAND_STYLES, evt.shiftKey);
      this.redraw_models();
    };
    // p
    kb[80] = function (this: Viewer, evt: KeyboardEvent) {
      evt.shiftKey ? this.permalink() : this.go_to_nearest_Ca();
    };
    // s
    kb[83] = function (this: Viewer, evt: KeyboardEvent) {
      this.select_next('rendering as', 'render_style', RENDER_STYLES, evt.shiftKey);
      this.redraw_models();
    };
    // t
    kb[84] = function (this: Viewer, evt: KeyboardEvent) {
      this.select_next('waters as', 'water_style', WATER_STYLES, evt.shiftKey);
      this.redraw_models();
    };
    // u
    kb[85] = function (this: Viewer) {
      this.hud('toggled unit cell box');
      this.toggle_cell_box();
    };
    // v
    kb[86] = function (this: Viewer) {
      this.toggle_inactive_models();
    };
    // y
    kb[89] = function (this: Viewer) {
      this.config.hydrogens = !this.config.hydrogens;
      this.hud((this.config.hydrogens ? 'show' : 'hide') +
               ' hydrogens (if any)');
      this.redraw_models();
    };
    // comma
    kb[188] = function (this: Viewer, evt: KeyboardEvent) {
      if (evt.shiftKey) this.shift_clip(1);
    };
    // period
    kb[190] = function (this: Viewer, evt: KeyboardEvent) {
      if (evt.shiftKey) this.shift_clip(-1);
    };
  }

  mousedown(event: MouseEvent) {
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

  dblclick(event: MouseEvent) {
    if (event.button !== 0) return;
    if (this.decor.selection) {
      this.remove_and_dispose(this.decor.selection);
      this.decor.selection = null;
    }
    const mouse: Num2 = [this.relX(event), this.relY(event)];
    const pick = this.pick_atom(mouse, this.camera);
    if (pick) {
      const atom = pick.atom;
      this.hud(pick.bag.label + ' ' + atom.long_label());
      this.dbl_click_callback(pick);
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

  touchstart(event: TouchEvent) {
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

  touchmove(event: TouchEvent) {
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

  wheel(evt: WheelEvent) {
    evt.preventDefault();
    evt.stopPropagation();
    this.mousewheel_action(evt.deltaY, evt);
    this.request_render();
  }

  // overrided in ReciprocalViewer
  mousewheel_action(delta: number, evt: WheelEvent) {
    const map_idx = evt.shiftKey ? 1 : 0;
    this.change_isolevel_by(map_idx, 0.0005 * delta);
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
  recenter(xyz?: Num3, cam?: Num3, steps?: number) {
    const bag = this.selected.bag;
    const new_up = new Vector3(0, 1, 0);
    let vec_cam;
    let vec_xyz;
    let eye;
    if (xyz != null && cam == null && bag != null) {
      // look from specified point toward the center of the molecule,
      // i.e. shift camera away from the molecule center.
      const mc = bag.model.get_center();
      eye = new Vector3(xyz[0] - mc[0], xyz[1] - mc[1], xyz[2] - mc[2]);
      eye.setLength(100);
      vec_xyz = new Vector3(xyz[0], xyz[1], xyz[2]);
      vec_cam = eye.clone().add(vec_xyz);
    } else {
      if (xyz == null) {
        if (bag != null) {
          xyz = bag.model.get_center();
        } else {
          const uc_func = this.get_cell_box_func();
          xyz = uc_func ? uc_func([0.5, 0.5, 0.5]) : [0, 0, 0];
        }
      }
      vec_xyz = new Vector3(xyz[0], xyz[1], xyz[2]);
      if (cam != null) {
        vec_cam = new Vector3(cam[0], cam[1], cam[2]);
        eye = vec_cam.clone().sub(vec_xyz);
        new_up.copy(this.camera.up); // preserve the up direction
      } else {
        const dc = this.default_camera_pos;
        vec_cam = new Vector3(xyz[0] + dc[0], xyz[1] + dc[1], xyz[2] + dc[2]);
      }
    }
    if (eye != null) {
      new_up.projectOnPlane(eye);
      if (new_up.lengthSq() < 0.0001) new_up.x += 1;
      new_up.normalize();
    }
    this.controls.go_to(vec_xyz, vec_cam, new_up, steps);
  }

  center_next_residue(back: boolean) {
    const bag = this.selected.bag;
    if (bag == null) return;
    const atom = bag.model.next_residue(this.selected.atom, back);
    if (atom != null) {
      this.select_atom({bag, atom}, {steps: 30});
    }
  }

  select_atom(pick: {bag: ModelBag, atom: Atom}, options: {steps?: number}={}) {
    this.hud('-> ' + pick.bag.label + ' ' + pick.atom.long_label());
    const xyz = pick.atom.xyz;
    this.controls.go_to(new Vector3(xyz[0], xyz[1], xyz[2]),
                        null, null, options.steps);
    this.toggle_label(this.selected, false);
    this.selected = pick;
    this.toggle_label(this.selected, true);
  }

  update_camera() {
    const dxyz = this.camera.position.distanceTo(this.target);
    const w = this.controls.slab_width;
    const scale = w[2] || this.camera.zoom;
    this.camera.near = dxyz * (1 - w[0] / scale);
    this.camera.far = dxyz * (1 + w[1] / scale);
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
    //if (this.nav) {
    //  this.nav.renderer.render(this.nav.scene, this.camera);
    //}
    this.scheduled = false;
    if (this.controls.is_moving()) {
      this.request_render();
    }
  }

  request_render() {
    if (typeof window !== 'undefined' && !this.scheduled) {
      this.scheduled = true;
      window.requestAnimationFrame(this.render.bind(this));
    }
  }

  add_model(model: Model, options: {hue_shift?: number}={}) {
    const model_bag = new ModelBag(model, this.config, this.window_size);
    model_bag.hue_shift = options.hue_shift || 0.06 * this.model_bags.length;
    this.model_bags.push(model_bag);
    this.set_model_objects(model_bag);
    this.request_render();
  }

  add_map(map: ElMap, is_diff_map: boolean) {
    //map.show_debug_info();
    const map_bag = new MapBag(map, this.config, is_diff_map);
    this.map_bags.push(map_bag);
    this.add_el_objects(map_bag);
    this.request_render();
  }

  load_file(url: string, options: Record<string, any>,
            callback: (arg: XMLHttpRequest) => void ) {
    if (this.renderer === null) return;  // no WebGL detected
    const req = new XMLHttpRequest();
    req.open('GET', url, true);
    if (options.binary) {
      req.responseType = 'arraybuffer';
    } else {
      // http://stackoverflow.com/questions/7374911/
      req.overrideMimeType('text/plain');
    }
    const self = this;
    Object.keys(this.xhr_headers).forEach(function (name) {
      req.setRequestHeader(name, self.xhr_headers[name]);
    });
    req.onreadystatechange = function () {
      if (req.readyState === 4) {
        // chrome --allow-file-access-from-files gives status 0
        if (req.status === 200 || (req.status === 0 && req.response !== null &&
                                                       req.response !== '')) {
          try {
            callback(req);
          } catch (e) {
            self.hud('Error: ' + e.message + '\nwhen processing ' + url, 'ERR');
          }
        } else {
          self.hud('Failed to fetch ' + url, 'ERR');
        }
      }
    };
    if (options.progress) {
      req.addEventListener('progress', function (evt: ProgressEvent) {
        if (evt.lengthComputable && evt.loaded && evt.total) {
          const fn = url.split('/').pop();
          self.hud('loading ' + fn + ' ... ' + ((evt.loaded / 1024) | 0) +
                   ' / ' + ((evt.total / 1024) | 0) + ' kB');
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

  set_dropzone(zone: HTMLElement, callback: (arg: File) => void) {
    const self = this;
    zone.addEventListener('dragover', function (e: DragEvent) {
      e.stopPropagation();
      e.preventDefault();
      if (e.dataTransfer != null) e.dataTransfer.dropEffect = 'copy';
      self.hud('ready for file drop...');
    });
    zone.addEventListener('drop', function (e: DragEvent) {
      e.stopPropagation();
      e.preventDefault();
      if (e.dataTransfer == null) return;
      const names = [];
      for (let i = 0; i < e.dataTransfer.files.length; i++) {
        const file = e.dataTransfer.files.item(i);
        try {
          self.hud('Loading ' + file.name);
          callback(file);
        } catch (e) {
          self.hud('Loading ' + file.name + ' failed.\n' + e.message, 'ERR');
          return;
        }
        names.push(file.name);
      }
      self.hud('loaded ' + names.join(', '));
    });
  }

  // for use with set_dropzone
  pick_pdb_and_map(file: File) {
    const self = this;
    const reader = new FileReader();
    if (/\.(pdb|ent)$/.test(file.name)) {
      reader.onload = function (evt) {
        self.load_pdb_from_text(evt.target.result as string);
        self.recenter();
      };
      reader.readAsText(file);
    } else if (/\.(map|ccp4|mrc|dsn6|omap)$/.test(file.name)) {
      const map_format = /\.(dsn6|omap)$/.test(file.name) ? 'dsn6' : 'ccp4';
      reader.onloadend = function (evt) {
        if (evt.target != null && evt.target.readyState == 2) {
          self.load_map_from_buffer(evt.target.result as ArrayBuffer,
                                    {format: map_format});
          if (self.model_bags.length === 0 && self.map_bags.length === 1) {
            self.recenter();
          }
        }
      };
      reader.readAsArrayBuffer(file);
    } else {
      throw Error('Unknown file extension. ' +
                  'Use: pdb, ent, ccp4, mrc, map, dsn6 or omap.');
    }
  }

  set_view(options?: Record<string, any>) {
    const frag = parse_url_fragment();
    if (frag.zoom) this.camera.zoom = frag.zoom;
    this.recenter(frag.xyz || (options && options.center), frag.eye, 1);
  }

  // Load molecular model from PDB file and centers the view
  load_pdb_from_text(text: string) {
    const len = this.model_bags.length;
    const models = modelsFromPDB(text);
    for (const model of models) {
      this.add_model(model);
    }
    this.selected.bag = this.model_bags[len];
  }

  load_pdb(url: string, options?: Record<string, any>,
           callback?: () => void) {
    const self = this;
    this.load_file(url, {binary: false, progress: true}, function (req) {
      self.load_pdb_from_text(req.responseText);
      if (options == null || !options.stay) self.set_view(options);
      if (callback) callback();
    });
  }

  load_map(url: string | null, options: Record<string, any>,
           callback?: () => void) {
    if (url == null) {
      if (callback) callback();
      return;
    }
    if (options.format !== 'ccp4' && options.format !== 'dsn6') {
      throw Error('Unknown map format.');
    }
    const self = this;
    this.load_file(url, {binary: true, progress: true}, function (req) {
      self.load_map_from_buffer(req.response, options);
      if (callback) callback();
    });
  }

  load_map_from_buffer(buffer: ArrayBuffer, options: Record<string, any>) {
    const map = new ElMap();
    if (options.format === 'dsn6') {
      map.from_dsn6(buffer);
    } else {
      map.from_ccp4(buffer, true);
    }
    this.add_map(map, options.diff_map);
  }

  // Load a normal map and a difference map.
  // To show the first map ASAP we do not download both maps in parallel.
  load_maps(url1: string, url2: string,
            options: Record<string, any>, callback?: () => void) {
    const format = options.format || 'ccp4';
    const self = this;
    this.load_map(url1, {diff_map: false, format: format}, function () {
      self.load_map(url2, {diff_map: true, format: format}, callback);
    });
  }

  // Load a model (PDB), normal map and a difference map - in this order.
  load_pdb_and_maps(pdb: string, map1: string, map2: string,
                    options: Record<string, any>, callback?: () => void) {
    const self = this;
    this.load_pdb(pdb, options, function () {
      self.load_maps(map1, map2, options, callback);
    });
  }

  // for backward compatibility:
  load_ccp4_maps(url1: string, url2: string, callback?: () => void) {
    this.load_maps(url1, url2, {format: 'ccp4'}, callback);
  }
  load_pdb_and_ccp4_maps(pdb: string, map1: string, map2: string,
                         callback?: () => void) {
    this.load_pdb_and_maps(pdb, map1, map2, {format: 'ccp4'}, callback);
  }

  // pdb_id here should be lowercase ('1abc')
  load_from_pdbe(pdb_id: string, callback?: () => void) {
    const id = pdb_id.toLowerCase();
    this.load_pdb_and_maps(
      'https://www.ebi.ac.uk/pdbe/entry-files/pdb' + id + '.ent',
      'https://www.ebi.ac.uk/pdbe/coordinates/files/' + id + '.ccp4',
      'https://www.ebi.ac.uk/pdbe/coordinates/files/' + id + '_diff.ccp4',
      {format: 'ccp4'}, callback);
  }
  load_from_rcsb(pdb_id: string, callback?: () => void) {
    const id = pdb_id.toLowerCase();
    this.load_pdb_and_maps(
      'https://files.rcsb.org/download/' + id + '.pdb',
      'https://edmaps.rcsb.org/maps/' + id + '_2fofc.dsn6',
      'https://edmaps.rcsb.org/maps/' + id + '_fofc.dsn6',
      {format: 'dsn6'}, callback);
  }

  // TODO: navigation window like in gimp and mifit
  /*
  show_nav(inset_id) {
    var inset = document.getElementById(inset_id);
    if (!inset) return;
    inset.style.display = 'block';
    var nav = {};
    nav.renderer = new WebGLRenderer();
    nav.renderer.setClearColor(0x555555, 1);
    nav.renderer.setSize(200, 200);
    inset.appendChild(nav.renderer.domElement);
    //nav.scene = new Scene();
    nav.scene = this.scene;
    //nav.scene.add(new AmbientLight(0xffffff));
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
  'S = general style',
  'L = ligand style',
  'T = water style',
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
  'V = inactive models',
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
  // @ts-expect-error Cannot find name 'VERSION'
  (typeof VERSION === 'string' ? VERSION : 'dev'); // eslint-disable-line

Viewer.prototype.ColorSchemes = ColorSchemes;

