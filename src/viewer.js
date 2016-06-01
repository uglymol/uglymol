
var THREE = THREE || require('three'); // eslint-disable-line
var ElMap = ElMap || require('./elmap'); // eslint-disable-line
var Model = Model || require('./model'); // eslint-disable-line
var isosurface = isosurface || require('./isosurface'); // eslint-disable-line

// colors are made global for easier tweaking in a browser
var ColorSchemes = {
  dark: {
    bg: 0x000000,
    map_den: 0x3362B2,
    map_pos: 0x298029,
    map_neg: 0x8B2E2E,
    center: 0xC997B0,
    cell_box: 0xFFFFFF,
    // atoms
    H: 0xf0f0f0, // H is normally invisible
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
    def: 0xa0a0a0 // default atom color
  },
  light: { // like in Coot after Edit > Background Color > White
    bg: 0xFFFFFF,
    map_den: 0x3362B2,
    map_pos: 0x298029,
    map_neg: 0x8B2E2E,
    center: 0xC7C769,
    cell_box: 0x000000,
    H: 0x404040,
    C: 0xA96464,
    N: 0x1C51B3,
    O: 0xC33869,
    S: 0x9E7B3D,
    def: 0x808080
  }
};

var auto_speed = 1.0;

var Viewer = (function () {
'use strict';

function set_colors(palette, o) {
  var scheme = ColorSchemes[palette];
  for (var key in scheme) {
    if (o[key]) {
      o[key].set(scheme[key]);
    } else {
      o[key] = new THREE.Color(scheme[key]);
    }
  }
  o.name = palette;
  return o;
}

var Colors = set_colors('dark', {});

var initial_hud_text = null;

function hud(text) {
  var el = document.getElementById('hud');
  if (el) {
    if (initial_hud_text === null) {
      initial_hud_text = el.textContent;
    }
    el.textContent = (text !== undefined ? text : initial_hud_text);
  } else {
    console.log('hud: ' + text);
  }
}

var target = new THREE.Vector3();
var raycaster = new THREE.Raycaster();
var pickable_model = null;

var selected_atom = null;

function raycast_intersect(coords, camera) {
  if (pickable_model === null) { return null; }
  raycaster.setFromCamera(coords, camera);
  // https://github.com/mrdoob/three.js/issues/9009 and 0.8 because fog
  raycaster.near = 1e-3;
  raycaster.far = 0.8 * (camera.far - camera.near);
  raycaster.linePrecision = 0.2;
  return raycaster.intersectObjects(pickable_model.atomic_objects);
}

function pick_atom(coords, camera) {
  var intersects = raycast_intersect(coords, camera);
  if (intersects.length < 1) { return null; }
  var p = intersects[0].point;
  return pickable_model.model.get_nearest_atom(p.x, p.y, p.z);
}

// relative position on canvas in normalized device coordinates [-1, +1]
function relX(evt) { return 2 * evt.pageX / window.innerWidth - 1; }
function relY(evt) { return 1 - 2 * evt.pageY / window.innerHeight; }

// map 2d position to sphere with radius 1.
function project_on_ball(x, y) {
  var z = 0;
  var length_sq = x * x + y * y;
  if (length_sq < 1) {  // in ellipse
    z = Math.sqrt(1.0 - length_sq);
  } else {  // in a corner
    var length = Math.sqrt(length_sq);
    x /= length;
    y /= length;
  }
  return [x, y, z];  // guaranteed to be normalized
}

// for two-finger touch events
function touch_info(evt) {
  var touches = evt.touches;
  var dx = touches[0].pageX - touches[1].pageX;
  var dy = touches[0].pageY - touches[1].pageY;
  return {pageX: (touches[0].pageX + touches[1].pageX) / 2,
          pageY: (touches[0].pageY + touches[1].pageY) / 2,
          dist: Math.sqrt(dx * dx + dy * dy)};
}

var STATE = {NONE: -1, ROTATE: 0, PAN: 1, ZOOM: 2, PAN_ZOOM: 3, SLAB: 4,
             ROLL: 5, AUTO_ROTATE: 6, GO: 7};

// based on three.js/examples/js/controls/OrthographicTrackballControls.js
var Controls = function (camera) {
  var _state = STATE.NONE;
  var _rotate_start = new THREE.Vector3();
  var _rotate_end = new THREE.Vector3();
  var _zoom_start = new THREE.Vector2();
  var _zoom_end = new THREE.Vector2();
  var _pinch_start = 0;
  var _pinch_end = 0;
  var _pan_start = new THREE.Vector2();
  var _pan_end = new THREE.Vector2();
  var _panned = true;
  var _slab_width = 10.0;
  var _rock_state = 0.0;
  var _go_func = null;

  function change_slab_width(delta) {
    _slab_width = Math.max(_slab_width + delta, 0.01);
  }

  function rotate_camera(eye) {
    var quat = new THREE.Quaternion();
    quat.setFromUnitVectors(_rotate_end, _rotate_start);
    eye.applyQuaternion(quat);
    camera.up.applyQuaternion(quat);
    _rotate_end.applyQuaternion(quat);
    _rotate_start.copy(_rotate_end);
  }

  function zoom_camera(eye) {
    var dx = _zoom_end.x - _zoom_start.x;
    var dy = _zoom_end.y - _zoom_start.y;
    if (_state === STATE.ZOOM) {
      camera.zoom /= (1 - dx + dy);
    } else if (_state === STATE.SLAB) {
      change_slab_width(10.0 * dx);
      target.addScaledVector(eye, -5.0 / eye.length() * dy);
    } else if (_state === STATE.ROLL) {
      camera.up.applyAxisAngle(eye, 0.05 * (dx - dy));
    }
    _zoom_start.copy(_zoom_end);
  }

  function pan_camera(eye) {
    var dx = _pan_end.x - _pan_start.x;
    var dy = _pan_end.y - _pan_start.y;
    dx *= 0.5 * (camera.right - camera.left) / camera.zoom;
    dy *= 0.5 * (camera.bottom - camera.top) / camera.zoom;
    var pan = eye.clone().cross(camera.up).setLength(dx);
    pan.addScaledVector(camera.up, dy / camera.up.length());
    camera.position.add(pan);
    target.add(pan);
    _pan_start.copy(_pan_end);
  }

  this.toggle_state = function (toggled, params) {
    _state = (_state === toggled ? STATE.NONE : toggled);
    _rock_state = params.rock ? 0.0 : null;
  };

  this.is_going = function () { return _state === STATE.GO; };

  function auto_rotate(eye) {
    _rotate_start.copy(eye).normalize();
    var speed = 3e-4 * auto_speed;
    if (_rock_state !== null) {
      _rock_state += 0.02;
      speed = 4e-5 * auto_speed * Math.cos(_rock_state);
    }
    _rotate_end.crossVectors(camera.up, eye).multiplyScalar(speed)
      .add(_rotate_start);
  }

  this.update = function () {
    var changed = false;
    var eye = camera.position.clone().sub(target);
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
      zoom_camera(eye);
      changed = true;
    }
    if (!_pan_end.equals(_pan_start)) {
      pan_camera(eye);
      _panned = true;
      changed = true;
    }
    camera.position.addVectors(target, eye);
    if (_state === STATE.GO) {
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
      case STATE.ROTATE:
        var xyz = project_on_ball(x, y);
        //console.log(this.camera.projectionMatrix);
        //console.log(this.camera.matrixWorld);
        // TODO maybe use project()/unproject()/applyProjection()
        var eye = camera.position.clone().sub(target);
        _rotate_end.crossVectors(camera.up, eye).setLength(xyz[0]);
        _rotate_end.addScaledVector(camera.up, xyz[1] / camera.up.length());
        _rotate_end.addScaledVector(eye, xyz[2] / eye.length());
        break;
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
    var atom = null;
    if (_state === STATE.PAN && !_panned) {
      atom = pick_atom(_pan_start, camera);
    }
    _state = STATE.NONE;
    _rotate_start.copy(_rotate_end);
    _pinch_start = _pinch_end;
    _pan_start.copy(_pan_end);
    if (atom !== null) { // center on atom
      this.go_to(new THREE.Vector3(atom.xyz[0], atom.xyz[1], atom.xyz[2]));
    }
  };

  this.slab_width = function () { return _slab_width; };
  this.change_slab_width = change_slab_width;

  this.go_to = function (targ, cam_pos, cam_up, steps) {
    if ((!targ || targ.distanceToSquared(target) < 0.1) &&
        (!cam_pos || cam_pos.distanceToSquared(camera.position) < 0.1) &&
        (!cam_up || cam_up.distanceToSquared(camera.up) < 0.1)) {
      return;
    }
    _state = STATE.GO;
    steps = steps || (60 / auto_speed);
    var alphas = [];
    var prev_pos = 0;
    for (var i = 1; i <= steps; ++i) {
      var pos = i / steps;
      // quadratic easing
      pos = pos < 0.5 ? 2 * pos * pos : -2 * pos * (pos-2) - 1;
      alphas.push((pos - prev_pos) / (1 - prev_pos));
      prev_pos = pos;
    }
    _go_func = function () {
      var a = alphas.shift();
      camera.position.sub(target); //XXX
      if (targ) { target.lerp(targ, a); }
      camera.position.add(target); //XXX
      if (cam_pos) { camera.position.lerp(cam_pos, a); }
      if (cam_up) { camera.up.lerp(cam_up, a); }
      if (alphas.length === 0) {
        _state = STATE.NONE;
        _go_func = null;
      }
    };
  };
};

var CUBE_EDGES = [[0, 0, 0], [1, 0, 0],
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

var COLOR_SCHEMES = ['element', 'B-factor', 'occupancy', 'index', 'chain'];
var RENDER_STYLES = ['lines', 'trace', 'lines+balls'];

function rainbow_value(v, vmin, vmax) {
  var c = new THREE.Color(0xe0e0e0);
  if (vmin < vmax) {
    var ratio = (v - vmin) / (vmax - vmin);
    var hue = (240 - (240 * ratio)) / 360;
    c.setHSL(hue, 1.0, 0.5);
  }
  return c;
}

function color_by(style, atoms) {
  var color_func;
  var i;
  var last_atom = atoms[atoms.length-1];
  if (style === 'index') {
    color_func = function (atom) {
      return rainbow_value(atom.i_seq, 0, last_atom.i_seq);
    };
  } else if (style === 'B-factor') {
    var vmin = Infinity;
    var vmax = -Infinity;
    for (i = 0; i < atoms.length; i++) {
      var v = atoms[i].b;
      if (v > vmax) { vmax = v; }
      if (v < vmin) { vmin = v; }
    }
    console.log('B-factors in [' + vmin + ', ' + vmax + ']');
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
      return Colors[atom.element] || Colors.def;
    };
  }
  var colors = [];
  for (i = 0; i < atoms.length; i++) {
    colors.push(color_func(atoms[i]));
  }
  return colors;
}

function make_balls(visible_atoms, colors, ball_size) {
  var ball_texture = new THREE.TextureLoader().load('src/ball.png');
  var pt_geometry = new THREE.Geometry();
  for (var i = 0; i < visible_atoms.length; i++) {
    var xyz = visible_atoms[i].xyz;
    pt_geometry.vertices.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]));
    pt_geometry.colors.push(colors[i]);
  }
  var pt_material = new THREE.PointsMaterial({
    vertexColors: THREE.VertexColors,
    map: ball_texture,
    size: ball_size,
    alphaTest: 0.5
  });
  return new THREE.Points(pt_geometry, pt_material);
}

function make_bonds(model, params, ligands_only, ball_size) {
  //console.time('bonds');
  var visible_atoms = [];
  var all_atoms = model.atoms;
  var i;
  if (!params.hydrogens && model.has_hydrogens) {
    for (i = 0; i < all_atoms.length; i++) {
      if (all_atoms[i].element !== 'H') {
        visible_atoms.push(all_atoms[i]);
      }
    }
  } else {
    visible_atoms = all_atoms;
  }
  var color_style = ligands_only ? 'element' : params.color_scheme;
  var colors = color_by(color_style, visible_atoms);
  var geometry = new THREE.Geometry();
  var opt = { hydrogens: params.hydrogens,
              ligands_only: ligands_only,
              balls: params.render_style === 'lines+balls' };
  for (i = 0; i < visible_atoms.length; i++) {
    var atom = visible_atoms[i];
    var color = colors[i];
    if (ligands_only && !atom.is_ligand) { continue; }
    if (atom.bonds.length === 0 && !opt.balls) { // nonbonded, draw star
      add_isolated_atom(geometry, atom, color);
    } else { // bonded, draw lines
      for (var j = 0; j < atom.bonds.length; j++) {
        // TODO: one line per bond (not trivial, because coloring)
        var other = model.atoms[atom.bonds[j]];
        if (!opt.hydrogens && other.element === 'H') { continue; }
        if (opt.ligands_only && !other.is_ligand) { continue; }
        var mid = atom.midpoint(other);
        var vmid = new THREE.Vector3(mid[0], mid[1], mid[2]);
        var vatom = new THREE.Vector3(atom.xyz[0], atom.xyz[1], atom.xyz[2]);
        if (opt.balls) {
          vatom.lerp(vmid, 0.3); // TODO: use ball_size
        }
        geometry.vertices.push(vatom, vmid);
        geometry.colors.push(color, color);
      }
    }
  }
  //console.timeEnd('bonds');
  var material = new THREE.LineBasicMaterial({
    vertexColors: THREE.VertexColors,
    linewidth: params.line_width
  });
  //console.log('make_bonds() vertex count: ' + geometry.vertices.length);
  var lines = new THREE.LineSegments(geometry, material);
  return opt.balls ? [lines, make_balls(visible_atoms, colors, ball_size)]
                   : [lines];
}

function make_trace(model, params) {
  var segments = model.extract_trace();
  var visible_atoms = [].concat.apply([], segments);
  var colors = color_by(params.color_scheme, visible_atoms);
  var k = 0;
  var geometry = new THREE.Geometry();
  for (var i = 0; i < segments.length; i++) {
    for (var j = 0; j < segments[i].length - 1; j++) {
      var xyz = segments[i][j].xyz;
      var next_xyz = segments[i][j+1].xyz;
      var color = colors[k];
      k++;
      geometry.vertices.push(new THREE.Vector3(xyz[0], xyz[1], xyz[2]),
                    new THREE.Vector3(next_xyz[0], next_xyz[1], next_xyz[2]));
      geometry.colors.push(color, color);
    }
    k++;  // for the last item of segments[i]
  }
  var material = new THREE.LineBasicMaterial({
    vertexColors: THREE.VertexColors,
    linewidth: params.line_width
  });
  return [new THREE.LineSegments(geometry, material)];
}

// Add a representation of an unbonded atom as a cross to geometry
function add_isolated_atom(geometry, atom, color) {
  var c = atom.xyz;
  var R = 0.7;
  geometry.vertices.push(new THREE.Vector3(c[0]-R, c[1], c[2]));
  geometry.vertices.push(new THREE.Vector3(c[0]+R, c[1], c[2]));
  geometry.vertices.push(new THREE.Vector3(c[0], c[1]-R, c[2]));
  geometry.vertices.push(new THREE.Vector3(c[0], c[1]+R, c[2]));
  geometry.vertices.push(new THREE.Vector3(c[0], c[1], c[2]-R));
  geometry.vertices.push(new THREE.Vector3(c[0], c[1], c[2]+R));
  for (var i = 0; i < 6; i++) {
    geometry.colors.push(color);
  }
}

function make_center_cube(size, ctr) {
  var geometry = new THREE.Geometry();
  for (var i = 0; i < CUBE_EDGES.length; i++) {
    var a = CUBE_EDGES[i];
    var x = ctr.x + size * (a[0] - 0.5);
    var y = ctr.y + size * (a[1] - 0.5);
    var z = ctr.z + size * (a[2] - 0.5);
    geometry.vertices.push(new THREE.Vector3(x, y, z));
  }
  var material = new THREE.LineBasicMaterial({color: Colors.center,
                                              linewidth: 2});
  return new THREE.LineSegments(geometry, material);
}

function make_unitcell_box(uc) {
  if (!uc) {
    throw Error('Unit cell not defined!');
  }
  var geometry = new THREE.Geometry();
  var i;
  for (i = 0; i < CUBE_EDGES.length; i += 2) {
    var xyz1 = uc.orthogonalize(CUBE_EDGES[i]);
    var xyz2 = uc.orthogonalize(CUBE_EDGES[i+1]);
    geometry.vertices.push(
      new THREE.Vector3(xyz1[0], xyz1[1], xyz1[2]),
      new THREE.Vector3(xyz2[0], xyz2[1], xyz2[2])
    );
  }
  geometry.colors.push(
    new THREE.Color(0xff0000), new THREE.Color(0xffaa00),
    new THREE.Color(0x00ff00), new THREE.Color(0xaaff00),
    new THREE.Color(0x0000ff), new THREE.Color(0x00aaff)
  );
  for (i = 6; i < CUBE_EDGES.length; i++) {
    geometry.colors.push(Colors.cell_box);
  }
  var material = new THREE.LineBasicMaterial({vertexColors:
                                                THREE.VertexColors});
  return new THREE.LineSegments(geometry, material);
}

function Viewer(element_id) {
  // options
  // on Linux the same width gives different results in Chrome and FF
  this.bond_line = 4.0; // for 700px height (in Coot it also depends on height)
  this.map_line = 1.25;  // for any height
  this.map_radius = 10.0;

  this.model_config = {
    render_style: 'lines',
    color_scheme: 'element',
    hydrogens: false,
    line_width: 0 // it will be set in resize()
  };

  // rendered objects
  this.model_bags = [];
  this.map_bags = [];
  this.decor = {cell_box: null, selection: null};
  this.nav = null;

  this.scene = new THREE.Scene();
  this.camera = new THREE.OrthographicCamera();
  this.scene.add(this.camera);
  this.scene.fog = new THREE.Fog(Colors.bg, 0, 1);
  if (typeof document === 'undefined') { // for testing on node
    return;
  }
  this.renderer = new THREE.WebGLRenderer({antialias: true});
  this.renderer.setClearColor(Colors.bg, 1);
  this.renderer.setPixelRatio(window.devicePixelRatio);
  this.resize();
  this.camera.zoom = this.camera.right / 35.0;
  var container = document.getElementById(element_id);
  container.appendChild(this.renderer.domElement);
  this.controls = new Controls(this.camera);

  this.light = new THREE.AmbientLight(0xffffff);
  this.scene.add(this.light);

  window.addEventListener('resize', this.resize.bind(this));
  window.addEventListener('keydown', this.keydown.bind(this));
  window.addEventListener('contextmenu', function (e) { e.preventDefault(); });
  window.addEventListener('mousewheel', this.mousewheel.bind(this));
  window.addEventListener('MozMousePixelScroll', this.mousewheel.bind(this));
  window.addEventListener('mousedown', this.mousedown.bind(this));
  window.addEventListener('touchstart', this.touchstart.bind(this));
  window.addEventListener('touchmove', this.touchmove.bind(this));
  window.addEventListener('touchend', this.touchend.bind(this));
  window.addEventListener('touchcancel', this.touchend.bind(this));
  window.addEventListener('dblclick', this.dblclick.bind(this));

  this.controls.update();

  var self = this;

  this.mousemove = function (event) {
    event.preventDefault();
    //event.stopPropagation();
    self.controls.move(relX(event), relY(event));
  };

  this.mouseup = function (event) {
    event.preventDefault();
    event.stopPropagation();
    self.controls.stop();
    document.removeEventListener('mousemove', self.mousemove);
    document.removeEventListener('mouseup', self.mouseup);
    self.redraw_maps();
  };
}

Viewer.prototype.toggle_help = function () {
  var el = document.getElementById('help');
  if (el) {
    el.style.display = el.style.display === 'block' ? 'none' : 'block';
  }
};

Viewer.prototype.redraw_center = (function () {
  var prev = new THREE.Vector3(Infinity, 0, 0);
  return function () {
    if (target.distanceToSquared(prev) > 0.0001) {
      prev.copy(target);
      if (this.mark) {
        this.scene.remove(this.mark);
      }
      this.mark = make_center_cube(0.1, target);
      this.scene.add(this.mark);
    }
  };
})();

Viewer.prototype.redraw_maps = function (force) {
  this.redraw_center();
  for (var i = 0; i < this.map_bags.length; i++) {
    var map_bag = this.map_bags[i];
    if (force || target.distanceToSquared(map_bag.block_ctr) > 0.01) {
      this.clear_el_objects(map_bag);
      map_bag.map.block = null;
      this.add_el_objects(map_bag);
    }
  }
};

Viewer.prototype.clear_el_objects = function (map_bag) {
  for (var i = 0; i < map_bag.el_objects.length; i++) {
    this.scene.remove(map_bag.el_objects[i]);
  }
  map_bag.el_objects = [];
};

Viewer.prototype.clear_atomic_objects = function (model) {
  if (model.atomic_objects) {
    for (var i = 0; i < model.atomic_objects.length; i++) {
      this.scene.remove(model.atomic_objects[i]);
    }
  }
  model.atomic_objects = null;
};

Viewer.prototype.set_atomic_objects = function (model_bag) {
  var model = model_bag.model;
  var conf = model_bag.conf;
  switch (conf.render_style) {
    case 'lines':
      model_bag.atomic_objects = make_bonds(model, conf);
      break;
    case 'lines+balls':
      var h_scale = this.camera.projectionMatrix.elements[5];
      var ball_size = Math.max(1, 80 * h_scale);
      model_bag.atomic_objects = make_bonds(model, conf, false, ball_size);
      break;
    case 'trace':  // + lines for ligands
      model_bag.atomic_objects = [].concat(make_trace(model, conf),
                                           make_bonds(model, conf, true));
      break;
  }
  for (var i = 0; i < model_bag.atomic_objects.length; i++) {
    this.scene.add(model_bag.atomic_objects[i]);
  }
  pickable_model = model_bag;
};

Viewer.prototype.toggle_map_visibility = function (map_bag, visible) {
  map_bag.visible = visible;
  if (visible) {
    map_bag.map.block = null;
    this.add_el_objects(map_bag);
  } else {
    this.clear_el_objects(map_bag);
  }
};

Viewer.prototype.toggle_model_visibility = function (model_bag, visible) {
  model_bag.visible = visible;
  if (visible) {
    this.set_atomic_objects(model_bag);
  } else {
    this.clear_atomic_objects(model_bag);
  }
};

Viewer.prototype.redraw_model = function (model_bag) {
  this.clear_atomic_objects(model_bag);
  if (model_bag.visible) {
    this.set_atomic_objects(model_bag);
  }
};

Viewer.prototype.redraw_models = function () {
  for (var i = 0; i < this.model_bags.length; i++) {
    this.redraw_model(this.model_bags[i]);
  }
};

Viewer.prototype.add_el_objects = function (map_bag) {
  if (!map_bag.visible) { return; }
  if (!map_bag.map.block) {
    map_bag.block_ctr.copy(target);
    map_bag.map.extract_block(this.map_radius, [target.x, target.y, target.z]);
  }
  for (var i = 0; i < map_bag.types.length; i++) {
    var mtype = map_bag.types[i];
    var isolevel = (mtype === 'map_neg' ? -1 : 1) * map_bag.isolevel;
    var abs_level = map_bag.map.abs_level(isolevel);
    var bl = map_bag.map.block;
    var geometry = isosurface(bl.points, bl.values, bl.size, abs_level);
    var material = new THREE.MeshBasicMaterial({
      color: Colors[mtype],
      wireframe: true,
      wireframeLinewidth: this.map_line
    });
    var mesh = new THREE.Mesh(geometry, material);
    map_bag.el_objects.push(mesh);
    this.scene.add(mesh);
  }
};

Viewer.prototype.change_isolevel_by = function (map_idx, delta) {
  if (map_idx >= this.map_bags.length) { return; }
  var map_bag = this.map_bags[map_idx];
  map_bag.isolevel += delta;
  var abs_level = map_bag.map.abs_level(map_bag.isolevel);
  hud('map ' + (map_idx+1) + ' level =  ' + abs_level.toFixed(4) +
      'e/A^3 (' + map_bag.isolevel.toFixed(2) + ' rmsd)');
  //TODO: move slow part into update()
  this.clear_el_objects(map_bag);
  this.add_el_objects(map_bag);
};

Viewer.prototype.toggle_cell_box = function () {
  if (this.decor.cell_box) {
    this.scene.remove(this.decor.cell_box);
    this.decor.cell_box = null;
  } else {
    var uc = null;
    if (this.model_bags.length > 0) {
      uc = this.model_bags[0].model.unit_cell;
    }
    // model may not have unit cell
    if (!uc && this.map_bags.length > 0) {
      uc = this.map_bags[0].map.unit_cell;
    }
    if (uc) {
      this.decor.cell_box = make_unitcell_box(uc);
      this.scene.add(this.decor.cell_box);
    }
  }
};

Viewer.prototype.redraw_all = function () {
  this.scene.fog.color = Colors.bg;
  this.renderer.setClearColor(Colors.bg, 1);
  this.redraw_models();
  this.redraw_maps(true);
};

function next(elem, arr) {
  return arr[(arr.indexOf(elem) + 1) % arr.length];
}

Viewer.prototype.keydown = function (evt) {
  var eye;
  var key = evt.keyCode;
  switch (key) {
    case 84:  // t
      this.model_config.render_style = next(this.model_config.render_style,
                                            RENDER_STYLES);
      hud('rendering as ' + this.model_config.render_style);
      this.redraw_models();
      break;
    case 67:  // c
      if (evt.shiftKey) {
        set_colors(next(Colors.name, Object.keys(ColorSchemes)), Colors);
        hud('color scheme: ' + Colors.name);
        this.redraw_all();
      } else { // color-by
        this.model_config.color_scheme = next(this.model_config.color_scheme,
                                              COLOR_SCHEMES);
        hud('coloring by ' + this.model_config.color_scheme);
        this.redraw_models();
      }
      break;
    case 107:  // add
    case 61:  // equals/firefox
    case 187:  // equal sign
      this.change_isolevel_by(evt.shiftKey ? 1 : 0, 0.1);
      break;
    case 109:  // subtract
    case 173:  // minus/firefox
    case 189:  // dash
      this.change_isolevel_by(evt.shiftKey ? 1 : 0, -0.1);
      break;
    case 68:  // d
    case 70:  // f
      this.controls.change_slab_width(key === 68 ? -0.1 : +0.1);
      this.update_camera();
      hud('clip width: ' + (this.camera.far - this.camera.near).toFixed(1));
      break;
    case 77:  // m
    case 78:  // n
      this.camera.zoom *= (key === 77 ? 1.03 : (1 / 1.03));
      this.update_camera();
      hud('zoom: ' + this.camera.zoom.toFixed(2));
      break;
    case 51:  // 3
    case 99:  // numpad 3
    case 108:  // numpad period (Linux)
    case 110:  // decimal point (Mac)
      eye = this.camera.position.clone().sub(target).setLength(1);
      if (key === 108 || key === 110) {
        eye.negate();
      }
      target.add(eye);
      this.camera.position.add(eye);
      this.update_camera();
      this.redraw_maps();
      hud('clip shifted by [' + eye.x.toFixed(2) + ' ' + eye.y.toFixed(2) +
          ' ' + eye.z.toFixed(2) + ']');
      break;
    case 85:  // u
      hud('toggled unit cell box');
      this.toggle_cell_box();
      break;
    case 73:  // i
      hud('toggled camera movement');
      this.controls.toggle_state(STATE.AUTO_ROTATE, {rock: evt.shiftKey});
      break;
    case 82:  // r
      if (evt.shiftKey) {
        hud('redraw!');
        this.redraw_all();
      } else {
        hud('model recentered');
        this.recenter();
      }
      break;
    case 72:  // h
      this.toggle_help();
      break;
    case 36: // Home
    case 35: // End
      this.bond_line += (key === 36 ? 0.2 : -0.2);
      this.bond_line = Math.max(this.bond_line, 0.1);
      this.resize(); // overkill
      hud('bond width: ' + this.model_config.line_width.toFixed(1));
      break;
    case 16: // shift
    case 17: // ctrl
    case 18: // alt
    case 225: // altgr
      break;
    case 32: // Space
      this.center_next_residue(evt.shiftKey);
      break;
    default:
      hud('Nothing here. Press H for help.');
      break;
  }
};

Viewer.prototype.mousedown = function (event) {
  event.preventDefault();
  event.stopPropagation();
  var state = STATE.NONE;
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
      state = STATE.ZOOM;
    }
  }
  this.controls.start(state, relX(event), relY(event));
  document.addEventListener('mousemove', this.mousemove);
  document.addEventListener('mouseup', this.mouseup);
};

Viewer.prototype.dblclick = function (event) {
  if (event.button !== 0) { return; }
  if (this.decor.selection) {
    this.scene.remove(this.decor.selection);
    this.decor.selection = null;
  }
  var mouse = new THREE.Vector2(relX(event), relY(event));
  var atom = pick_atom(mouse, this.camera);
  if (atom !== null) {
    hud(atom.long_label());
    this.set_selection(atom);
  } else {
    hud();
  }
};

Viewer.prototype.set_selection = function (atom) {
  var geometry = new THREE.Geometry();
  geometry.vertices.push(new THREE.Vector3(atom.xyz[0], atom.xyz[1],
                                           atom.xyz[2]));
  var color = Colors[atom.element] || Colors.def;
  var material = new THREE.PointsMaterial({size: 3, color: color});
  this.decor.selection = new THREE.Points(geometry, material);
  this.scene.add(this.decor.selection);
};

Viewer.prototype.touchstart = function (event) {
  var touches = event.touches;
  if (touches.length === 1) {
    this.controls.start(STATE.ROTATE, relX(touches[0]), relY(touches[0]));
  } else { // for now using only two touches
    var info = touch_info(event);
    this.controls.start(STATE.PAN_ZOOM, relX(info), relY(info), info.dist);
  }
};

Viewer.prototype.touchmove = function (event) {
  event.preventDefault();
  event.stopPropagation();
  var touches = event.touches;
  if (touches.length === 1) {
    this.controls.move(relX(touches[0]), relY(touches[0]));
  } else { // for now using only two touches
    var info = touch_info(event);
    this.controls.move(relX(info), relY(info), info.dist);
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
  var delta = evt.wheelDelta ? evt.wheelDelta / 2000
                             : (evt.detail || 0) / -1000;
  this.change_isolevel_by(evt.shiftKey ? 1 : 0, delta);
};

Viewer.prototype.resize = function (/*evt*/) {
  var width = window.innerWidth;
  var height = window.innerHeight;
  this.camera.left = -width;
  this.camera.right = width;
  this.camera.top = height;
  this.camera.bottom = -height;
  this.camera.updateProjectionMatrix();
  this.renderer.setSize(width, height);
  var line_width = this.bond_line * height / 700;
  if (line_width !== this.model_config.line_width) {
    this.model_config.line_width = line_width;
    this.redraw_models();
  }
};

Viewer.prototype.recenter = function (xyz, steps) {
  if (!xyz) { // recenter on the last model
    var len = this.model_bags.length;
    if (len === 0) { return; }
    xyz = this.model_bags[len - 1].model.get_center();
  }
  this.controls.go_to(new THREE.Vector3(xyz[0], xyz[1], xyz[2]),
                      new THREE.Vector3(xyz[0], xyz[1], xyz[2] + 100),
                      new THREE.Vector3(0, 1, 0),
                      steps);
};

Viewer.prototype.center_next_residue = function (back) {
  if (!pickable_model) { return; }
  var a = pickable_model.model.next_residue(selected_atom, back);
  if (a) {
    hud('-> ' + a.long_label());
    this.controls.go_to(new THREE.Vector3(a.xyz[0], a.xyz[1], a.xyz[2]),
                        null, null, 30 / auto_speed);
    selected_atom = a;
  }
};

Viewer.prototype.update_camera = function () {
  var dxyz = this.camera.position.distanceTo(target);
  // the far plane is more distant from the target than the near plane (3:1)
  var w = 0.25 * this.controls.slab_width() / this.camera.zoom;
  this.camera.near = dxyz * (1 - w);
  this.camera.far = dxyz * (1 + 3 * w);
  //this.light.position.copy(this.camera.position);
  var h_scale = this.camera.projectionMatrix.elements[5];
  this.camera.updateProjectionMatrix();
  // temporary hack - scaling balls
  if (h_scale !== this.camera.projectionMatrix.elements[5]) {
    var ball_size = Math.max(1, 80 * this.camera.projectionMatrix.elements[5]);
    for (var i = 0; i < this.model_bags.length; i++) {
      var obj = this.model_bags[i].atomic_objects;
      if (obj.length === 2 && obj[1].material.size) {
        obj[1].material.size = ball_size;
      }
    }
  }
};

Viewer.prototype.render = function render() {
  if (this.controls.update()) {
    this.update_camera();
  }
  if (!this.controls.is_going()) {
    this.redraw_maps();
  }
  this.renderer.render(this.scene, this.camera);
  if (this.nav) {
    this.nav.renderer.render(this.nav.scene, this.camera);
  }
  if (true) { // TODO
    window.requestAnimationFrame(render.bind(this));
  }
};

Viewer.prototype.add_map_bag = function (map, name, is_diff_map) {
  var map_bag = {
    map: map,
    name: name,
    isolevel: is_diff_map ? 3.0 : 1.5,
    visible: true,
    types: is_diff_map ? ['map_pos', 'map_neg'] : ['map_den'],
    block_ctr: new THREE.Vector3(Infinity, 0, 0),
    el_objects: []  // three.js objects
  };
  this.map_bags.push(map_bag);
  //this.add_el_objects(map_bag);
  //this.update_camera();
};

Viewer.prototype.add_model_bag = function (model, name) {
  var model_bag = {
    model: model,
    name: name,
    visible: true,
    conf: this.model_config,
    atomic_objects: null // list of three.js objects
  };
  this.model_bags.push(model_bag);
  this.set_atomic_objects(model_bag);
};

Viewer.prototype.load_pdb = function (url, name) {
  var req = new XMLHttpRequest();
  req.open('GET', url, true);
  var self = this;
  req.onreadystatechange = function () {
    if (req.readyState === 4) {
      // chrome --allow-file-access-from-files gives status 0
      if (req.status === 200 || req.status === 0) {
        var model = new Model();
        model.from_pdb(req.responseText);
        self.add_model_bag(model, name);
        self.recenter(null, 1);
      } else {
        console.log('Error fetching ' + url);
      }
    }
  };
  req.send(null);
};

Viewer.prototype.load_map = function (url, map_name, is_diff_map, filetype) {
  var req = new XMLHttpRequest();
  req.responseType = 'arraybuffer';
  req.open('GET', url, true);
  var self = this;
  req.onreadystatechange = function () {
    if (req.readyState === 4) {
      if (req.status === 200 || req.status === 0) {
        var map = new ElMap();
        if (filetype === 'ccp4') {
          map.from_ccp4(req.response);
        } else if (filetype === 'dsn6') {
          map.from_dsn6(req.response);
        } else {
          throw Error('Unknown map filetype.');
        }
        //map.show_debug_info();
        self.add_map_bag(map, map_name, is_diff_map);
        self.redraw_maps();
      } else {
        console.log('Error fetching ' + url);
      }
    }
  };
  req.send(null);
};

// TODO: navigation window like in gimp and mifit
/*
Viewer.prototype.show_nav = function (inset_id) {
  var inset = document.getElementById(inset_id);
  if (!inset) { return; }
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

return Viewer;
})();

if (typeof module !== 'undefined') { module.exports = Viewer; }
