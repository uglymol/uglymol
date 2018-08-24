// @flow

import { Vector3, Quaternion } from './fromthree.js';

/*:: import type {OrthographicCamera} from './fromthree.js' */

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

export const STATE = { NONE: -1, ROTATE: 0, PAN: 1, ZOOM: 2, PAN_ZOOM: 3,
                       SLAB: 4, ROLL: 5, AUTO_ROTATE: 6, GO: 7 };

const auto_speed = 1.0;

// based on three.js/examples/js/controls/OrthographicTrackballControls.js
export class Controls {
  /*::
    _camera: OrthographicCamera
    _target: Vector3
    _state: number
    _rotate_start: Vector3
    _rotate_end: Vector3
    _zoom_start: [number, number]
    _zoom_end: [number, number]
    _pinch_start: number
    _pinch_end: number
    _pan_start: [number, number]
    _pan_end: [number, number]
    _panned: boolean
    _rotating: number | boolean
    _auto_stamp: number | null
    _go_func:  ?Function
    slab_width: [number, number, ?number];
   */
  constructor(camera /*:OrthographicCamera*/, target /*:Vector3*/) {
    this._camera = camera;
    this._target = target;
    this._state = STATE.NONE;
    this._rotate_start = new Vector3();
    this._rotate_end = new Vector3();
    this._zoom_start = [0, 0];
    this._zoom_end = [0, 0];
    this._pinch_start = 0;
    this._pinch_end = 0;
    this._pan_start = [0, 0];
    this._pan_end = [0, 0];
    this._panned = true;
    this._rotating = 0.0;
    this._auto_stamp = null;
    this._go_func = null;

    // the far plane is more distant from the target than the near plane (3:1)
    this.slab_width = [2.5, 7.5, null];
  }

  _rotate_camera(eye /*:Vector3*/) {
    let quat = new Quaternion();
    quat.setFromUnitVectors(this._rotate_end, this._rotate_start);
    eye.applyQuaternion(quat);
    this._camera.up.applyQuaternion(quat);
    this._rotate_end.applyQuaternion(quat);
    this._rotate_start.copy(this._rotate_end);
  }

  _zoom_camera(eye /*:Vector3*/) {
    const dx = this._zoom_end[0] - this._zoom_start[0];
    const dy = this._zoom_end[1] - this._zoom_start[1];
    if (this._state === STATE.ZOOM) {
      this._camera.zoom /= (1 - dx + dy);
    } else if (this._state === STATE.SLAB) {
      this._target.addScaledVector(eye, -5.0 / eye.length() * dy);
    } else if (this._state === STATE.ROLL) {
      this._camera.up.applyAxisAngle(eye, 0.05 * (dx - dy));
    }
    this._zoom_start[0] = this._zoom_end[0];
    this._zoom_start[1] = this._zoom_end[1];
    return this._state === STATE.SLAB ? 10*dx : null;
  }

  _pan_camera(eye /*:Vector3*/) {
    let dx = this._pan_end[0] - this._pan_start[0];
    let dy = this._pan_end[1] - this._pan_start[1];
    dx *= 0.5 * (this._camera.right - this._camera.left) / this._camera.zoom;
    dy *= 0.5 * (this._camera.bottom - this._camera.top) / this._camera.zoom;
    let pan = eye.clone().cross(this._camera.up).setLength(dx);
    pan.addScaledVector(this._camera.up, dy / this._camera.up.length());
    this._camera.position.add(pan);
    this._target.add(pan);
    this._pan_start[0] = this._pan_end[0];
    this._pan_start[1] = this._pan_end[1];
  }

  _auto_rotate(eye /*:Vector3*/) {
    this._rotate_start.copy(eye).normalize();
    const now = Date.now();
    const elapsed = (this._auto_stamp !== null ? now - this._auto_stamp : 16.7);
    let speed = 1.8e-5 * elapsed * auto_speed;
    this._auto_stamp = now;
    if (this._rotating === true) {
      speed = -speed;
    } else if (this._rotating !== false) {
      this._rotating += 0.02;
      speed = 4e-5 * auto_speed * Math.cos(this._rotating);
    }
    this._rotate_end.crossVectors(this._camera.up, eye).multiplyScalar(speed)
      .add(this._rotate_start);
  }

  toggle_auto(param /*:number|boolean*/) {
    if (this._state === STATE.AUTO_ROTATE &&
        typeof param === typeof this._rotating) {
      this._state = STATE.NONE;
    } else {
      this._state = STATE.AUTO_ROTATE;
      this._auto_stamp = null;
      this._rotating = param;
    }
  }

  is_going() { return this._state === STATE.GO; }

  is_moving() { return this._state !== STATE.NONE; }

  update() {
    let changed = false;
    let eye = this._camera.position.clone().sub(this._target);
    if (this._state === STATE.AUTO_ROTATE) {
      this._auto_rotate(eye);
    }
    if (!this._rotate_start.equals(this._rotate_end)) {
      this._rotate_camera(eye);
      changed = true;
    }
    if (this._pinch_end !== this._pinch_start) {
      this._camera.zoom *= this._pinch_end / this._pinch_start;
      this._pinch_start = this._pinch_end;
      changed = true;
    }
    if (this._zoom_end[0] !== this._zoom_start[0] ||
        this._zoom_end[1] !== this._zoom_start[1]) {
      const dslab = this._zoom_camera(eye);
      if (dslab) {
        this.slab_width[0] = Math.max(this.slab_width[0] + dslab, 0.01);
        this.slab_width[1] = Math.max(this.slab_width[1] + dslab, 0.01);
      }
      changed = true;
    }
    if (this._pan_end[0] !== this._pan_start[0] ||
        this._pan_end[1] !== this._pan_start[1]) {
      this._pan_camera(eye);
      this._panned = true;
      changed = true;
    }
    this._camera.position.addVectors(this._target, eye);
    if (this._state === STATE.GO && this._go_func) {
      this._go_func();
      changed = true;
    }
    this._camera.lookAt(this._target);
    return changed;
  }

  start(new_state /*:number*/, x /*:number*/, y /*:number*/, dist/*:?number*/) {
    if (this._state === STATE.NONE || this._state === STATE.AUTO_ROTATE) {
      this._state = new_state;
    }
    this.move(x, y, dist);
    switch (this._state) {
      case STATE.ROTATE:
        this._rotate_start.copy(this._rotate_end);
        break;
      case STATE.ZOOM:
      case STATE.SLAB:
      case STATE.ROLL:
        this._zoom_start[0] = this._zoom_end[0];
        this._zoom_start[1] = this._zoom_end[1];
        break;
      case STATE.PAN:
        this._pan_start[0] = this._pan_end[0];
        this._pan_start[1] = this._pan_end[1];
        this._panned = false;
        break;
      case STATE.PAN_ZOOM:
        this._pinch_start = this._pinch_end;
        this._pan_start[0] = this._pan_end[0];
        this._pan_start[1] = this._pan_end[1];
        break;
    }
  }

  move(x /*:number*/, y /*:number*/, dist /*:?number*/) {
    switch (this._state) {
      case STATE.ROTATE: {
        const xyz = project_on_ball(x, y);
        //console.log(this._camera.projectionMatrix);
        //console.log(this._camera.matrixWorld);
        // TODO maybe use project()/unproject()/applyProjection()
        const eye = this._camera.position.clone().sub(this._target);
        const up = this._camera.up;
        this._rotate_end.crossVectors(up, eye).setLength(xyz[0]);
        this._rotate_end.addScaledVector(up, xyz[1] / up.length());
        this._rotate_end.addScaledVector(eye, xyz[2] / eye.length());
        break;
      }
      case STATE.ZOOM:
      case STATE.SLAB:
      case STATE.ROLL:
        this._zoom_end = [x, y];
        break;
      case STATE.PAN:
        this._pan_end = [x, y];
        break;
      case STATE.PAN_ZOOM:
        if (dist == null) return; // should not happen
        this._pan_end = [x, y];
        this._pinch_end = dist;
        break;
    }
  }

  // returned coordinates can be used for atom picking
  stop() {
    let ret = null;
    if (this._state === STATE.PAN && !this._panned) ret = this._pan_end;
    this._state = STATE.NONE;
    this._rotate_start.copy(this._rotate_end);
    this._pinch_start = this._pinch_end;
    this._pan_start[0] = this._pan_end[0];
    this._pan_start[1] = this._pan_end[1];
    return ret;
  }

  go_to(targ /*:Vector3*/, cam_pos /*:Vector3*/, cam_up /*:Vector3*/,
        steps /*:?number*/) {
    if (targ instanceof Array) {
      targ = new Vector3(targ[0], targ[1], targ[2]);
    }
    if ((!targ || targ.distanceToSquared(this._target) < 0.001) &&
        (!cam_pos || cam_pos.distanceToSquared(this._camera.position) < 0.1) &&
        (!cam_up || cam_up.distanceToSquared(this._camera.up) < 0.1)) {
      return;
    }
    this._state = STATE.GO;
    steps = (steps || 60) / auto_speed;
    let alphas = [];
    let prev_pos = 0;
    for (let i = 1; i <= steps; ++i) {
      let pos = i / steps;
      // quadratic easing
      pos = pos < 0.5 ? 2 * pos * pos : -2 * pos * (pos-2) - 1;
      alphas.push((pos - prev_pos) / (1 - prev_pos));
      prev_pos = pos;
    }
    this._go_func = function () {
      const a = alphas.shift();
      if (targ) {
        // unspecified cam_pos - _camera stays in the same distance to _target
        if (!cam_pos) this._camera.position.sub(this._target);
        this._target.lerp(targ, a);
        if (!cam_pos) this._camera.position.add(this._target);
      }
      if (cam_pos) this._camera.position.lerp(cam_pos, a);
      if (cam_up) this._camera.up.lerp(cam_up, a);
      if (alphas.length === 0) {
        this._state = STATE.NONE;
        this._go_func = null;
      }
    };
  }
}
