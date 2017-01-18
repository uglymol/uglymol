// @flow

import * as THREE from 'three';

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

// based on three.js/examples/js/controls/OrthographicTrackballControls.js
export class Controls {
  /*::
    slab_width: [number, number]
    go_to: Function
    toggle_auto: Function
    is_going: () => boolean
    is_moving: () => boolean
    update: () => boolean
    start: Function
    move: (number, number, number) => void
    stop: () => any
   */
  constructor(camera /*:THREE.Camera*/, target /*:THREE.Vector3*/) {
    const auto_speed = 1.0;
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
      if ((!targ || targ.distanceToSquared(target) < 0.001) &&
          (!cam_pos || cam_pos.distanceToSquared(camera.position) < 0.1) &&
          (!cam_up || cam_up.distanceToSquared(camera.up) < 0.1)) {
        return;
      }
      _state = STATE.GO;
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
}
