// a small subset of THREE.js v83 that is used by UglyMol
// modified with eslint --fix and manually,
// Copyright 2010-2023 Three.js Authors
// SPDX-License-Identifier: MIT

/* eslint-disable max-len, one-var, guard-for-in */
/* eslint-disable prefer-rest-params, no-invalid-this, no-useless-escape */
/* eslint-disable new-cap, no-extend-native */

// Polyfills

if ( Object.assign === undefined ) {
  // Missing in IE.
  // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Object/assign
  ( function () {
    Object.assign = function ( target ) {
      if ( target === undefined || target === null ) {
        throw new TypeError( 'Cannot convert undefined or null to object' );
      }
      let output = Object( target );
      for ( let index = 1; index < arguments.length; index ++ ) {
        let source = arguments[index];
        if ( source !== undefined && source !== null ) {
          for ( let nextKey in source ) {
            if ( Object.prototype.hasOwnProperty.call( source, nextKey ) ) {
              output[nextKey] = source[nextKey];
            }
          }
        }
      }
      return output;
    };
  } )();
}

let NoBlending = 0;
let NormalBlending = 1;
let TrianglesDrawMode = 0;
let TriangleStripDrawMode = 1;
let TriangleFanDrawMode = 2;

const _lut = [
  '00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '0a', '0b', '0c',
  '0d', '0e', '0f', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
  '1a', '1b', '1c', '1d', '1e', '1f', '20', '21', '22', '23', '24', '25', '26',
  '27', '28', '29', '2a', '2b', '2c', '2d', '2e', '2f', '30', '31', '32', '33',
  '34', '35', '36', '37', '38', '39', '3a', '3b', '3c', '3d', '3e', '3f', '40',
  '41', '42', '43', '44', '45', '46', '47', '48', '49', '4a', '4b', '4c', '4d',
  '4e', '4f', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '5a',
  '5b', '5c', '5d', '5e', '5f', '60', '61', '62', '63', '64', '65', '66', '67',
  '68', '69', '6a', '6b', '6c', '6d', '6e', '6f', '70', '71', '72', '73', '74',
  '75', '76', '77', '78', '79', '7a', '7b', '7c', '7d', '7e', '7f', '80', '81',
  '82', '83', '84', '85', '86', '87', '88', '89', '8a', '8b', '8c', '8d', '8e',
  '8f', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '9a', '9b',
  '9c', '9d', '9e', '9f', 'a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8',
  'a9', 'aa', 'ab', 'ac', 'ad', 'ae', 'af', 'b0', 'b1', 'b2', 'b3', 'b4', 'b5',
  'b6', 'b7', 'b8', 'b9', 'ba', 'bb', 'bc', 'bd', 'be', 'bf', 'c0', 'c1', 'c2',
  'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'ca', 'cb', 'cc', 'cd', 'ce', 'cf',
  'd0', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'da', 'db', 'dc',
  'dd', 'de', 'df', 'e0', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9',
  'ea', 'eb', 'ec', 'ed', 'ee', 'ef', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6',
  'f7', 'f8', 'f9', 'fa', 'fb', 'fc', 'fd', 'fe', 'ff'
];

// http://stackoverflow.com/questions/105034/how-to-create-a-guid-uuid-in-javascript/21963136#21963136
function generateUUID() {
  const d0 = (Math.random() * 0xffffffff) | 0;
  const d1 = (Math.random() * 0xffffffff) | 0;
  const d2 = (Math.random() * 0xffffffff) | 0;
  const d3 = (Math.random() * 0xffffffff) | 0;
  const uuid =
    _lut[d0 & 0xff] + _lut[(d0 >> 8) & 0xff] +
    _lut[(d0 >> 16) & 0xff] + _lut[(d0 >> 24) & 0xff] + '-' +
    _lut[d1 & 0xff] + _lut[(d1 >> 8) & 0xff] + '-' +
    _lut[((d1 >> 16) & 0x0f) | 0x40] + _lut[(d1 >> 24) & 0xff] + '-' +
    _lut[(d2 & 0x3f) | 0x80] + _lut[(d2 >> 8) & 0xff] + '-' +
    _lut[(d2 >> 16) & 0xff] + _lut[(d2 >> 24) & 0xff] + _lut[d3 & 0xff] +
    _lut[(d3 >> 8) & 0xff] + _lut[(d3 >> 16) & 0xff] + _lut[(d3 >> 24) & 0xff];
  // .toLowerCase() here flattens concatenated strings to save heap memory space.
  return uuid.toLowerCase();
}

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

// compute euclidean modulo of m % n
// https://en.wikipedia.org/wiki/Modulo_operation
function euclideanModulo(n, m) {
  return ((n % m) + m) % m;
}


class Quaternion {
  constructor(x = 0, y = 0, z = 0, w = 1) {
    this._x = x;
    this._y = y;
    this._z = z;
    this._w = w;
  }

  get x() {
    return this._x;
  }

  get y() {
    return this._y;
  }

  get z() {
    return this._z;
  }

  get w() {
    return this._w;
  }

  setFromAxisAngle(axis, angle) {
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
    // assumes axis is normalized
    let halfAngle = angle / 2, s = Math.sin(halfAngle);

    this._x = axis.x * s;
    this._y = axis.y * s;
    this._z = axis.z * s;
    this._w = Math.cos(halfAngle);

    return this;
  }

  setFromRotationMatrix(m) {
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
    // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

    const te = m.elements,

      m11 = te[0], m12 = te[4], m13 = te[8],
      m21 = te[1], m22 = te[5], m23 = te[9],
      m31 = te[2], m32 = te[6], m33 = te[10],

      trace = m11 + m22 + m33;

    if (trace > 0) {
      const s = 0.5 / Math.sqrt(trace + 1.0);

      this._w = 0.25 / s;
      this._x = (m32 - m23) * s;
      this._y = (m13 - m31) * s;
      this._z = (m21 - m12) * s;
    } else if (m11 > m22 && m11 > m33) {
      const s = 2.0 * Math.sqrt(1.0 + m11 - m22 - m33);

      this._w = (m32 - m23) / s;
      this._x = 0.25 * s;
      this._y = (m12 + m21) / s;
      this._z = (m13 + m31) / s;
    } else if (m22 > m33) {
      const s = 2.0 * Math.sqrt(1.0 + m22 - m11 - m33);

      this._w = (m13 - m31) / s;
      this._x = (m12 + m21) / s;
      this._y = 0.25 * s;
      this._z = (m23 + m32) / s;
    } else {
      const s = 2.0 * Math.sqrt(1.0 + m33 - m11 - m22);

      this._w = (m21 - m12) / s;
      this._x = (m13 + m31) / s;
      this._y = (m23 + m32) / s;
      this._z = 0.25 * s;
    }

    return this;
  }

  setFromUnitVectors(vFrom, vTo) {
    // assumes direction vectors vFrom and vTo are normalized
    let r = vFrom.dot(vTo) + 1;
    if (r < Number.EPSILON) {
      // vFrom and vTo point in opposite directions
      r = 0;
      if (Math.abs(vFrom.x) > Math.abs(vFrom.z)) {
        this._x = -vFrom.y;
        this._y = vFrom.x;
        this._z = 0;
        this._w = r;
      } else {
        this._x = 0;
        this._y = -vFrom.z;
        this._z = vFrom.y;
        this._w = r;
      }
    } else {
      // crossVectors( vFrom, vTo ); // inlined to avoid cyclic dependency on Vector3
      this._x = vFrom.y * vTo.z - vFrom.z * vTo.y;
      this._y = vFrom.z * vTo.x - vFrom.x * vTo.z;
      this._z = vFrom.x * vTo.y - vFrom.y * vTo.x;
      this._w = r;
    }
    return this.normalize();
  }

  length() {
    return Math.sqrt(this._x * this._x + this._y * this._y +
                     this._z * this._z + this._w * this._w);
  }

  normalize() {
    let l = this.length();
    if (l === 0) {
      this._x = 0;
      this._y = 0;
      this._z = 0;
      this._w = 1;
    } else {
      l = 1 / l;
      this._x = this._x * l;
      this._y = this._y * l;
      this._z = this._z * l;
      this._w = this._w * l;
    }
    return this;
  }
}


class Vector3 {
  constructor(x = 0, y = 0, z = 0) {
    Vector3.prototype.isVector3 = true;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  set(x, y, z) {
    if (z === undefined) z = this.z; // sprite.scale.set(x,y)
    this.x = x;
    this.y = y;
    this.z = z;
    return this;
  }

  clone() {
    return new this.constructor(this.x, this.y, this.z);
  }

  copy(v) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;
    return this;
  }

  add(v) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;
    return this;
  }

  addVectors(a, b) {
    this.x = a.x + b.x;
    this.y = a.y + b.y;
    this.z = a.z + b.z;
    return this;
  }

  addScaledVector(v, s) {
    this.x += v.x * s;
    this.y += v.y * s;
    this.z += v.z * s;
    return this;
  }

  sub(v) {
    this.x -= v.x;
    this.y -= v.y;
    this.z -= v.z;
    return this;
  }

  subVectors(a, b) {
    this.x = a.x - b.x;
    this.y = a.y - b.y;
    this.z = a.z - b.z;
    return this;
  }

  multiplyScalar(scalar) {
    this.x *= scalar;
    this.y *= scalar;
    this.z *= scalar;
    return this;
  }

  applyMatrix4(m) {
    const x = this.x, y = this.y, z = this.z;
    const e = m.elements;
    const w = 1 / (e[3] * x + e[7] * y + e[11] * z + e[15]);
    this.x = (e[0] * x + e[4] * y + e[8] * z + e[12]) * w;
    this.y = (e[1] * x + e[5] * y + e[9] * z + e[13]) * w;
    this.z = (e[2] * x + e[6] * y + e[10] * z + e[14]) * w;
    return this;
  }

  applyQuaternion(q) {
    // quaternion q is assumed to have unit length
    const vx = this.x, vy = this.y, vz = this.z;
    const qx = q.x, qy = q.y, qz = q.z, qw = q.w;
    // t = 2 * cross( q.xyz, v );
    const tx = 2 * (qy * vz - qz * vy);
    const ty = 2 * (qz * vx - qx * vz);
    const tz = 2 * (qx * vy - qy * vx);
    // v + q.w * t + cross( q.xyz, t );
    this.x = vx + qw * tx + qy * tz - qz * ty;
    this.y = vy + qw * ty + qz * tx - qx * tz;
    this.z = vz + qw * tz + qx * ty - qy * tx;
    return this;
  }

  unproject(camera) {
    return this.applyMatrix4(camera.projectionMatrixInverse).applyMatrix4(camera.matrixWorld);
  }

  transformDirection(m) {
    // input: THREE.Matrix4 affine matrix
    // vector interpreted as a direction

    const x = this.x, y = this.y, z = this.z;
    const e = m.elements;

    this.x = e[0] * x + e[4] * y + e[8] * z;
    this.y = e[1] * x + e[5] * y + e[9] * z;
    this.z = e[2] * x + e[6] * y + e[10] * z;

    return this.normalize();
  }

  divideScalar(scalar) {
    return this.multiplyScalar(1 / scalar);
  }

  dot(v) {
    return this.x * v.x + this.y * v.y + this.z * v.z;
  }

  lengthSq() {
    return this.x * this.x + this.y * this.y + this.z * this.z;
  }

  length() {
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
  }

  normalize() {
    return this.divideScalar(this.length() || 1);
  }

  setLength(length) {
    return this.normalize().multiplyScalar(length);
  }

  lerp(v, alpha) {
    this.x += (v.x - this.x) * alpha;
    this.y += (v.y - this.y) * alpha;
    this.z += (v.z - this.z) * alpha;
    return this;
  }

  cross(v) {
    return this.crossVectors(this, v);
  }

  crossVectors(a, b) {
    const ax = a.x, ay = a.y, az = a.z;
    const bx = b.x, by = b.y, bz = b.z;
    this.x = ay * bz - az * by;
    this.y = az * bx - ax * bz;
    this.z = ax * by - ay * bx;
    return this;
  }

  projectOnVector(v) {
    const denominator = v.lengthSq();
    if (denominator === 0) return this.set(0, 0, 0);
    const scalar = v.dot(this) / denominator;
    return this.copy(v).multiplyScalar(scalar);
  }

  projectOnPlane(planeNormal) {
    _vector.copy(this).projectOnVector(planeNormal);
    return this.sub(_vector);
  }

  distanceTo(v) {
    return Math.sqrt(this.distanceToSquared(v));
  }

  distanceToSquared(v) {
    const dx = this.x - v.x, dy = this.y - v.y, dz = this.z - v.z;
    return dx * dx + dy * dy + dz * dz;
  }

  setFromMatrixPosition(m) {
    const e = m.elements;
    this.x = e[12];
    this.y = e[13];
    this.z = e[14];
    return this;
  }

  setFromMatrixColumn(m, index) {
    return this.fromArray(m.elements, index * 4);
  }

  equals(v) {
    return v.x === this.x && v.y === this.y && v.z === this.z;
  }

  fromArray(array, offset = 0) {
    this.x = array[offset];
    this.y = array[offset + 1];
    this.z = array[offset + 2];
    return this;
  }
}
const _vector = /*@__PURE__*/ new Vector3();


class Vector4 {
  constructor(x = 0, y = 0, z = 0, w = 1) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.w = w;
  }

  set(x, y, z, w) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.w = w;
    return this;
  }

  copy(v) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;
    this.w = v.w !== undefined ? v.w : 1;
    return this;
  }

  multiplyScalar(scalar) {
    this.x *= scalar;
    this.y *= scalar;
    this.z *= scalar;
    this.w *= scalar;
    return this;
  }

  equals(v) {
    return v.x === this.x && v.y === this.y && v.z === this.z && v.w === this.w;
  }
}


/**
* @author mrdoob / http://mrdoob.com/
* @author supereggbert / http://www.paulbrunt.co.uk/
* @author philogb / http://blog.thejit.org/
* @author jordi_ros / http://plattsoft.com
* @author D1plo1d / http://github.com/D1plo1d
* @author alteredq / http://alteredqualia.com/
* @author mikael emtinger / http://gomo.se/
* @author timknip / http://www.floorplanner.com/
* @author bhouston / http://clara.io
* @author WestLangley / http://github.com/WestLangley
*/

function Matrix4() {
  this.elements = new Float32Array( [
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1,
  ] );
}

Matrix4.prototype = {

  constructor: Matrix4,

  isMatrix4: true,


  copy: function ( m ) {
    this.elements.set( m.elements );

    return this;
  },

  makeRotationFromQuaternion: function ( q ) {
    let te = this.elements;

    let x = q.x, y = q.y, z = q.z, w = q.w;
    let x2 = x + x, y2 = y + y, z2 = z + z;
    let xx = x * x2, xy = x * y2, xz = x * z2;
    let yy = y * y2, yz = y * z2, zz = z * z2;
    let wx = w * x2, wy = w * y2, wz = w * z2;

    te[0] = 1 - ( yy + zz );
    te[4] = xy - wz;
    te[8] = xz + wy;

    te[1] = xy + wz;
    te[5] = 1 - ( xx + zz );
    te[9] = yz - wx;

    te[2] = xz - wy;
    te[6] = yz + wx;
    te[10] = 1 - ( xx + yy );

    // last column
    te[3] = 0;
    te[7] = 0;
    te[11] = 0;

    // bottom row
    te[12] = 0;
    te[13] = 0;
    te[14] = 0;
    te[15] = 1;

    return this;
  },

  lookAt: function () {
    let x, y, z;

    return function lookAt( eye, target, up ) {
      if ( x === undefined ) {
        x = new Vector3();
        y = new Vector3();
        z = new Vector3();
      }

      let te = this.elements;

      z.subVectors( eye, target ).normalize();

      if ( z.lengthSq() === 0 ) {
        z.z = 1;
      }

      x.crossVectors( up, z ).normalize();

      if ( x.lengthSq() === 0 ) {
        z.z += 0.0001;
        x.crossVectors( up, z ).normalize();
      }

      y.crossVectors( z, x );


      te[0] = x.x; te[4] = y.x; te[8] = z.x;
      te[1] = x.y; te[5] = y.y; te[9] = z.y;
      te[2] = x.z; te[6] = y.z; te[10] = z.z;

      return this;
    };
  }(),

  multiplyMatrices: function ( a, b ) {
    let ae = a.elements;
    let be = b.elements;
    let te = this.elements;

    let a11 = ae[0], a12 = ae[4], a13 = ae[8], a14 = ae[12];
    let a21 = ae[1], a22 = ae[5], a23 = ae[9], a24 = ae[13];
    let a31 = ae[2], a32 = ae[6], a33 = ae[10], a34 = ae[14];
    let a41 = ae[3], a42 = ae[7], a43 = ae[11], a44 = ae[15];

    let b11 = be[0], b12 = be[4], b13 = be[8], b14 = be[12];
    let b21 = be[1], b22 = be[5], b23 = be[9], b24 = be[13];
    let b31 = be[2], b32 = be[6], b33 = be[10], b34 = be[14];
    let b41 = be[3], b42 = be[7], b43 = be[11], b44 = be[15];

    te[0] = a11 * b11 + a12 * b21 + a13 * b31 + a14 * b41;
    te[4] = a11 * b12 + a12 * b22 + a13 * b32 + a14 * b42;
    te[8] = a11 * b13 + a12 * b23 + a13 * b33 + a14 * b43;
    te[12] = a11 * b14 + a12 * b24 + a13 * b34 + a14 * b44;

    te[1] = a21 * b11 + a22 * b21 + a23 * b31 + a24 * b41;
    te[5] = a21 * b12 + a22 * b22 + a23 * b32 + a24 * b42;
    te[9] = a21 * b13 + a22 * b23 + a23 * b33 + a24 * b43;
    te[13] = a21 * b14 + a22 * b24 + a23 * b34 + a24 * b44;

    te[2] = a31 * b11 + a32 * b21 + a33 * b31 + a34 * b41;
    te[6] = a31 * b12 + a32 * b22 + a33 * b32 + a34 * b42;
    te[10] = a31 * b13 + a32 * b23 + a33 * b33 + a34 * b43;
    te[14] = a31 * b14 + a32 * b24 + a33 * b34 + a34 * b44;

    te[3] = a41 * b11 + a42 * b21 + a43 * b31 + a44 * b41;
    te[7] = a41 * b12 + a42 * b22 + a43 * b32 + a44 * b42;
    te[11] = a41 * b13 + a42 * b23 + a43 * b33 + a44 * b43;
    te[15] = a41 * b14 + a42 * b24 + a43 * b34 + a44 * b44;

    return this;
  },

  setPosition: function ( v ) {
    let te = this.elements;

    te[12] = v.x;
    te[13] = v.y;
    te[14] = v.z;

    return this;
  },

  invert: function () {
    // based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
    const te = this.elements,

      n11 = te[0], n21 = te[1], n31 = te[2], n41 = te[3],
      n12 = te[4], n22 = te[5], n32 = te[6], n42 = te[7],
      n13 = te[8], n23 = te[9], n33 = te[10], n43 = te[11],
      n14 = te[12], n24 = te[13], n34 = te[14], n44 = te[15],

      t11 = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44,
      t12 = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44,
      t13 = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44,
      t14 = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;

    const det = n11 * t11 + n21 * t12 + n31 * t13 + n41 * t14;

    if (det === 0) return this.set(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    const detInv = 1 / det;

    te[0] = t11 * detInv;
    te[1] = ( n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 + n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44 ) * detInv;
    te[2] = ( n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 - n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44 ) * detInv;
    te[3] = ( n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 + n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43 ) * detInv;

    te[4] = t12 * detInv;
    te[5] = ( n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 - n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44 ) * detInv;
    te[6] = ( n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 + n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44 ) * detInv;
    te[7] = ( n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 - n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43 ) * detInv;

    te[8] = t13 * detInv;
    te[9] = ( n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 + n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44 ) * detInv;
    te[10] = ( n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 - n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44 ) * detInv;
    te[11] = ( n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 + n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43 ) * detInv;

    te[12] = t14 * detInv;
    te[13] = ( n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 - n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34 ) * detInv;
    te[14] = ( n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 + n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34 ) * detInv;
    te[15] = ( n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33 ) * detInv;

    return this;
  },

  scale: function ( v ) {
    let te = this.elements;
    let x = v.x, y = v.y, z = v.z;

    te[0] *= x; te[4] *= y; te[8] *= z;
    te[1] *= x; te[5] *= y; te[9] *= z;
    te[2] *= x; te[6] *= y; te[10] *= z;
    te[3] *= x; te[7] *= y; te[11] *= z;

    return this;
  },

  getMaxScaleOnAxis: function () {
    let te = this.elements;

    let scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
    let scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
    let scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];

    return Math.sqrt( Math.max( scaleXSq, scaleYSq, scaleZSq ) );
  },

  compose: function ( position, quaternion, scale ) {
    this.makeRotationFromQuaternion( quaternion );
    this.scale( scale );
    this.setPosition( position );

    return this;
  },

  makeOrthographic: function ( left, right, top, bottom, near, far ) {
    let te = this.elements;
    let w = 1.0 / ( right - left );
    let h = 1.0 / ( top - bottom );
    let p = 1.0 / ( far - near );

    let x = ( right + left ) * w;
    let y = ( top + bottom ) * h;
    let z = ( far + near ) * p;

    te[0] = 2 * w; te[4] = 0; te[8] = 0; te[12] = - x;
    te[1] = 0; te[5] = 2 * h; te[9] = 0; te[13] = - y;
    te[2] = 0; te[6] = 0; te[10] = - 2 * p; te[14] = - z;
    te[3] = 0; te[7] = 0; te[11] = 0; te[15] = 1;

    return this;
  },
};

/**
* https://github.com/mrdoob/eventdispatcher.js/
*/

class EventDispatcher {
  addEventListener(type, listener) {
    if (this._listeners === undefined) this._listeners = {};

    const listeners = this._listeners;

    if (listeners[type] === undefined) {
      listeners[type] = [];
    }

    if (listeners[type].indexOf(listener) === -1) {
      listeners[type].push(listener);
    }
  }

  removeEventListener(type, listener) {
    if (this._listeners === undefined) return;

    const listeners = this._listeners;
    const listenerArray = listeners[type];

    if (listenerArray !== undefined) {
      const index = listenerArray.indexOf(listener);

      if (index !== -1) {
        listenerArray.splice(index, 1);
      }
    }
  }

  dispatchEvent(event) {
    if (this._listeners === undefined) return;

    const listeners = this._listeners;
    const listenerArray = listeners[event.type];

    if (listenerArray !== undefined) {
      event.target = this;

      // Make a copy, in case listeners are removed while iterating.
      const array = listenerArray.slice(0);

      for (let i = 0, l = array.length; i < l; i++) {
        array[i].call(this, event);
      }

      event.target = null;
    }
  }
}

let _sourceId = 0;
class Source {
  constructor(data = null) {
    Object.defineProperty(this, 'id', { value: _sourceId++ });
    this.uuid = generateUUID();
    this.data = data;
    this.dataReady = true;
    this.version = 0;
  }
  set needsUpdate( value ) {
    if (value === true) this.version++;
  }
}

let _textureId = 0;
class Texture extends EventDispatcher {
  constructor(image) {
    super();
    Object.defineProperty(this, 'id', { value: _textureId++ });
    this.uuid = generateUUID();
    this.name = '';
    this.source = new Source(image);
    this.version = 0;
  }

  get image() {
    return this.source.data;
  }

  dispose() {
    this.dispatchEvent({ type: 'dispose' });
  }

  set needsUpdate(value) {
    if (value === true) {
      this.version++;
      this.source.needsUpdate = true;
    }
  }
}


/**
 * Uniforms of a program.
 * Those form a tree structure with a special top-level container for the root,
 * which you get by calling 'new WebGLUniforms( gl, program )'.
 *
 *
 * Properties of inner nodes including the top-level container:
 *
 * .seq - array of nested uniforms
 * .map - nested uniforms by name
 *
 *
 * Methods of all nodes except the top-level container:
 *
 * .setValue( gl, value, [textures] )
 *
 *              uploads a uniform value(s)
 *      the 'textures' parameter is needed for sampler uniforms
 *
 *
 * Static methods of the top-level container (textures factorizations):
 *
 * .upload( gl, seq, values, textures )
 *
 *              sets uniforms in 'seq' to 'values[id].value'
 *
 * .seqWithValue( seq, values ) : filteredSeq
 *
 *              filters 'seq' entries with corresponding entry in values
 *
 *
 * Methods of the top-level container (textures factorizations):
 *
 * .setValue( gl, name, value, textures )
 *
 *              sets uniform with  name 'name' to 'value'
 *
 * .setOptional( gl, obj, prop )
 *
 *              like .set for an optional property of the object
 *
 */
const emptyTexture = /*@__PURE__*/ new Texture();

// --- Base for inner nodes (including the root) ---

// TODO
function UniformContainer() {
  this.seq = [];
  this.map = {};
}

// --- Setters ---

// Note: Defining these methods externally, because they come in a bunch
// and this way their names minify.

// Single scalar

function setValue1f( gl, v ) { gl.uniform1f( this.addr, v ); }
function setValue1i( gl, v ) { gl.uniform1i( this.addr, v ); }

// Single float vector (from flat array or VectorN)

function setValue2fv( gl, v ) {
  if ( v.x === undefined ) gl.uniform2fv( this.addr, v );
  else gl.uniform2f( this.addr, v.x, v.y );
}

function setValue3fv( gl, v ) {
  if ( v.x !== undefined ) {
    gl.uniform3f( this.addr, v.x, v.y, v.z );
  } else if ( v.r !== undefined ) {
    gl.uniform3f( this.addr, v.r, v.g, v.b );
  } else {
    gl.uniform3fv( this.addr, v );
  }
}

function setValue4fv( gl, v ) {
  if ( v.x === undefined ) gl.uniform4fv( this.addr, v );
  else gl.uniform4f( this.addr, v.x, v.y, v.z, v.w );
}

// Single matrix (from flat array or MatrixN)

function setValue2fm( gl, v ) {
  gl.uniformMatrix2fv( this.addr, false, v.elements || v );
}

function setValue3fm( gl, v ) {
  gl.uniformMatrix3fv( this.addr, false, v.elements || v );
}

function setValue4fm( gl, v ) {
  gl.uniformMatrix4fv( this.addr, false, v.elements || v );
}

// Single texture (2D / Cube)

function setValueT1( gl, v, renderer ) {
  let unit = renderer.allocTextureUnit();
  gl.uniform1i( this.addr, unit );
  renderer.setTexture2D( v || emptyTexture, unit );
}

// Integer / Boolean vectors or arrays thereof (always flat arrays)

function setValue2iv( gl, v ) { gl.uniform2iv( this.addr, v ); }
function setValue3iv( gl, v ) { gl.uniform3iv( this.addr, v ); }
function setValue4iv( gl, v ) { gl.uniform4iv( this.addr, v ); }

// Helper to pick the right setter for the singular case

function getSingularSetter( type ) {
  switch ( type ) {
    case 0x1406: return setValue1f; // FLOAT
    case 0x8b50: return setValue2fv; // _VEC2
    case 0x8b51: return setValue3fv; // _VEC3
    case 0x8b52: return setValue4fv; // _VEC4

    case 0x8b5a: return setValue2fm; // _MAT2
    case 0x8b5b: return setValue3fm; // _MAT3
    case 0x8b5c: return setValue4fm; // _MAT4

    case 0x8b5e: return setValueT1; // SAMPLER_2D

    case 0x1404: case 0x8b56: return setValue1i; // INT, BOOL
    case 0x8b53: case 0x8b57: return setValue2iv; // _VEC2
    case 0x8b54: case 0x8b58: return setValue3iv; // _VEC3
    case 0x8b55: case 0x8b59: return setValue4iv; // _VEC4
  }
}


// --- Uniform Classes ---

class SingleUniform {
  constructor(id, activeInfo, addr) {
    this.id = id;
    this.addr = addr;
    //this.cache = [];
    //this.type = activeInfo.type;
    this.setValue = getSingularSetter(activeInfo.type);
    // this.path = activeInfo.name; // DEBUG
  }
}

// --- Top-level ---

// Parser - builds up the property tree from the path strings

const RePathPart = /(\w+)(\])?(\[|\.)?/g;

// extracts
//  - the identifier (member name or array index)
//  - followed by an optional right bracket (found when array index)
//  - followed by an optional left bracket or dot (type of subscript)
//
// Note: These portions can be read in a non-overlapping fashion and
// allow straightforward parsing of the hierarchy that WebGL encodes
// in the uniform names.

function addUniform(container, uniformObject) {
  container.seq.push( uniformObject );
  container.map[uniformObject.id] = uniformObject;
}

function parseUniform( activeInfo, addr, container ) {
  let path = activeInfo.name;

  // reset RegExp object, because of the early exit of a previous run
  RePathPart.lastIndex = 0;

  for (; ;) {
    let match = RePathPart.exec( path ),
      id = match[1],
      idIsIndex = match[2] === ']',
      subscript = match[3];

    if ( idIsIndex ) id = id | 0; // convert to integer

    if ( subscript === undefined ) {
      addUniform( container, new SingleUniform( id, activeInfo, addr ) );
      break;
    }
  }
}

// Root Container

function WebGLUniforms( gl, program, renderer ) {
  UniformContainer.call( this );

  this.renderer = renderer;

  let n = gl.getProgramParameter( program, gl.ACTIVE_UNIFORMS );

  for ( let i = 0; i !== n; ++ i ) {
    let info = gl.getActiveUniform( program, i ),
      path = info.name,
      addr = gl.getUniformLocation( program, path );

    parseUniform( info, addr, this );
  }
}

WebGLUniforms.prototype.setValue = function ( gl, name, value ) {
  let u = this.map[name];
  if ( u !== undefined ) u.setValue( gl, value, this.renderer );
};

WebGLUniforms.prototype.set = function ( gl, object, name ) {
  let u = this.map[name];
  if ( u !== undefined ) u.setValue( gl, object[name], this.renderer );
};

// Static interface

WebGLUniforms.upload = function ( gl, seq, values, renderer ) {
  for ( let i = 0, n = seq.length; i !== n; ++ i ) {
    let u = seq[i],
      v = values[u.id];

    if ( v.needsUpdate !== false ) {
      // note: always updating when .needsUpdate is undefined

      u.setValue( gl, v.value, renderer );
    }
  }
};

WebGLUniforms.seqWithValue = function ( seq, values ) {
  let r = [];

  for ( let i = 0, n = seq.length; i !== n; ++ i ) {
    let u = seq[i];
    if ( u.id in values ) r.push( u );
  }

  return r;
};


function hue2rgb(p, q, t) {
  if (t < 0) t += 1;
  if (t > 1) t -= 1;
  if (t < 1 / 6) return p + (q - p) * 6 * t;
  if (t < 1 / 2) return q;
  if (t < 2 / 3) return p + (q - p) * 6 * (2 / 3 - t);
  return p;
}

class Color {
  constructor(r, g, b) {
    this.isColor = true;
    this.r = 1;
    this.g = 1;
    this.b = 1;
    return this.set(r, g, b);
  }

  set(r, g, b) {
    if (g === undefined && b === undefined) {
      // r is THREE.Color, hex or string
      const value = r;
      if (value && value.isColor) {
        this.copy(value);
      } else if (typeof value === 'number') {
        this.setHex(value);
      }
    } else {
      this.setRGB(r, g, b);
    }
    return this;
  }

  setHex(hex) {
    hex = Math.floor(hex);
    this.r = ((hex >> 16) & 255) / 255;
    this.g = ((hex >> 8) & 255) / 255;
    this.b = (hex & 255) / 255;
    return this;
  }

  setRGB(r, g, b) {
    this.r = r;
    this.g = g;
    this.b = b;
    return this;
  }

  setHSL(h, s, l) {
    // h,s,l ranges are in 0.0 - 1.0
    h = euclideanModulo(h, 1);
    s = clamp(s, 0, 1);
    l = clamp(l, 0, 1);

    if (s === 0) {
      this.r = this.g = this.b = l;
    } else {
      const p = l <= 0.5 ? l * (1 + s) : l + s - l * s;
      const q = 2 * l - p;

      this.r = hue2rgb(q, p, h + 1 / 3);
      this.g = hue2rgb(q, p, h);
      this.b = hue2rgb(q, p, h - 1 / 3);
    }
    return this;
  }

  getHSL(target) {
    // h,s,l ranges are in 0.0 - 1.0

    const r = this.r, g = this.g, b = this.b;

    const max = Math.max(r, g, b);
    const min = Math.min(r, g, b);

    let hue, saturation;
    const lightness = (min + max) / 2.0;

    if (min === max) {
      hue = 0;
      saturation = 0;
    } else {
      const delta = max - min;

      saturation = lightness <= 0.5 ? delta / (max + min) : delta / (2 - max - min);

      switch (max) {
        case r:
          hue = (g - b) / delta + (g < b ? 6 : 0);
          break;
        case g:
          hue = (b - r) / delta + 2;
          break;
        case b:
          hue = (r - g) / delta + 4;
          break;
      }

      hue /= 6;
    }
    target.h = hue;
    target.s = saturation;
    target.l = lightness;

    return target;
  }

  clone() {
    return new this.constructor(this.r, this.g, this.b);
  }

  copy(color) {
    this.r = color.r;
    this.g = color.g;
    this.b = color.b;

    return this;
  }

  getHex() {
    return ( this.r * 255 ) << 16 ^ ( this.g * 255 ) << 8 ^ ( this.b * 255 ) << 0;
  }

  getHexString() {
    return ('000000' + this.getHex().toString(16)).slice(-6);
  }
}

/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
*/

let materialId = 0;

function Material() {
  Object.defineProperty( this, 'id', { value: materialId ++ } );

  this.uuid = generateUUID();

  this.name = '';
  this.type = 'Material';

  this.fog = true;

  this.opacity = 1;
  this.transparent = false;

  this.depthTest = true;
  this.depthWrite = true;

  this.precision = null; // override the renderer's default precision for this material

  this.premultipliedAlpha = false;

  this.visible = true;

  this._needsUpdate = true;
}

Material.prototype = {

  constructor: Material,

  isMaterial: true,

  get needsUpdate() {
    return this._needsUpdate;
  },

  set needsUpdate( value ) {
    if ( value === true ) this.update();
    this._needsUpdate = value;
  },

  setValues: function ( values ) {
    for ( let key in values ) {
      let newValue = values[key];
      this[key] = newValue;
    }
  },

  update: function () {
    this.dispatchEvent( { type: 'update' } );
  },

  dispose: function () {
    this.dispatchEvent( { type: 'dispose' } );
  },

};

Object.assign( Material.prototype, EventDispatcher.prototype );

/**
* @author alteredq / http://alteredqualia.com/
*
* parameters = {
*  uniforms: { "parameter1": { value: 1.0 }, "parameter2": { value2: 2 } },
*
*  fragmentShader: <string>,
*  vertexShader: <string>,
* }
*/

function ShaderMaterial( parameters ) {
  Material.call( this );

  this.type = 'ShaderMaterial';

  this.uniforms = {};

  this.vertexShader = '';
  this.fragmentShader = '';

  this.linewidth = 1;

  this.fog = false; // set to use scene fog

  this.extensions = {
    fragDepth: false, // set to use fragment depth values
  };

  this.setValues( parameters );
}

ShaderMaterial.prototype = Object.create( Material.prototype );
ShaderMaterial.prototype.constructor = ShaderMaterial;

ShaderMaterial.prototype.isShaderMaterial = true;


/**
* @author bhouston / http://clara.io
*/

function Ray( origin, direction ) {
  this.origin = ( origin !== undefined ) ? origin : new Vector3();
  this.direction = ( direction !== undefined ) ? direction : new Vector3();
}

Ray.prototype = {

  constructor: Ray,

  copy: function ( ray ) {
    this.origin.copy( ray.origin );
    this.direction.copy( ray.direction );

    return this;
  },

  distanceSqToPoint: function () {
    let v1 = new Vector3();

    return function distanceSqToPoint( point ) {
      let directionDistance = v1.subVectors( point, this.origin ).dot( this.direction );
      // point behind the ray
      if ( directionDistance < 0 ) {
        return this.origin.distanceToSquared( point );
      }
      v1.copy( this.direction ).multiplyScalar( directionDistance ).add( this.origin );
      return v1.distanceToSquared( point );
    };
  }(),

  distanceSqToSegment: function () {
    let segCenter = new Vector3();
    let segDir = new Vector3();
    let diff = new Vector3();

    return function distanceSqToSegment( v0, v1, optionalPointOnRay, optionalPointOnSegment ) {
      // from http://www.geometrictools.com/GTEngine/Include/Mathematics/GteDistRaySegment.h
      // It returns the min distance between the ray and the segment
      // defined by v0 and v1
      // It can also set two optional targets :
      // - The closest point on the ray
      // - The closest point on the segment

      segCenter.copy( v0 ).add( v1 ).multiplyScalar( 0.5 );
      segDir.copy( v1 ).sub( v0 ).normalize();
      diff.copy( this.origin ).sub( segCenter );

      let segExtent = v0.distanceTo( v1 ) * 0.5;
      let a01 = - this.direction.dot( segDir );
      let b0 = diff.dot( this.direction );
      let b1 = - diff.dot( segDir );
      let c = diff.lengthSq();
      let det = Math.abs( 1 - a01 * a01 );
      let s0, s1, sqrDist, extDet;

      if ( det > 0 ) {
        // The ray and segment are not parallel.

        s0 = a01 * b1 - b0;
        s1 = a01 * b0 - b1;
        extDet = segExtent * det;

        if ( s0 >= 0 ) {
          if ( s1 >= - extDet ) {
            if ( s1 <= extDet ) {
              // region 0
              // Minimum at interior points of ray and segment.

              let invDet = 1 / det;
              s0 *= invDet;
              s1 *= invDet;
              sqrDist = s0 * ( s0 + a01 * s1 + 2 * b0 ) + s1 * ( a01 * s0 + s1 + 2 * b1 ) + c;
            } else {
              // region 1

              s1 = segExtent;
              s0 = Math.max( 0, - ( a01 * s1 + b0 ) );
              sqrDist = - s0 * s0 + s1 * ( s1 + 2 * b1 ) + c;
            }
          } else {
            // region 5

            s1 = - segExtent;
            s0 = Math.max( 0, - ( a01 * s1 + b0 ) );
            sqrDist = - s0 * s0 + s1 * ( s1 + 2 * b1 ) + c;
          }
        } else {
          if ( s1 <= - extDet ) {
            // region 4

            s0 = Math.max( 0, - ( - a01 * segExtent + b0 ) );
            s1 = ( s0 > 0 ) ? - segExtent : Math.min( Math.max( - segExtent, - b1 ), segExtent );
            sqrDist = - s0 * s0 + s1 * ( s1 + 2 * b1 ) + c;
          } else if ( s1 <= extDet ) {
            // region 3

            s0 = 0;
            s1 = Math.min( Math.max( - segExtent, - b1 ), segExtent );
            sqrDist = s1 * ( s1 + 2 * b1 ) + c;
          } else {
            // region 2

            s0 = Math.max( 0, - ( a01 * segExtent + b0 ) );
            s1 = ( s0 > 0 ) ? segExtent : Math.min( Math.max( - segExtent, - b1 ), segExtent );
            sqrDist = - s0 * s0 + s1 * ( s1 + 2 * b1 ) + c;
          }
        }
      } else {
        // Ray and segment are parallel.

        s1 = ( a01 > 0 ) ? - segExtent : segExtent;
        s0 = Math.max( 0, - ( a01 * s1 + b0 ) );
        sqrDist = - s0 * s0 + s1 * ( s1 + 2 * b1 ) + c;
      }

      if ( optionalPointOnRay ) {
        optionalPointOnRay.copy( this.direction ).multiplyScalar( s0 ).add( this.origin );
      }

      if ( optionalPointOnSegment ) {
        optionalPointOnSegment.copy( segDir ).multiplyScalar( s1 ).add( segCenter );
      }

      return sqrDist;
    };
  }(),

  applyMatrix4: function ( matrix4 ) {
    this.direction.add( this.origin ).applyMatrix4( matrix4 );
    this.origin.applyMatrix4( matrix4 );
    this.direction.sub( this.origin );
    this.direction.normalize();

    return this;
  },

};

/**
* @author mrdoob / http://mrdoob.com/
* @author mikael emtinger / http://gomo.se/
* @author alteredq / http://alteredqualia.com/
* @author WestLangley / http://github.com/WestLangley
* @author elephantatwork / www.elephantatwork.ch
*/

let object3DId = 0;

function Object3D() {
  Object.defineProperty( this, 'id', { value: object3DId ++ } );

  this.uuid = generateUUID();

  this.name = '';
  this.type = 'Object3D';

  this.parent = null;
  this.children = [];

  this.up = Object3D.DefaultUp.clone();

  let position = new Vector3();
  let quaternion = new Quaternion();
  let scale = new Vector3( 1, 1, 1 );


  Object.defineProperties( this, {
    position: {
      enumerable: true,
      value: position,
    },
    quaternion: {
      enumerable: true,
      value: quaternion,
    },
    scale: {
      enumerable: true,
      value: scale,
    },
    modelViewMatrix: {
      value: new Matrix4(),
    },
  } );

  this.matrix = new Matrix4();
  this.matrixWorld = new Matrix4();

  this.matrixAutoUpdate = Object3D.DefaultMatrixAutoUpdate;
  this.matrixWorldNeedsUpdate = false;

  this.visible = true;

  this.frustumCulled = true;
  this.renderOrder = 0;

  this.userData = {};
}

Object3D.DefaultUp = new Vector3( 0, 1, 0 );
Object3D.DefaultMatrixAutoUpdate = true;

Object.assign( Object3D.prototype, EventDispatcher.prototype, {

  isObject3D: true,

  add: function ( object ) {
    if ( arguments.length > 1 ) {
      for ( let i = 0; i < arguments.length; i ++ ) {
        this.add( arguments[i] );
      }

      return this;
    }

    if ( ( object && object.isObject3D ) ) {
      if ( object.parent !== null ) {
        object.parent.remove( object );
      }

      object.parent = this;
      object.dispatchEvent( { type: 'added' } );

      this.children.push( object );
    }

    return this;
  },

  remove: function ( object ) {
    if ( arguments.length > 1 ) {
      for ( let i = 0; i < arguments.length; i ++ ) {
        this.remove( arguments[i] );
      }
    }

    let index = this.children.indexOf( object );

    if ( index !== - 1 ) {
      object.parent = null;

      object.dispatchEvent( { type: 'removed' } );

      this.children.splice( index, 1 );
    }
  },

  updateMatrix: function () {
    this.matrix.compose( this.position, this.quaternion, this.scale );

    this.matrixWorldNeedsUpdate = true;
  },

  updateMatrixWorld: function ( force ) {
    if ( this.matrixAutoUpdate === true ) this.updateMatrix();

    if ( this.matrixWorldNeedsUpdate === true || force === true ) {
      if ( this.parent === null ) {
        this.matrixWorld.copy( this.matrix );
      } else {
        this.matrixWorld.multiplyMatrices( this.parent.matrixWorld, this.matrix );
      }

      this.matrixWorldNeedsUpdate = false;

      force = true;
    }

    // update children

    let children = this.children;

    for ( let i = 0, l = children.length; i < l; i ++ ) {
      children[i].updateMatrixWorld( force );
    }
  },
} );


/**
* @author mrdoob / http://mrdoob.com/
*/

function BufferAttribute( array, itemSize, normalized ) {
  if ( Array.isArray( array ) ) {
    throw new TypeError( 'BufferAttribute: array should be a Typed Array.' );
  }

  this.uuid = generateUUID();

  this.array = array;
  this.itemSize = itemSize;
  this.count = array !== undefined ? array.length / itemSize : 0;
  this.normalized = normalized === true;

  this.dynamic = false;
  this.updateRange = { offset: 0, count: - 1 };

  this.onUploadCallback = function () {};

  this.version = 0;
}

BufferAttribute.prototype = {
  constructor: BufferAttribute,
  isBufferAttribute: true,
};


let count = 0;
function GeometryIdCount() { return count++; }

/**
* @author alteredq / http://alteredqualia.com/
* @author mrdoob / http://mrdoob.com/
*/

function BufferGeometry() {
  Object.defineProperty( this, 'id', { value: GeometryIdCount() } );

  this.uuid = generateUUID();

  this.name = '';
  this.type = 'BufferGeometry';

  this.index = null;
  this.attributes = {};

  this.groups = [];

  this.boundingBox = null;
  this.boundingSphere = null;

  this.drawRange = { start: 0, count: Infinity };
}

Object.assign( BufferGeometry.prototype, EventDispatcher.prototype, {

  isBufferGeometry: true,

  setIndex: function ( index ) {
    this.index = index;
  },

  setAttribute: function ( name, attribute ) {
    this.attributes[name] = attribute;
    return this;
  },

  dispose: function () {
    this.dispatchEvent( { type: 'dispose' } );
  },

} );

BufferGeometry.MaxIndex = 65535;

/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
* @author mikael emtinger / http://gomo.se/
* @author jonobr1 / http://jonobr1.com/
*/

function Mesh( geometry, material ) {
  Object3D.call( this );

  this.type = 'Mesh';

  this.geometry = geometry !== undefined ? geometry : new BufferGeometry();
  this.material = material;

  this.drawMode = TrianglesDrawMode;
}

Mesh.prototype = Object.assign( Object.create( Object3D.prototype ), {

  constructor: Mesh,

  isMesh: true,
} );


/**
* @author mrdoob / http://mrdoob.com/
* @author mikael emtinger / http://gomo.se/
* @author WestLangley / http://github.com/WestLangley
*/

function Camera() {
  Object3D.call( this );

  this.type = 'Camera';

  this.matrixWorldInverse = new Matrix4();
  this.projectionMatrix = new Matrix4();
  this.projectionMatrixInverse = new Matrix4();
}

Camera.prototype = Object.create( Object3D.prototype );
Camera.prototype.constructor = Camera;

Camera.prototype.isCamera = true;

Camera.prototype.lookAt = function () {
// This routine does not support cameras with rotated and/or translated parent(s)

let m1 = new Matrix4();

return function lookAt( vector ) {
  m1.lookAt( this.position, vector, this.up );

  this.quaternion.setFromRotationMatrix( m1 );
};
}();

/**
* @author alteredq / http://alteredqualia.com/
* @author arose / http://github.com/arose
*/

function OrthographicCamera( left, right, top, bottom, near, far ) {
  Camera.call( this );

  this.type = 'OrthographicCamera';

  this.zoom = 1;

  this.left = left;
  this.right = right;
  this.top = top;
  this.bottom = bottom;

  this.near = ( near !== undefined ) ? near : 0.1;
  this.far = ( far !== undefined ) ? far : 2000;

  this.updateProjectionMatrix();
}

OrthographicCamera.prototype = Object.assign( Object.create( Camera.prototype ), {

  constructor: OrthographicCamera,

  updateProjectionMatrix: function () {
    let dx = ( this.right - this.left ) / ( 2 * this.zoom );
    let dy = ( this.top - this.bottom ) / ( 2 * this.zoom );
    let cx = ( this.right + this.left ) / 2;
    let cy = ( this.top + this.bottom ) / 2;

    let left = cx - dx;
    let right = cx + dx;
    let top = cy + dy;
    let bottom = cy - dy;

    this.projectionMatrix.makeOrthographic( left, right, top, bottom, this.near, this.far );
    this.projectionMatrixInverse.copy( this.projectionMatrix ).invert();
  },
} );

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLIndexedBufferRenderer( gl, extensions, infoRender ) {
  let mode;

  function setMode( value ) {
    mode = value;
  }

  let type, size;

  function setIndex( index ) {
    if ( index.array instanceof Uint32Array && extensions.get( 'OES_element_index_uint' ) ) {
      type = gl.UNSIGNED_INT;
      size = 4;
    } else if ( index.array instanceof Uint16Array ) {
      type = gl.UNSIGNED_SHORT;
      size = 2;
    } else {
      type = gl.UNSIGNED_BYTE;
      size = 1;
    }
  }

  function render( start, count ) {
    gl.drawElements( mode, count, type, start * size );

    infoRender.calls ++;
    infoRender.vertices += count;

    if ( mode === gl.TRIANGLES ) infoRender.faces += count / 3;
  }

  return {
    setMode: setMode,
    setIndex: setIndex,
    render: render,
  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLBufferRenderer( gl, extensions, infoRender ) {
  let mode;

  function setMode( value ) {
    mode = value;
  }

  function render( start, count ) {
    gl.drawArrays( mode, start, count );

    infoRender.calls ++;
    infoRender.vertices += count;

    if ( mode === gl.TRIANGLES ) infoRender.faces += count / 3;
  }

  return {
    setMode: setMode,
    render: render,
  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/


function WebGLShader( gl, type, string ) {
  let shader = gl.createShader( type );

  gl.shaderSource( shader, string );
  gl.compileShader( shader );

  if ( gl.getShaderParameter( shader, gl.COMPILE_STATUS ) === false ) {
    console.error( 'WebGLShader: Shader couldn\'t compile.' );
  }

  if ( gl.getShaderInfoLog( shader ) !== '' ) {
    let info = gl.getShaderInfoLog( shader );
    // workaround for https://github.com/mrdoob/three.js/issues/9716
    if (info.indexOf('GL_ARB_gpu_shader5') === -1) {
      console.warn( 'WebGLShader: gl.getShaderInfoLog()', type === gl.VERTEX_SHADER ? 'vertex' : 'fragment', info, string );
    }
  }

  return shader;
}

/**
* @author mrdoob / http://mrdoob.com/
*/

let programIdCount = 0;

function generateExtensions( extensions, parameters, rendererExtensions ) {
  extensions = extensions || {};

  let chunks = [
  ( extensions.fragDepth ) && rendererExtensions.get( 'EXT_frag_depth' ) ? '#extension GL_EXT_frag_depth : enable' : '',
  ];

  return chunks.join( '\n' );
}

function fetchAttributeLocations( gl, program ) {
  let attributes = {};

  let n = gl.getProgramParameter( program, gl.ACTIVE_ATTRIBUTES );

  for ( let i = 0; i < n; i ++ ) {
    let info = gl.getActiveAttrib( program, i );
    let name = info.name;

    // console.log("WebGLProgram: ACTIVE VERTEX ATTRIBUTE:", name, i );

    attributes[name] = gl.getAttribLocation( program, name );
  }

  return attributes;
}

function WebGLProgram( renderer, code, material, parameters ) {
  let gl = renderer.context;

  let extensions = material.extensions;

  let vertexShader = material.__webglShader.vertexShader;
  let fragmentShader = material.__webglShader.fragmentShader;

  // console.log( 'building new program ' );

  //

  let customExtensions = generateExtensions( extensions, parameters, renderer.extensions );

  //

  let program = gl.createProgram();

  let prefixVertex, prefixFragment;

  prefixVertex = [
    'precision ' + parameters.precision + ' float;',
    'precision ' + parameters.precision + ' int;',
    '#define SHADER_NAME ' + material.__webglShader.name,
    'uniform mat4 modelMatrix;',
    'uniform mat4 modelViewMatrix;',
    'uniform mat4 projectionMatrix;',
    'uniform mat4 viewMatrix;',
    'attribute vec3 position;',
    '',
  ].join( '\n' );

  prefixFragment = [
    customExtensions,
    'precision ' + parameters.precision + ' float;',
    'precision ' + parameters.precision + ' int;',
    '#define SHADER_NAME ' + material.__webglShader.name,
    ( parameters.useFog && parameters.fog ) ? '#define USE_FOG' : '',
    '',
  ].join( '\n' );

  let vertexGlsl = prefixVertex + vertexShader;
  let fragmentGlsl = prefixFragment + fragmentShader;

  // console.log( '*VERTEX*', vertexGlsl );
  // console.log( '*FRAGMENT*', fragmentGlsl );

  let glVertexShader = WebGLShader( gl, gl.VERTEX_SHADER, vertexGlsl );
  let glFragmentShader = WebGLShader( gl, gl.FRAGMENT_SHADER, fragmentGlsl );

  gl.attachShader( program, glVertexShader );
  gl.attachShader( program, glFragmentShader );

  gl.linkProgram( program );

  let programLog = gl.getProgramInfoLog( program );
  let vertexLog = gl.getShaderInfoLog( glVertexShader );
  let fragmentLog = gl.getShaderInfoLog( glFragmentShader );

  // console.log( '**VERTEX**', gl.getExtension( 'WEBGL_debug_shaders' ).getTranslatedShaderSource( glVertexShader ) );
  // console.log( '**FRAGMENT**', gl.getExtension( 'WEBGL_debug_shaders' ).getTranslatedShaderSource( glFragmentShader ) );

  if ( gl.getProgramParameter( program, gl.LINK_STATUS ) === false ) {
    console.error( 'WebGLProgram: shader error: ', gl.getError(), 'gl.VALIDATE_STATUS', gl.getProgramParameter( program, gl.VALIDATE_STATUS ), 'gl.getProgramInfoLog', programLog, vertexLog, fragmentLog );
  } else if ( programLog !== '' ) {
    console.warn( 'WebGLProgram: gl.getProgramInfoLog()', programLog );
  }

  // clean up

  gl.deleteShader( glVertexShader );
  gl.deleteShader( glFragmentShader );

  // set up caching for uniform locations

  let cachedUniforms;

  this.getUniforms = function () {
    if ( cachedUniforms === undefined ) {
      cachedUniforms =
      new WebGLUniforms( gl, program, renderer );
    }

    return cachedUniforms;
  };

  // set up caching for attribute locations

  let cachedAttributes;

  this.getAttributes = function () {
    if ( cachedAttributes === undefined ) {
      cachedAttributes = fetchAttributeLocations( gl, program );
    }

    return cachedAttributes;
  };

  // free resource

  this.destroy = function () {
    gl.deleteProgram( program );
    this.program = undefined;
  };

  //

  this.id = programIdCount ++;
  this.code = code;
  this.usedTimes = 1;
  this.program = program;
  this.vertexShader = glVertexShader;
  this.fragmentShader = glFragmentShader;

  return this;
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLPrograms( renderer, capabilities ) {
  let programs = [];

  let parameterNames = [
    'precision',
    'fog', 'useFog',
    'premultipliedAlpha',
  ];

  this.getParameters = function ( material, fog ) {
    let precision = renderer.getPrecision();

    if ( material.precision !== null ) {
      precision = capabilities.getMaxPrecision( material.precision );

      if ( precision !== material.precision ) {
        console.warn( 'WebGLProgram.getParameters:', material.precision, 'not supported, using', precision, 'instead.' );
      }
    }

    let parameters = {
      precision: precision,
      fog: !! fog,
      useFog: material.fog,
      premultipliedAlpha: material.premultipliedAlpha,
    };

    return parameters;
  };

  this.getProgramCode = function ( material, parameters ) {
    let array = [];

    array.push( material.fragmentShader );
    array.push( material.vertexShader );

    for ( let i = 0; i < parameterNames.length; i ++ ) {
      array.push( parameters[parameterNames[i]] );
    }

    return array.join();
  };

  this.acquireProgram = function ( material, parameters, code ) {
    let program;

    // Check if code has been already compiled
    for ( let p = 0, pl = programs.length; p < pl; p ++ ) {
      let programInfo = programs[p];

      if ( programInfo.code === code ) {
        program = programInfo;
        ++ program.usedTimes;

        break;
      }
    }

    if ( program === undefined ) {
      program = new WebGLProgram( renderer, code, material, parameters );
      programs.push( program );
    }

    return program;
  };

  this.releaseProgram = function ( program ) {
    if ( -- program.usedTimes === 0 ) {
      // Remove from unordered set
      let i = programs.indexOf( program );
      programs[i] = programs[programs.length - 1];
      programs.pop();

      // Free WebGL resources
      program.destroy();
    }
  };

  // Exposed for resource monitoring & error feedback via renderer.info:
  this.programs = programs;
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLGeometries( gl, properties ) {
  let geometries = {};

  function onGeometryDispose( event ) {
    let geometry = event.target;
    let buffergeometry = geometries[geometry.id];

    if ( buffergeometry.index !== null ) {
      deleteAttribute( buffergeometry.index );
    }

    deleteAttributes( buffergeometry.attributes );

    geometry.removeEventListener( 'dispose', onGeometryDispose );

    delete geometries[geometry.id];

    properties.delete( geometry );

    properties.delete( buffergeometry );
  }

  function getAttributeBuffer( attribute ) {
    return properties.get( attribute ).__webglBuffer;
  }

  function deleteAttribute( attribute ) {
    let buffer = getAttributeBuffer( attribute );

    if ( buffer !== undefined ) {
      gl.deleteBuffer( buffer );
      removeAttributeBuffer( attribute );
    }
  }

  function deleteAttributes( attributes ) {
    for ( let name in attributes ) {
      deleteAttribute( attributes[name] );
    }
  }

  function removeAttributeBuffer( attribute ) {
    properties.delete( attribute );
  }

  return {

    get: function ( object ) {
      let geometry = object.geometry;

      if ( geometries[geometry.id] !== undefined ) {
        return geometries[geometry.id];
      }

      geometry.addEventListener( 'dispose', onGeometryDispose );

      let buffergeometry;

      if ( geometry.isBufferGeometry ) {
        buffergeometry = geometry;
      }

      geometries[geometry.id] = buffergeometry;

      return buffergeometry;
    },

  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLObjects( gl, properties, info ) {
  let geometries = new WebGLGeometries( gl, properties, info );

  //

  function update( object ) {
    let geometry = geometries.get( object );

    let index = geometry.index;
    let attributes = geometry.attributes;

    if ( index !== null ) {
      updateAttribute( index, gl.ELEMENT_ARRAY_BUFFER );
    }

    for ( let name in attributes ) {
      updateAttribute( attributes[name], gl.ARRAY_BUFFER );
    }

    return geometry;
  }

  function updateAttribute( attribute, bufferType ) {
    let data = attribute;

    let attributeProperties = properties.get( data );

    if ( attributeProperties.__webglBuffer === undefined ) {
      createBuffer( attributeProperties, data, bufferType );
    } else if ( attributeProperties.version !== data.version ) {
      updateBuffer( attributeProperties, data, bufferType );
    }
  }

  function createBuffer( attributeProperties, data, bufferType ) {
    attributeProperties.__webglBuffer = gl.createBuffer();
    gl.bindBuffer( bufferType, attributeProperties.__webglBuffer );

    let usage = data.dynamic ? gl.DYNAMIC_DRAW : gl.STATIC_DRAW;

    gl.bufferData( bufferType, data.array, usage );

    let type = gl.FLOAT;
    let array = data.array;

    if ( array instanceof Float32Array ) {
      type = gl.FLOAT;
    } else if ( array instanceof Float64Array ) {
      console.warn( 'Unsupported data buffer format: Float64Array' );
    } else if ( array instanceof Uint16Array ) {
      type = gl.UNSIGNED_SHORT;
    } else if ( array instanceof Int16Array ) {
      type = gl.SHORT;
    } else if ( array instanceof Uint32Array ) {
      type = gl.UNSIGNED_INT;
    } else if ( array instanceof Int32Array ) {
      type = gl.INT;
    } else if ( array instanceof Int8Array ) {
      type = gl.BYTE;
    } else if ( array instanceof Uint8Array ) {
      type = gl.UNSIGNED_BYTE;
    }

    attributeProperties.bytesPerElement = array.BYTES_PER_ELEMENT;
    attributeProperties.type = type;
    attributeProperties.version = data.version;

    data.onUploadCallback();
  }

  function updateBuffer( attributeProperties, data, bufferType ) {
    gl.bindBuffer( bufferType, attributeProperties.__webglBuffer );

    if ( data.dynamic === false ) {
      gl.bufferData( bufferType, data.array, gl.STATIC_DRAW );
    } else if ( data.updateRange.count === - 1 ) {
      // Not using update ranges

      gl.bufferSubData( bufferType, 0, data.array );
    } else if ( data.updateRange.count === 0 ) {
      console.error( 'WebGLObjects.updateBuffer: updateRange.count is 0.' );
    } else {
      gl.bufferSubData( bufferType, data.updateRange.offset * data.array.BYTES_PER_ELEMENT,
                        data.array.subarray( data.updateRange.offset, data.updateRange.offset + data.updateRange.count ) );

      data.updateRange.count = 0; // reset range
    }

    attributeProperties.version = data.version;
  }

  function getAttributeBuffer( attribute ) {
    return properties.get( attribute ).__webglBuffer;
  }

  function getAttributeProperties( attribute ) {
    return properties.get( attribute );
  }


  return {

    getAttributeBuffer: getAttributeBuffer,
    getAttributeProperties: getAttributeProperties,

    update: update,

  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLTextures( _gl, extensions, state, properties ) {
  function onTextureDispose( event ) {
    let texture = event.target;

    texture.removeEventListener( 'dispose', onTextureDispose );

    deallocateTexture( texture );
  }

  //

  function deallocateTexture( texture ) {
    let textureProperties = properties.get( texture );

    // 2D texture

    if ( textureProperties.__webglInit === undefined ) return;

    _gl.deleteTexture( textureProperties.__webglTexture );

    // remove all webgl properties
    properties.delete( texture );
  }


  function setTexture2D( texture, slot ) {
    let textureProperties = properties.get( texture );

    if ( texture.version > 0 && textureProperties.__version !== texture.version ) {
      let image = texture.image;

      if ( image === undefined ) {
        console.warn( 'WebGLRenderer: Texture marked for update but image is undefined', texture );
      } else if ( image.complete === false ) {
        console.warn( 'WebGLRenderer: Texture marked for update but image is incomplete', texture );
      } else {
        uploadTexture( textureProperties, texture, slot );
        return;
      }
    }

    state.activeTexture( _gl.TEXTURE0 + slot );
    state.bindTexture( _gl.TEXTURE_2D, textureProperties.__webglTexture );
  }

  function setTextureParameters( textureType ) {
    _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_S, _gl.CLAMP_TO_EDGE );
    _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_T, _gl.CLAMP_TO_EDGE );
    _gl.texParameteri( textureType, _gl.TEXTURE_MAG_FILTER, _gl.LINEAR );
    _gl.texParameteri( textureType, _gl.TEXTURE_MIN_FILTER, _gl.LINEAR_MIPMAP_LINEAR );
  }

  function uploadTexture( textureProperties, texture, slot ) {
    if ( textureProperties.__webglInit === undefined ) {
      textureProperties.__webglInit = true;

      texture.addEventListener( 'dispose', onTextureDispose );

      textureProperties.__webglTexture = _gl.createTexture();
    }

    state.activeTexture( _gl.TEXTURE0 + slot );
    state.bindTexture( _gl.TEXTURE_2D, textureProperties.__webglTexture );

    _gl.pixelStorei( _gl.UNPACK_FLIP_Y_WEBGL, true );
    _gl.pixelStorei( _gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, false );
    _gl.pixelStorei( _gl.UNPACK_ALIGNMENT, 4 );

    let image = texture.image;

    let glFormat = _gl.RGBA;
    let glType = _gl.UNSIGNED_BYTE;

    setTextureParameters( _gl.TEXTURE_2D );

    state.texImage2D( _gl.TEXTURE_2D, 0, glFormat, glFormat, glType, image );

    _gl.generateMipmap( _gl.TEXTURE_2D );

    textureProperties.__version = texture.version;
  }

  this.setTexture2D = setTexture2D;
}

/**
* @author fordacious / fordacious.github.io
*/

function WebGLProperties() {
  let properties = {};

  return {

    get: function ( object ) {
      let uuid = object.uuid;
      let map = properties[uuid];

      if ( map === undefined ) {
        map = {};
        properties[uuid] = map;
      }

      return map;
    },

    delete: function ( object ) {
      delete properties[object.uuid];
    },

  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLState( gl ) {
  function ColorBuffer() {
    let color = new Vector4();
    let currentColorClear = new Vector4();

    return {

      setClear: function ( r, g, b, a, premultipliedAlpha ) {
        if ( premultipliedAlpha === true ) {
          r *= a; g *= a; b *= a;
        }

        color.set( r, g, b, a );

        if ( currentColorClear.equals( color ) === false ) {
          gl.clearColor( r, g, b, a );
          currentColorClear.copy( color );
        }
      },

    };
  }

  function DepthBuffer() {
    let currentDepthMask = null;
    let currentDepthClear = null;

    return {

      setTest: function ( depthTest ) {
        if ( depthTest ) {
          enable( gl.DEPTH_TEST );
        } else {
          disable( gl.DEPTH_TEST );
        }
      },

      setMask: function ( depthMask ) {
        if ( currentDepthMask !== depthMask ) {
          gl.depthMask( depthMask );
          currentDepthMask = depthMask;
        }
      },

      setClear: function ( depth ) {
        if ( currentDepthClear !== depth ) {
          gl.clearDepth( depth );
          currentDepthClear = depth;
        }
      },

      reset: function () {
        currentDepthMask = null;
        currentDepthClear = null;
      },

    };
  }


  //

  let colorBuffer = new ColorBuffer();
  let depthBuffer = new DepthBuffer();

  let maxVertexAttributes = gl.getParameter( gl.MAX_VERTEX_ATTRIBS );
  let newAttributes = new Uint8Array( maxVertexAttributes );
  let enabledAttributes = new Uint8Array( maxVertexAttributes );

  let capabilities = {};

  let currentBlending = null;
  let currentPremultipledAlpha = false;

  let currentLineWidth = null;

  let version = parseFloat( /^WebGL\ ([0-9])/.exec( gl.getParameter( gl.VERSION ) )[1] );
  let lineWidthAvailable = parseFloat( version ) >= 1.0;

  let currentTextureSlot = null;
  let currentBoundTextures = {};

  let currentViewport = new Vector4();

  function createTexture( type, target, count ) {
    let data = new Uint8Array( 4 ); // 4 is required to match default unpack alignment of 4.
    let texture = gl.createTexture();

    gl.bindTexture( type, texture );
    gl.texParameteri( type, gl.TEXTURE_MIN_FILTER, gl.NEAREST );
    gl.texParameteri( type, gl.TEXTURE_MAG_FILTER, gl.NEAREST );

    for ( let i = 0; i < count; i ++ ) {
      gl.texImage2D( target + i, 0, gl.RGBA, 1, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, data );
    }

    return texture;
  }

  let emptyTextures = {};
  emptyTextures[gl.TEXTURE_2D] = createTexture( gl.TEXTURE_2D, gl.TEXTURE_2D, 1 );

  //

  function init() {
    colorBuffer.setClear( 0, 0, 0, 1 );
    depthBuffer.setClear( 1 );

    enable( gl.DEPTH_TEST );
    gl.depthFunc( gl.LEQUAL );

    enable( gl.BLEND );
    setBlending( NormalBlending );
  }

  function initAttributes() {
    for ( let i = 0, l = newAttributes.length; i < l; i ++ ) {
      newAttributes[i] = 0;
    }
  }

  function enableAttribute( attribute ) {
    newAttributes[attribute] = 1;

    if ( enabledAttributes[attribute] === 0 ) {
      gl.enableVertexAttribArray( attribute );
      enabledAttributes[attribute] = 1;
    }
  }

  function disableUnusedAttributes() {
    for ( let i = 0, l = enabledAttributes.length; i !== l; ++ i ) {
      if ( enabledAttributes[i] !== newAttributes[i] ) {
        gl.disableVertexAttribArray( i );
        enabledAttributes[i] = 0;
      }
    }
  }

  function enable( id ) {
    if ( capabilities[id] !== true ) {
      gl.enable( id );
      capabilities[id] = true;
    }
  }

  function disable( id ) {
    if ( capabilities[id] !== false ) {
      gl.disable( id );
      capabilities[id] = false;
    }
  }

  function setBlending( blending, premultipliedAlpha ) {
    if ( blending !== NoBlending ) {
      enable( gl.BLEND );
    } else {
      disable( gl.BLEND );
    }

    if ( blending !== currentBlending || premultipliedAlpha !== currentPremultipledAlpha ) {
      if ( premultipliedAlpha ) {
        gl.blendEquationSeparate( gl.FUNC_ADD, gl.FUNC_ADD );
        gl.blendFuncSeparate( gl.ONE, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA );
      } else {
        gl.blendEquationSeparate( gl.FUNC_ADD, gl.FUNC_ADD );
        gl.blendFuncSeparate( gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA );
      }

      currentBlending = blending;
      currentPremultipledAlpha = premultipliedAlpha;
    }
  }

  function setDepthTest( depthTest ) {
    depthBuffer.setTest( depthTest );
  }

  function setDepthWrite( depthWrite ) {
    depthBuffer.setMask( depthWrite );
  }

  //

  function setLineWidth( width ) {
    if ( width !== currentLineWidth ) {
      if ( lineWidthAvailable ) gl.lineWidth( width );

      currentLineWidth = width;
    }
  }

  // texture

  function activeTexture( webglSlot ) {
    if ( currentTextureSlot !== webglSlot ) {
      gl.activeTexture( webglSlot );
      currentTextureSlot = webglSlot;
    }
  }

  function bindTexture( webglType, webglTexture ) {
    let boundTexture = currentBoundTextures[currentTextureSlot];
    if ( boundTexture === undefined ) {
      boundTexture = { type: undefined, texture: undefined };
      currentBoundTextures[currentTextureSlot] = boundTexture;
    }

    if ( boundTexture.type !== webglType || boundTexture.texture !== webglTexture ) {
      gl.bindTexture( webglType, webglTexture || emptyTextures[webglType] );

      boundTexture.type = webglType;
      boundTexture.texture = webglTexture;
    }
  }

  function texImage2D() {
    try {
      gl.texImage2D.apply( gl, arguments );
    } catch ( error ) {
      console.error( error );
    }
  }

  //

  function viewport( viewport ) {
    if ( currentViewport.equals( viewport ) === false ) {
      gl.viewport( viewport.x, viewport.y, viewport.z, viewport.w );
      currentViewport.copy( viewport );
    }
  }

  //

  return {
    buffers: {
      color: colorBuffer,
      depth: depthBuffer,
    },

    init: init,
    initAttributes: initAttributes,
    enableAttribute: enableAttribute,
    disableUnusedAttributes: disableUnusedAttributes,
    enable: enable,
    disable: disable,

    setBlending: setBlending,

    setDepthTest: setDepthTest,
    setDepthWrite: setDepthWrite,

    setLineWidth: setLineWidth,

    activeTexture: activeTexture,
    bindTexture: bindTexture,
    texImage2D: texImage2D,

    viewport: viewport,
  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLCapabilities( gl, extensions, parameters ) {
  function getMaxPrecision( precision ) {
    if ( precision === 'highp' ) {
      if ( gl.getShaderPrecisionFormat( gl.VERTEX_SHADER, gl.HIGH_FLOAT ).precision > 0 &&
         gl.getShaderPrecisionFormat( gl.FRAGMENT_SHADER, gl.HIGH_FLOAT ).precision > 0 ) {
        return 'highp';
      }

      precision = 'mediump';
    }

    if ( precision === 'mediump' ) {
      if ( gl.getShaderPrecisionFormat( gl.VERTEX_SHADER, gl.MEDIUM_FLOAT ).precision > 0 &&
         gl.getShaderPrecisionFormat( gl.FRAGMENT_SHADER, gl.MEDIUM_FLOAT ).precision > 0 ) {
        return 'mediump';
      }
    }

    return 'lowp';
  }

  let precision = parameters.precision !== undefined ? parameters.precision : 'highp';
  let maxPrecision = getMaxPrecision( precision );

  if ( maxPrecision !== precision ) {
    console.warn( 'WebGLRenderer:', precision, 'not supported, using', maxPrecision, 'instead.' );
    precision = maxPrecision;
  }

  let maxTextures = gl.getParameter( gl.MAX_TEXTURE_IMAGE_UNITS );

  return {
    getMaxPrecision: getMaxPrecision,
    precision: precision,
    maxTextures: maxTextures,
  };
}

/**
* @author mrdoob / http://mrdoob.com/
*/

function WebGLExtensions( gl ) {
  let extensions = {};

  return {

    get: function ( name ) {
      if ( extensions[name] !== undefined ) {
        return extensions[name];
      }

      let extension;

      switch ( name ) {
        case 'WEBGL_depth_texture':
          extension = gl.getExtension( 'WEBGL_depth_texture' ) || gl.getExtension( 'MOZ_WEBGL_depth_texture' ) || gl.getExtension( 'WEBKIT_WEBGL_depth_texture' );
          break;
        default:
          extension = gl.getExtension( name );
      }

      if ( extension === null ) {
        console.warn( 'WebGLRenderer: ' + name + ' extension not supported.' );
      }

      extensions[name] = extension;

      return extension;
    },

  };
}

/**
* @author supereggbert / http://www.paulbrunt.co.uk/
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
* @author szimek / https://github.com/szimek/
* @author tschw
*/

function WebGLRenderer( parameters ) {
  parameters = parameters || {};

  let _canvas = parameters.canvas !== undefined ? parameters.canvas : document.createElementNS( 'http://www.w3.org/1999/xhtml', 'canvas' ),
    _context = parameters.context !== undefined ? parameters.context : null,
    _antialias = parameters.antialias !== undefined ? parameters.antialias : false,
    _premultipliedAlpha = parameters.premultipliedAlpha !== undefined ? parameters.premultipliedAlpha : true;

  let opaqueObjects = [];
  let opaqueObjectsLastIndex = - 1;
  let transparentObjects = [];
  let transparentObjectsLastIndex = - 1;

  // public properties

  this.domElement = _canvas;
  this.context = null;

  // clearing

  this.autoClear = true;
  this.autoClearColor = true;
  this.autoClearDepth = true;

  // scene graph

  this.sortObjects = true;

  // internal properties

  let _this = this,

    // internal state cache

    _currentProgram = null,
    _currentRenderTarget = null,
    _currentFramebuffer = null,
    _currentMaterialId = - 1,
    _currentGeometryProgram = '',
    _currentCamera = null,

    _currentViewport = new Vector4(),

    //

    _usedTextureUnits = 0,

    //

    _clearColor = new Color( 0x000000 ),
    _clearAlpha = 0,

    _width = _canvas.width,
    _height = _canvas.height,

    _pixelRatio = 1,

    _viewport = new Vector4( 0, 0, _width, _height ),

    // camera matrices cache

    _projScreenMatrix = new Matrix4(),

    _vector3 = new Vector3(),

    // info

    _infoRender = {
      calls: 0,
      vertices: 0,
      faces: 0,
      points: 0,
    };

  this.info = {
    render: _infoRender,
    programs: null,
  };


  // initialize

  let _gl;

  try {
    let attributes = {
      alpha: false,
      depth: true,
      antialias: _antialias,
      premultipliedAlpha: _premultipliedAlpha,
      preserveDrawingBuffer: false,
    };

    _gl = _context || _canvas.getContext( 'webgl', attributes ) || _canvas.getContext( 'experimental-webgl', attributes );

    if ( _gl === null ) {
      if ( _canvas.getContext( 'webgl' ) !== null ) {
        throw new Error('Error creating WebGL context with your selected attributes.');
      } else {
        throw new Error('Error creating WebGL context.');
      }
    }

    // Some experimental-webgl implementations do not have getShaderPrecisionFormat

    if ( _gl.getShaderPrecisionFormat === undefined ) {
      _gl.getShaderPrecisionFormat = function () {
        return { 'rangeMin': 1, 'rangeMax': 1, 'precision': 1 };
      };
    }

    _canvas.addEventListener( 'webglcontextlost', onContextLost, false );
  } catch ( error ) {
    console.error( 'WebGLRenderer: ' + error );
  }

  let extensions = new WebGLExtensions( _gl );
  extensions.get( 'WEBGL_depth_texture' );
  if ( extensions.get( 'OES_element_index_uint' ) ) {
    BufferGeometry.MaxIndex = 4294967296;
  }

  let capabilities = new WebGLCapabilities( _gl, extensions, parameters );

  let state = new WebGLState( _gl, extensions );
  let properties = new WebGLProperties();
  let textures = new WebGLTextures( _gl, extensions, state, properties, capabilities, this.info );
  let objects = new WebGLObjects( _gl, properties, this.info );
  let programCache = new WebGLPrograms( this, capabilities );

  this.info.programs = programCache.programs;

  let bufferRenderer = new WebGLBufferRenderer( _gl, extensions, _infoRender );
  let indexedBufferRenderer = new WebGLIndexedBufferRenderer( _gl, extensions, _infoRender );

  //

  function getTargetPixelRatio() {
    return _currentRenderTarget === null ? _pixelRatio : 1;
  }

  function setDefaultGLState() {
    state.init();

    state.viewport( _currentViewport.copy( _viewport ).multiplyScalar( _pixelRatio ) );

    state.buffers.color.setClear( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha, _premultipliedAlpha );
  }

  function resetGLState() {
    _currentProgram = null;
    _currentCamera = null;

    _currentGeometryProgram = '';
    _currentMaterialId = - 1;

    state.reset();
  }

  setDefaultGLState();

  this.context = _gl;
  this.capabilities = capabilities;
  this.extensions = extensions;
  this.properties = properties;
  this.state = state;

  // API

  this.getContext = function () {
    return _gl;
  };

  this.getPrecision = function () {
    return capabilities.precision;
  };

  this.getPixelRatio = function () {
    return _pixelRatio;
  };

  this.setPixelRatio = function ( value ) {
    if ( value === undefined ) return;

    _pixelRatio = value;

    this.setSize( _viewport.z, _viewport.w, false );
  };

  this.setSize = function ( width, height, updateStyle ) {
    _width = width;
    _height = height;

    _canvas.width = width * _pixelRatio;
    _canvas.height = height * _pixelRatio;

    if ( updateStyle !== false ) {
      _canvas.style.width = width + 'px';
      _canvas.style.height = height + 'px';
    }

    this.setViewport( 0, 0, width, height );
  };

  this.setViewport = function ( x, y, width, height ) {
    state.viewport( _viewport.set( x, y, width, height ) );
  };

  // Clearing

  this.getClearColor = function () {
    return _clearColor;
  };

  this.setClearColor = function ( color, alpha ) {
    _clearColor.set( color );

    _clearAlpha = alpha !== undefined ? alpha : 1;

    state.buffers.color.setClear( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha, _premultipliedAlpha );
  };

  this.clear = function ( color, depth ) {
    let bits = 0;

    if ( color === undefined || color ) bits |= _gl.COLOR_BUFFER_BIT;
    if ( depth === undefined || depth ) bits |= _gl.DEPTH_BUFFER_BIT;

    _gl.clear( bits );
  };

  this.clearColor = function () {
    this.clear( true, false );
  };

  this.clearDepth = function () {
    this.clear( false, true );
  };

  // Reset

  this.resetGLState = resetGLState;

  this.dispose = function () {
    transparentObjects = [];
    transparentObjectsLastIndex = -1;
    opaqueObjects = [];
    opaqueObjectsLastIndex = -1;

    _canvas.removeEventListener( 'webglcontextlost', onContextLost, false );
  };

  // Events

  function onContextLost( event ) {
    event.preventDefault();

    resetGLState();
    setDefaultGLState();

    properties.clear();
  }

  function onMaterialDispose( event ) {
    let material = event.target;

    material.removeEventListener( 'dispose', onMaterialDispose );

    deallocateMaterial( material );
  }

  // Buffer deallocation

  function deallocateMaterial( material ) {
    releaseMaterialProgramReference( material );

    properties.delete( material );
  }


  function releaseMaterialProgramReference( material ) {
    let programInfo = properties.get( material ).program;

    material.program = undefined;

    if ( programInfo !== undefined ) {
      programCache.releaseProgram( programInfo );
    }
  }

  // Buffer rendering


  this.renderBufferDirect = function ( camera, fog, geometry, material, object, group ) {
    setMaterial( material );

    let program = setProgram( camera, fog, material, object );

    let updateBuffers = false;
    let geometryProgram = geometry.id + '_' + program.id;

    if ( geometryProgram !== _currentGeometryProgram ) {
      _currentGeometryProgram = geometryProgram;
      updateBuffers = true;
    }

    //

    let index = geometry.index;
    let position = geometry.attributes.position;
    let rangeFactor = 1;

    let renderer;

    if ( index !== null ) {
      renderer = indexedBufferRenderer;
      renderer.setIndex( index );
    } else {
      renderer = bufferRenderer;
    }

    if ( updateBuffers ) {
      setupVertexAttributes( material, program, geometry );

      if ( index !== null ) {
        _gl.bindBuffer( _gl.ELEMENT_ARRAY_BUFFER, objects.getAttributeBuffer( index ) );
      }
    }

    //

    let dataCount = 0;

    if ( index !== null ) {
      dataCount = index.count;
    } else if ( position !== undefined ) {
      dataCount = position.count;
    }

    let rangeStart = geometry.drawRange.start * rangeFactor;
    let rangeCount = geometry.drawRange.count * rangeFactor;

    let groupStart = group !== null ? group.start * rangeFactor : 0;
    let groupCount = group !== null ? group.count * rangeFactor : Infinity;

    let drawStart = Math.max( rangeStart, groupStart );
    let drawEnd = Math.min( dataCount, rangeStart + rangeCount, groupStart + groupCount ) - 1;

    let drawCount = Math.max( 0, drawEnd - drawStart + 1 );

    if ( drawCount === 0 ) return;

    //

    if ( object.isMesh ) {
      switch ( object.drawMode ) {
        case TrianglesDrawMode:
          renderer.setMode( _gl.TRIANGLES );
          break;

        case TriangleStripDrawMode:
          renderer.setMode( _gl.TRIANGLE_STRIP );
          break;

        case TriangleFanDrawMode:
          renderer.setMode( _gl.TRIANGLE_FAN );
          break;
      }
    } else if ( object.isLine ) {
      let lineWidth = material.linewidth;
      if ( lineWidth === undefined ) lineWidth = 1; // Not using Line*Material
      state.setLineWidth( lineWidth * getTargetPixelRatio() );

      if ( object.isLineSegments ) {
        renderer.setMode( _gl.LINES );
      } else {
        renderer.setMode( _gl.LINE_STRIP );
      }
    } else if ( object.isPoints ) {
      renderer.setMode( _gl.POINTS );
    }

    renderer.render( drawStart, drawCount );
  };

  function setupVertexAttributes( material, program, geometry, startIndex ) {
    if ( startIndex === undefined ) startIndex = 0;

    state.initAttributes();

    let geometryAttributes = geometry.attributes;

    let programAttributes = program.getAttributes();

    for ( let name in programAttributes ) {
      let programAttribute = programAttributes[name];

      if ( programAttribute >= 0 ) {
        let geometryAttribute = geometryAttributes[name];

        if ( geometryAttribute !== undefined ) {
          let normalized = geometryAttribute.normalized;
          let size = geometryAttribute.itemSize;

          let attributeProperties = objects.getAttributeProperties( geometryAttribute );

          let buffer = attributeProperties.__webglBuffer;
          let type = attributeProperties.type;
          let bytesPerElement = attributeProperties.bytesPerElement;

          state.enableAttribute( programAttribute );

          _gl.bindBuffer( _gl.ARRAY_BUFFER, buffer );
          _gl.vertexAttribPointer( programAttribute, size, type, normalized, 0, startIndex * size * bytesPerElement );
        } else {
          console.error( 'undefined geometryAttribute' );
        }
      }
    }
    state.disableUnusedAttributes();
  }

  // Sorting

  function painterSortStable( a, b ) {
    if ( a.object.renderOrder !== b.object.renderOrder ) {
      return a.object.renderOrder - b.object.renderOrder;
    } else if ( a.material.program && b.material.program && a.material.program !== b.material.program ) {
      return a.material.program.id - b.material.program.id;
    } else if ( a.material.id !== b.material.id ) {
      return a.material.id - b.material.id;
    } else if ( a.z !== b.z ) {
      return a.z - b.z;
    } else {
      return a.id - b.id;
    }
  }

  function reversePainterSortStable( a, b ) {
    if ( a.object.renderOrder !== b.object.renderOrder ) {
      return a.object.renderOrder - b.object.renderOrder;
    } if ( a.z !== b.z ) {
      return b.z - a.z;
    } else {
      return a.id - b.id;
    }
  }

  // Rendering

  this.render = function ( scene, camera, renderTarget, forceClear ) {
    if ( camera !== undefined && camera.isCamera !== true ) {
      console.error( 'camera is not an instance of Camera.' );
      return;
    }

    // reset caching for this frame

    _currentGeometryProgram = '';
    _currentMaterialId = - 1;
    _currentCamera = null;

    // update scene graph

    if ( scene.autoUpdate === true ) scene.updateMatrixWorld();

    // update camera matrices and frustum

    if ( camera.parent === null ) camera.updateMatrixWorld();

    camera.matrixWorldInverse.copy( camera.matrixWorld ).invert();

    _projScreenMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse );

    opaqueObjectsLastIndex = - 1;
    transparentObjectsLastIndex = - 1;

    projectObject( scene );

    opaqueObjects.length = opaqueObjectsLastIndex + 1;
    transparentObjects.length = transparentObjectsLastIndex + 1;

    if ( _this.sortObjects === true ) {
      opaqueObjects.sort( painterSortStable );
      transparentObjects.sort( reversePainterSortStable );
    }

    //

    _infoRender.calls = 0;
    _infoRender.vertices = 0;
    _infoRender.faces = 0;
    _infoRender.points = 0;

    if ( renderTarget === undefined ) {
      renderTarget = null;
    }

    this.setRenderTarget( renderTarget );

    state.buffers.color.setClear( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha, _premultipliedAlpha );

    if ( this.autoClear || forceClear ) {
      this.clear( this.autoClearColor, this.autoClearDepth );
    }

    // opaque pass (front-to-back order)

    state.setBlending( NoBlending );
    renderObjects( opaqueObjects, scene, camera );

    // transparent pass (back-to-front order)

    renderObjects( transparentObjects, scene, camera );

    // Generate mipmap if we're using any kind of mipmap filtering

    if ( renderTarget ) {
      textures.updateRenderTargetMipmap( renderTarget );
    }

    // Ensure depth buffer writing is enabled so it can be cleared on next render

    state.setDepthTest( true );
    state.setDepthWrite( true );

    // _gl.finish();
  };

  function pushRenderItem( object, geometry, material, z, group ) {
    let array, index;

    // allocate the next position in the appropriate array

    if ( material.transparent ) {
      array = transparentObjects;
      index = ++ transparentObjectsLastIndex;
    } else {
      array = opaqueObjects;
      index = ++ opaqueObjectsLastIndex;
    }

    // recycle existing render item or grow the array

    let renderItem = array[index];

    if ( renderItem !== undefined ) {
      renderItem.id = object.id;
      renderItem.object = object;
      renderItem.geometry = geometry;
      renderItem.material = material;
      renderItem.z = _vector3.z;
      renderItem.group = group;
    } else {
      renderItem = {
        id: object.id,
        object: object,
        geometry: geometry,
        material: material,
        z: _vector3.z,
        group: group,
      };

      // assert( index === array.length );
      array.push( renderItem );
    }
  }

  function projectObject( object ) {
    if ( object.visible === false ) return;

    if ( object.isMesh || object.isLine || object.isPoints ) {
      let material = object.material;

      if ( material.visible === true ) {
        if ( _this.sortObjects === true ) {
          _vector3.setFromMatrixPosition( object.matrixWorld );
          _vector3.applyMatrix4( _projScreenMatrix );
        }

        let geometry = objects.update( object );

        pushRenderItem( object, geometry, material, _vector3.z, null );
      }
    }

    let children = object.children;

    for ( let i = 0, l = children.length; i < l; i ++ ) {
      projectObject( children[i] );
    }
  }

  function renderObjects( renderList, scene, camera, overrideMaterial ) {
    for ( let i = 0, l = renderList.length; i < l; i ++ ) {
      let renderItem = renderList[i];
      let object = renderItem.object;
      let geometry = renderItem.geometry;
      let material = overrideMaterial === undefined ? renderItem.material : overrideMaterial;
      let group = renderItem.group;

      object.modelViewMatrix.multiplyMatrices( camera.matrixWorldInverse, object.matrixWorld );
      _this.renderBufferDirect( camera, scene.fog, geometry, material, object, group );
    }
  }

  function initMaterial( material, fog, object ) {
    let materialProperties = properties.get( material );

    let parameters = programCache.getParameters(
      material, fog, object );

    let code = programCache.getProgramCode( material, parameters );

    let program = materialProperties.program;
    let programChange = true;

    if ( program === undefined ) {
      // new material
      material.addEventListener( 'dispose', onMaterialDispose );
    } else if ( program.code !== code ) {
      // changed glsl or parameters
      releaseMaterialProgramReference( material );
    } else {
      // only rebuild uniform list
      programChange = false;
    }

    if ( programChange ) {
      materialProperties.__webglShader = {
        name: material.type,
        uniforms: material.uniforms,
        vertexShader: material.vertexShader,
        fragmentShader: material.fragmentShader,
      };

      material.__webglShader = materialProperties.__webglShader;

      program = programCache.acquireProgram( material, parameters, code );

      materialProperties.program = program;
      material.program = program;
    }

    let uniforms = materialProperties.__webglShader.uniforms;

    materialProperties.fog = fog;

    let progUniforms = materialProperties.program.getUniforms(),
      uniformsList =
      WebGLUniforms.seqWithValue( progUniforms.seq, uniforms );

    materialProperties.uniformsList = uniformsList;
  }

  function setMaterial( material ) {
    material.transparent === true ?
      state.setBlending( NormalBlending, material.premultipliedAlpha )
      : state.setBlending( NoBlending );

    state.setDepthTest( material.depthTest );
    state.setDepthWrite( material.depthWrite );
  }

  function setProgram( camera, fog, material, object ) {
    _usedTextureUnits = 0;

    let materialProperties = properties.get( material );

    if ( material.needsUpdate === false ) {
      if ( materialProperties.program === undefined ) {
        material.needsUpdate = true;
      } else if ( material.fog && materialProperties.fog !== fog ) {
        material.needsUpdate = true;
      }
    }

    if ( material.needsUpdate ) {
      initMaterial( material, fog, object );
      material.needsUpdate = false;
    }

    let refreshProgram = false;
    let refreshMaterial = false;

    let program = materialProperties.program,
      p_uniforms = program.getUniforms(),
      m_uniforms = materialProperties.__webglShader.uniforms;

    if ( program.id !== _currentProgram ) {
      _gl.useProgram( program.program );
      _currentProgram = program.id;

      refreshProgram = true;
      refreshMaterial = true;
    }

    if ( material.id !== _currentMaterialId ) {
      _currentMaterialId = material.id;

      refreshMaterial = true;
    }

    if ( refreshProgram || camera !== _currentCamera ) {
      p_uniforms.set( _gl, camera, 'projectionMatrix' );

      if ( camera !== _currentCamera ) {
        _currentCamera = camera;

        // lighting uniforms depend on the camera so enforce an update
        // now, in case this material supports lights - or later, when
        // the next material that does gets activated:

        refreshMaterial = true;   // set to true on material change
      }

      // load material specific uniforms
      // (shader material also gets them for the sake of genericity)

      if ( material.isShaderMaterial ) {
        p_uniforms.setValue( _gl, 'viewMatrix', camera.matrixWorldInverse );
      }
    }

    if ( refreshMaterial ) {
      // refresh uniforms common to several materials

      if ( fog && material.fog ) {
        refreshUniformsFog( m_uniforms, fog );
      }

      // refresh single material specific uniforms

      WebGLUniforms.upload(
        _gl, materialProperties.uniformsList, m_uniforms, _this );
    }


    // common matrices

    p_uniforms.set( _gl, object, 'modelViewMatrix' );
    p_uniforms.setValue( _gl, 'modelMatrix', object.matrixWorld );

    return program;
  }

  // Uniforms (refresh uniforms objects)

  function refreshUniformsFog( uniforms, fog ) {
    uniforms.fogColor.value = fog.color;

    if ( fog.isFog ) {
      uniforms.fogNear.value = fog.near;
      uniforms.fogFar.value = fog.far;
    }
  }


  // Textures

  function allocTextureUnit() {
    let textureUnit = _usedTextureUnits;
    _usedTextureUnits += 1;
    return textureUnit;
  }

  this.allocTextureUnit = allocTextureUnit;

  // this.setTexture2D = setTexture2D;
  this.setTexture2D = ( function () {
    let warned = false;

    // backwards compatibility: peel texture.texture
    return function setTexture2D( texture, slot ) {
      if ( texture && texture.isWebGLRenderTarget ) {
        if ( ! warned ) {
          console.warn( 'WebGLRenderer.setTexture2D: don\'t use render targets as textures. Use their .texture property instead.' );
          warned = true;
        }

        texture = texture.texture;
      }

      textures.setTexture2D( texture, slot );
    };
  }() );

  this.setRenderTarget = function ( renderTarget ) {
    _currentRenderTarget = renderTarget;

    if ( renderTarget && properties.get( renderTarget ).__webglFramebuffer === undefined ) {
      textures.setupRenderTarget( renderTarget );
    }

    let framebuffer;

    if ( renderTarget ) {
      let renderTargetProperties = properties.get( renderTarget );

      framebuffer = renderTargetProperties.__webglFramebuffer;

      _currentViewport.copy( renderTarget.viewport );
    } else {
      framebuffer = null;

      _currentViewport.copy( _viewport ).multiplyScalar( _pixelRatio );
    }

    if ( _currentFramebuffer !== framebuffer ) {
      _gl.bindFramebuffer( _gl.FRAMEBUFFER, framebuffer );
      _currentFramebuffer = framebuffer;
    }

    state.viewport( _currentViewport );
  };
}

/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
*/

function Fog( color, near, far ) {
  this.name = '';

  this.color = new Color( color );

  this.near = ( near !== undefined ) ? near : 1;
  this.far = ( far !== undefined ) ? far : 1000;
}

Fog.prototype.isFog = true;


/**
* @author mrdoob / http://mrdoob.com/
*/

function Scene() {
  Object3D.call( this );

  this.type = 'Scene';

  this.fog = null;

  this.autoUpdate = true; // checked by the renderer
}

Scene.prototype = Object.create( Object3D.prototype );

Scene.prototype.constructor = Scene;


/**
* @author mrdoob / http://mrdoob.com/
*/

function Line( geometry, material ) {
  Object3D.call( this );

  this.type = 'Line';

  this.geometry = geometry !== undefined ? geometry : new BufferGeometry();
  this.material = material;
}

Line.prototype = Object.assign( Object.create( Object3D.prototype ), {

  constructor: Line,

  isLine: true,

} );

/**
* @author mrdoob / http://mrdoob.com/
*/

function LineSegments( geometry, material ) {
  Line.call( this, geometry, material );

  this.type = 'LineSegments';
}

LineSegments.prototype = Object.assign( Object.create( Line.prototype ), {

  constructor: LineSegments,

  isLineSegments: true,

} );


/**
* @author alteredq / http://alteredqualia.com/
*/

function Points( geometry, material ) {
  Object3D.call( this );

  this.type = 'Points';

  this.geometry = geometry !== undefined ? geometry : new BufferGeometry();
  this.material = material;
}

Points.prototype = Object.assign( Object.create( Object3D.prototype ), {

  constructor: Points,

  isPoints: true,

} );

// kept for compatibility with THREE
function AmbientLight() {}


/**
* @author zz85 / http://www.lab4games.net/zz85/blog
* Extensible curve object
*
**/

/**************************************************************
* Abstract Curve base class
**************************************************************/

function Curve() {}

Curve.prototype = {

  constructor: Curve,

  // Get sequence of points using getPoint( t )

  getPoints: function ( divisions ) {
    let points = [];
    for ( let d = 0; d <= divisions; d ++ ) {
      points.push( this.getPoint( d / divisions ) );
    }
    return points;
  },

};

// TODO: Transformation for Curves?

/**************************************************************
* 3D Curves
**************************************************************/

// A Factory method for creating new curve subclasses

Curve.create = function ( constructor, getPointFunc ) {
  constructor.prototype = Object.create( Curve.prototype );
  constructor.prototype.constructor = constructor;
  constructor.prototype.getPoint = getPointFunc;

  return constructor;
};

/**
* @author zz85 https://github.com/zz85
*
* Centripetal CatmullRom Curve - which is useful for avoiding
* cusps and self-intersections in non-uniform catmull rom curves.
* http://www.cemyuksel.com/research/catmullrom_param/catmullrom.pdf
*
* curve.type accepts centripetal(default), chordal and catmullrom
* curve.tension is used for catmullrom which defaults to 0.5
*/

let CatmullRomCurve3 = ( function () {
let
  tmp = new Vector3(),
  px = new CubicPoly(),
  py = new CubicPoly(),
  pz = new CubicPoly();

/*
Based on an optimized c++ solution in
 - http://stackoverflow.com/questions/9489736/catmull-rom-curve-with-no-cusps-and-no-self-intersections/
 - http://ideone.com/NoEbVM

This CubicPoly class could be used for reusing some variables and calculations,
but for three.js curve use, it could be possible inlined and flatten into a single function call
which can be placed in CurveUtils.
*/

function CubicPoly() {}

/*
 * Compute coefficients for a cubic polynomial
 *   p(s) = c0 + c1*s + c2*s^2 + c3*s^3
 * such that
 *   p(0) = x0, p(1) = x1
 *  and
 *   p'(0) = t0, p'(1) = t1.
 */
CubicPoly.prototype.init = function ( x0, x1, t0, t1 ) {
  this.c0 = x0;
  this.c1 = t0;
  this.c2 = - 3 * x0 + 3 * x1 - 2 * t0 - t1;
  this.c3 = 2 * x0 - 2 * x1 + t0 + t1;
};

CubicPoly.prototype.initNonuniformCatmullRom = function ( x0, x1, x2, x3, dt0, dt1, dt2 ) {
  // compute tangents when parameterized in [t1,t2]
  let t1 = ( x1 - x0 ) / dt0 - ( x2 - x0 ) / ( dt0 + dt1 ) + ( x2 - x1 ) / dt1;
  let t2 = ( x2 - x1 ) / dt1 - ( x3 - x1 ) / ( dt1 + dt2 ) + ( x3 - x2 ) / dt2;

  // rescale tangents for parametrization in [0,1]
  t1 *= dt1;
  t2 *= dt1;

  // initCubicPoly
  this.init( x1, x2, t1, t2 );
};

CubicPoly.prototype.calc = function ( t ) {
  let t2 = t * t;
  let t3 = t2 * t;
  return this.c0 + this.c1 * t + this.c2 * t2 + this.c3 * t3;
};

// Subclass Three.js curve
return Curve.create(

  function ( p /* array of Vector3 */ ) {
    this.points = p || [];
    this.closed = false;
  },

  function ( t ) {
    let points = this.points,
      point, intPoint, weight, l;

    l = points.length;

    if ( l < 2 ) console.log( 'duh, you need at least 2 points' );

    point = ( l - ( this.closed ? 0 : 1 ) ) * t;
    intPoint = Math.floor( point );
    weight = point - intPoint;

    if ( this.closed ) {
      intPoint += intPoint > 0 ? 0 : ( Math.floor( Math.abs( intPoint ) / points.length ) + 1 ) * points.length;
    } else if ( weight === 0 && intPoint === l - 1 ) {
      intPoint = l - 2;
      weight = 1;
    }

    let p0, p1, p2, p3; // 4 points

    if ( this.closed || intPoint > 0 ) {
      p0 = points[( intPoint - 1 ) % l];
    } else {
      // extrapolate first point
      tmp.subVectors( points[0], points[1] ).add( points[0] );
      p0 = tmp;
    }

    p1 = points[intPoint % l];
    p2 = points[( intPoint + 1 ) % l];

    if ( this.closed || intPoint + 2 < l ) {
      p3 = points[( intPoint + 2 ) % l];
    } else {
      // extrapolate last point
      tmp.subVectors( points[l - 1], points[l - 2] ).add( points[l - 1] );
      p3 = tmp;
    }

    // init Centripetal Catmull-Rom
    let pow = 0.25;
    let dt0 = Math.pow( p0.distanceToSquared( p1 ), pow );
    let dt1 = Math.pow( p1.distanceToSquared( p2 ), pow );
    let dt2 = Math.pow( p2.distanceToSquared( p3 ), pow );

    // safety check for repeated points
    if ( dt1 < 1e-4 ) dt1 = 1.0;
    if ( dt0 < 1e-4 ) dt0 = dt1;
    if ( dt2 < 1e-4 ) dt2 = dt1;

    px.initNonuniformCatmullRom( p0.x, p1.x, p2.x, p3.x, dt0, dt1, dt2 );
    py.initNonuniformCatmullRom( p0.y, p1.y, p2.y, p3.y, dt0, dt1, dt2 );
    pz.initNonuniformCatmullRom( p0.z, p1.z, p2.z, p3.z, dt0, dt1, dt2 );

    let v = new Vector3(
      px.calc( weight ),
      py.calc( weight ),
      pz.calc( weight )
    );

    return v;
  }

);
} )();

export {
  WebGLRenderer,
  Fog,
  Scene,
  Mesh,
  LineSegments,
  Line,
  Points,
  ShaderMaterial,
  AmbientLight,
  OrthographicCamera,
  BufferGeometry,
  BufferAttribute,
  Object3D,
  Ray,
  Matrix4,
  Vector3,
  Quaternion,
  Color,
  CatmullRomCurve3,
  TriangleStripDrawMode,
  Texture,
};
