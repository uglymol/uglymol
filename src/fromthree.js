// a small subset of THREE.js v83 that is used by UglyMol
// (reformatted with eslint --fix)
// threejs.org/license

/* eslint-disable max-len, one-var, guard-for-in */
/* eslint-disable prefer-rest-params, no-invalid-this */
/* eslint-disable new-cap, no-unused-vars, no-undef */
/* eslint-disable no-extend-native, no-useless-escape, no-throw-literal */

// Polyfills

if ( Function.prototype.name === undefined ) {
  // Missing in IE9-11.
  // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Function/name
  Object.defineProperty( Function.prototype, 'name', {
    get: function () {
      return this.toString().match( /^\s*function\s*([^\(\s]*)/ )[1];
    },
  } );
}

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

let REVISION = '83mod';
let SmoothShading = 2;
let NoColors = 0;
let VertexColors = 2;
let NoBlending = 0;
let NormalBlending = 1;
let AdditiveBlending = 2;
let SubtractiveBlending = 3;
let MultiplyBlending = 4;
let CustomBlending = 5;
let AddEquation = 100;
let SubtractEquation = 101;
let ReverseSubtractEquation = 102;
let MinEquation = 103;
let MaxEquation = 104;
let ZeroFactor = 200;
let OneFactor = 201;
let SrcColorFactor = 202;
let OneMinusSrcColorFactor = 203;
let SrcAlphaFactor = 204;
let OneMinusSrcAlphaFactor = 205;
let DstAlphaFactor = 206;
let OneMinusDstAlphaFactor = 207;
let DstColorFactor = 208;
let OneMinusDstColorFactor = 209;
let SrcAlphaSaturateFactor = 210;
let NeverDepth = 0;
let AlwaysDepth = 1;
let LessDepth = 2;
let LessEqualDepth = 3;
let EqualDepth = 4;
let GreaterEqualDepth = 5;
let GreaterDepth = 6;
let NotEqualDepth = 7;
let UVMapping = 300;
let RepeatWrapping = 1000;
let ClampToEdgeWrapping = 1001;
let MirroredRepeatWrapping = 1002;
let NearestFilter = 1003;
let NearestMipMapNearestFilter = 1004;
let NearestMipMapLinearFilter = 1005;
let LinearFilter = 1006;
let LinearMipMapNearestFilter = 1007;
let LinearMipMapLinearFilter = 1008;
let UnsignedByteType = 1009;
let ByteType = 1010;
let ShortType = 1011;
let UnsignedShortType = 1012;
let IntType = 1013;
let UnsignedIntType = 1014;
let FloatType = 1015;
let HalfFloatType = 1016;
let UnsignedShort4444Type = 1017;
let UnsignedShort5551Type = 1018;
let UnsignedShort565Type = 1019;
let UnsignedInt248Type = 1020;
let AlphaFormat = 1021;
let RGBFormat = 1022;
let RGBAFormat = 1023;
let LuminanceFormat = 1024;
let LuminanceAlphaFormat = 1025;
let DepthFormat = 1026;
let DepthStencilFormat = 1027;
let RGB_S3TC_DXT1_Format = 2001;
let RGBA_S3TC_DXT1_Format = 2002;
let RGBA_S3TC_DXT3_Format = 2003;
let RGBA_S3TC_DXT5_Format = 2004;
let RGB_PVRTC_4BPPV1_Format = 2100;
let RGB_PVRTC_2BPPV1_Format = 2101;
let RGBA_PVRTC_4BPPV1_Format = 2102;
let RGBA_PVRTC_2BPPV1_Format = 2103;
let RGB_ETC1_Format = 2151;
let TrianglesDrawMode = 0;
let TriangleStripDrawMode = 1;
let TriangleFanDrawMode = 2;
let BasicDepthPacking = 3200;
let RGBADepthPacking = 3201;

/**
* @author alteredq / http://alteredqualia.com/
* @author mrdoob / http://mrdoob.com/
*/

let _Math = {

  DEG2RAD: Math.PI / 180,
  RAD2DEG: 180 / Math.PI,

  generateUUID: function () {
    // http://www.broofa.com/Tools/Math.uuid.htm

    let chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'.split( '' );
    let uuid = new Array( 36 );
    let rnd = 0, r;

    return function generateUUID() {
      for ( let i = 0; i < 36; i ++ ) {
        if ( i === 8 || i === 13 || i === 18 || i === 23 ) {
          uuid[i] = '-';
        } else if ( i === 14 ) {
          uuid[i] = '4';
        } else {
          if ( rnd <= 0x02 ) rnd = 0x2000000 + ( Math.random() * 0x1000000 ) | 0;
          r = rnd & 0xf;
          rnd = rnd >> 4;
          uuid[i] = chars[( i === 19 ) ? ( r & 0x3 ) | 0x8 : r];
        }
      }

      return uuid.join( '' );
    };
  }(),

  clamp: function ( value, min, max ) {
    return Math.max( min, Math.min( max, value ) );
  },

  // compute euclidian modulo of m % n
  // https://en.wikipedia.org/wiki/Modulo_operation

  euclideanModulo: function ( n, m ) {
    return ( ( n % m ) + m ) % m;
  },

  isPowerOfTwo: function ( value ) {
    return ( value & ( value - 1 ) ) === 0 && value !== 0;
  },

};

/**
* @author mikael emtinger / http://gomo.se/
* @author alteredq / http://alteredqualia.com/
* @author WestLangley / http://github.com/WestLangley
* @author bhouston / http://clara.io
*/

function Quaternion( x, y, z, w ) {
  this._x = x || 0;
  this._y = y || 0;
  this._z = z || 0;
  this._w = ( w !== undefined ) ? w : 1;
}

Quaternion.prototype = {

  constructor: Quaternion,

  get x() {
    return this._x;
  },

  get y() {
    return this._y;
  },

  get z() {
    return this._z;
  },

  get w() {
    return this._w;
  },

  setFromRotationMatrix: function ( m ) {
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm

    // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

    let te = m.elements,

      m11 = te[0], m12 = te[4], m13 = te[8],
      m21 = te[1], m22 = te[5], m23 = te[9],
      m31 = te[2], m32 = te[6], m33 = te[10],

      trace = m11 + m22 + m33,
      s;

    if ( trace > 0 ) {
      s = 0.5 / Math.sqrt( trace + 1.0 );

      this._w = 0.25 / s;
      this._x = ( m32 - m23 ) * s;
      this._y = ( m13 - m31 ) * s;
      this._z = ( m21 - m12 ) * s;
    } else if ( m11 > m22 && m11 > m33 ) {
      s = 2.0 * Math.sqrt( 1.0 + m11 - m22 - m33 );

      this._w = ( m32 - m23 ) / s;
      this._x = 0.25 * s;
      this._y = ( m12 + m21 ) / s;
      this._z = ( m13 + m31 ) / s;
    } else if ( m22 > m33 ) {
      s = 2.0 * Math.sqrt( 1.0 + m22 - m11 - m33 );

      this._w = ( m13 - m31 ) / s;
      this._x = ( m12 + m21 ) / s;
      this._y = 0.25 * s;
      this._z = ( m23 + m32 ) / s;
    } else {
      s = 2.0 * Math.sqrt( 1.0 + m33 - m11 - m22 );

      this._w = ( m21 - m12 ) / s;
      this._x = ( m13 + m31 ) / s;
      this._y = ( m23 + m32 ) / s;
      this._z = 0.25 * s;
    }

    return this;
  },

  setFromUnitVectors: function () {
    // http://lolengine.net/blog/2014/02/24/quaternion-from-two-vectors-final

    // assumes direction vectors vFrom and vTo are normalized

    let v1, r;

    let EPS = 0.000001;

    return function setFromUnitVectors( vFrom, vTo ) {
      if ( v1 === undefined ) v1 = new Vector3();

      r = vFrom.dot( vTo ) + 1;

      if ( r < EPS ) {
        r = 0;

        if ( Math.abs( vFrom.x ) > Math.abs( vFrom.z ) ) {
          v1.set( - vFrom.y, vFrom.x, 0 );
        } else {
          v1.set( 0, - vFrom.z, vFrom.y );
        }
      } else {
        v1.crossVectors( vFrom, vTo );
      }

      this._x = v1.x;
      this._y = v1.y;
      this._z = v1.z;
      this._w = r;

      return this.normalize();
    };
  }(),

  length: function () {
    return Math.sqrt( this._x * this._x + this._y * this._y + this._z * this._z + this._w * this._w );
  },

  normalize: function () {
    let l = this.length();

    if ( l === 0 ) {
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
  },

};

/**
* @author mrdoob / http://mrdoob.com/
* @author philogb / http://blog.thejit.org/
* @author egraether / http://egraether.com/
* @author zz85 / http://www.lab4games.net/zz85/blog
*/

function Vector2( x, y ) {
  this.x = x || 0;
  this.y = y || 0;
}

Vector2.prototype = {

  constructor: Vector2,

  isVector2: true,

  set: function ( x, y ) {
    this.x = x;
    this.y = y;

    return this;
  },

  clone: function () {
    return new this.constructor( this.x, this.y );
  },

  copy: function ( v ) {
    this.x = v.x;
    this.y = v.y;

    return this;
  },

  equals: function ( v ) {
    return ( ( v.x === this.x ) && ( v.y === this.y ) );
  },
};

/**
* @author mrdoob / http://mrdoob.com/
* @author *kile / http://kile.stravaganza.org/
* @author philogb / http://blog.thejit.org/
* @author mikael emtinger / http://gomo.se/
* @author egraether / http://egraether.com/
* @author WestLangley / http://github.com/WestLangley
*/

function Vector3( x, y, z ) {
  this.x = x || 0;
  this.y = y || 0;
  this.z = z || 0;
}

Vector3.prototype = {

  constructor: Vector3,

  isVector3: true,

  set: function ( x, y, z ) {
    this.x = x;
    this.y = y;
    this.z = z;

    return this;
  },

  clone: function () {
    return new this.constructor( this.x, this.y, this.z );
  },

  copy: function ( v ) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;

    return this;
  },

  add: function ( v ) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;

    return this;
  },

  addVectors: function ( a, b ) {
    this.x = a.x + b.x;
    this.y = a.y + b.y;
    this.z = a.z + b.z;

    return this;
  },

  addScaledVector: function ( v, s ) {
    this.x += v.x * s;
    this.y += v.y * s;
    this.z += v.z * s;

    return this;
  },

  sub: function ( v ) {
    this.x -= v.x;
    this.y -= v.y;
    this.z -= v.z;

    return this;
  },

  subVectors: function ( a, b ) {
    this.x = a.x - b.x;
    this.y = a.y - b.y;
    this.z = a.z - b.z;

    return this;
  },

  multiplyScalar: function ( scalar ) {
    if ( isFinite( scalar ) ) {
      this.x *= scalar;
      this.y *= scalar;
      this.z *= scalar;
    } else {
      this.x = 0;
      this.y = 0;
      this.z = 0;
    }

    return this;
  },

  applyMatrix4: function ( m ) {
    // input: THREE.Matrix4 affine matrix

    let x = this.x, y = this.y, z = this.z;
    let e = m.elements;

    this.x = e[0] * x + e[4] * y + e[8] * z + e[12];
    this.y = e[1] * x + e[5] * y + e[9] * z + e[13];
    this.z = e[2] * x + e[6] * y + e[10] * z + e[14];

    return this;
  },

  applyProjection: function ( m ) {
    // input: THREE.Matrix4 projection matrix

    let x = this.x, y = this.y, z = this.z;
    let e = m.elements;
    let d = 1 / ( e[3] * x + e[7] * y + e[11] * z + e[15] ); // perspective divide

    this.x = ( e[0] * x + e[4] * y + e[8] * z + e[12] ) * d;
    this.y = ( e[1] * x + e[5] * y + e[9] * z + e[13] ) * d;
    this.z = ( e[2] * x + e[6] * y + e[10] * z + e[14] ) * d;

    return this;
  },

  applyQuaternion: function ( q ) {
    let x = this.x, y = this.y, z = this.z;
    let qx = q.x, qy = q.y, qz = q.z, qw = q.w;

    // calculate quat * vector

    let ix = qw * x + qy * z - qz * y;
    let iy = qw * y + qz * x - qx * z;
    let iz = qw * z + qx * y - qy * x;
    let iw = - qx * x - qy * y - qz * z;

    // calculate result * inverse quat

    this.x = ix * qw + iw * - qx + iy * - qz - iz * - qy;
    this.y = iy * qw + iw * - qy + iz * - qx - ix * - qz;
    this.z = iz * qw + iw * - qz + ix * - qy - iy * - qx;

    return this;
  },

  unproject: function () {
    let matrix;

    return function unproject( camera ) {
      if ( matrix === undefined ) matrix = new Matrix4();

      matrix.multiplyMatrices( camera.matrixWorld, matrix.getInverse( camera.projectionMatrix ) );
      return this.applyProjection( matrix );
    };
  }(),

  transformDirection: function ( m ) {
    // input: THREE.Matrix4 affine matrix
    // vector interpreted as a direction

    let x = this.x, y = this.y, z = this.z;
    let e = m.elements;

    this.x = e[0] * x + e[4] * y + e[8] * z;
    this.y = e[1] * x + e[5] * y + e[9] * z;
    this.z = e[2] * x + e[6] * y + e[10] * z;

    return this.normalize();
  },

  divideScalar: function ( scalar ) {
    return this.multiplyScalar( 1 / scalar );
  },

  dot: function ( v ) {
    return this.x * v.x + this.y * v.y + this.z * v.z;
  },

  lengthSq: function () {
    return this.x * this.x + this.y * this.y + this.z * this.z;
  },

  length: function () {
    return Math.sqrt( this.x * this.x + this.y * this.y + this.z * this.z );
  },

  normalize: function () {
    return this.divideScalar( this.length() );
  },

  setLength: function ( length ) {
    return this.multiplyScalar( length / this.length() );
  },

  lerp: function ( v, alpha ) {
    this.x += ( v.x - this.x ) * alpha;
    this.y += ( v.y - this.y ) * alpha;
    this.z += ( v.z - this.z ) * alpha;

    return this;
  },

  cross: function ( v ) {
    let x = this.x, y = this.y, z = this.z;

    this.x = y * v.z - z * v.y;
    this.y = z * v.x - x * v.z;
    this.z = x * v.y - y * v.x;

    return this;
  },

  crossVectors: function ( a, b ) {
    let ax = a.x, ay = a.y, az = a.z;
    let bx = b.x, by = b.y, bz = b.z;

    this.x = ay * bz - az * by;
    this.y = az * bx - ax * bz;
    this.z = ax * by - ay * bx;

    return this;
  },

  distanceTo: function ( v ) {
    return Math.sqrt( this.distanceToSquared( v ) );
  },

  distanceToSquared: function ( v ) {
    let dx = this.x - v.x, dy = this.y - v.y, dz = this.z - v.z;
    return dx * dx + dy * dy + dz * dz;
  },

  setFromMatrixPosition: function ( m ) {
    return this.setFromMatrixColumn( m, 3 );
  },

  setFromMatrixColumn: function ( m, index ) {
    return this.fromArray( m.elements, index * 4 );
  },

  equals: function ( v ) {
    return ( ( v.x === this.x ) && ( v.y === this.y ) && ( v.z === this.z ) );
  },

  fromArray: function ( array, offset ) {
    if ( offset === undefined ) offset = 0;

    this.x = array[offset];
    this.y = array[offset + 1];
    this.z = array[offset + 2];

    return this;
  },
};

/**
* @author alteredq / http://alteredqualia.com/
* @author WestLangley / http://github.com/WestLangley
* @author bhouston / http://clara.io
* @author tschw
*/

function Matrix3() {
  this.elements = new Float32Array( [

    1, 0, 0,
    0, 1, 0,
    0, 0, 1,

  ] );

  if ( arguments.length > 0 ) {
    console.error( 'THREE.Matrix3: the constructor no longer reads arguments. use .set() instead.' );
  }
}

Matrix3.prototype = {

  constructor: Matrix3,

  isMatrix3: true,

  set: function ( n11, n12, n13, n21, n22, n23, n31, n32, n33 ) {
    let te = this.elements;

    te[0] = n11; te[1] = n21; te[2] = n31;
    te[3] = n12; te[4] = n22; te[5] = n32;
    te[6] = n13; te[7] = n23; te[8] = n33;

    return this;
  },

  setFromMatrix4: function ( m ) {
    let me = m.elements;

    this.set(

      me[0], me[4], me[8],
      me[1], me[5], me[9],
      me[2], me[6], me[10]

    );

    return this;
  },

  getInverse: function ( matrix, throwOnDegenerate ) {
    if ( matrix && matrix.isMatrix4 ) {
      console.error( 'THREE.Matrix3.getInverse no longer takes a Matrix4 argument.' );
    }

    let me = matrix.elements,
      te = this.elements,

      n11 = me[0], n21 = me[1], n31 = me[2],
      n12 = me[3], n22 = me[4], n32 = me[5],
      n13 = me[6], n23 = me[7], n33 = me[8],

      t11 = n33 * n22 - n32 * n23,
      t12 = n32 * n13 - n33 * n12,
      t13 = n23 * n12 - n22 * n13,

      det = n11 * t11 + n21 * t12 + n31 * t13;

    if ( det === 0 ) {
      let msg = 'THREE.Matrix3.getInverse(): can\'t invert matrix, determinant is 0';

      if ( throwOnDegenerate === true ) {
        throw new Error( msg );
      } else {
        console.warn( msg );
      }

      return this.identity();
    }

    let detInv = 1 / det;

    te[0] = t11 * detInv;
    te[1] = ( n31 * n23 - n33 * n21 ) * detInv;
    te[2] = ( n32 * n21 - n31 * n22 ) * detInv;

    te[3] = t12 * detInv;
    te[4] = ( n33 * n11 - n31 * n13 ) * detInv;
    te[5] = ( n31 * n12 - n32 * n11 ) * detInv;

    te[6] = t13 * detInv;
    te[7] = ( n21 * n13 - n23 * n11 ) * detInv;
    te[8] = ( n22 * n11 - n21 * n12 ) * detInv;

    return this;
  },

  transpose: function () {
    let tmp, m = this.elements;

    tmp = m[1]; m[1] = m[3]; m[3] = tmp;
    tmp = m[2]; m[2] = m[6]; m[6] = tmp;
    tmp = m[5]; m[5] = m[7]; m[7] = tmp;

    return this;
  },

  getNormalMatrix: function ( matrix4 ) {
    return this.setFromMatrix4( matrix4 ).getInverse( this ).transpose();
  },
};

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

  if ( arguments.length > 0 ) {
    console.error( 'THREE.Matrix4: the constructor no longer reads arguments. use .set() instead.' );
  }
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

  getInverse: function ( m, throwOnDegenerate ) {
    // based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
    let te = this.elements,
      me = m.elements,

      n11 = me[0], n21 = me[1], n31 = me[2], n41 = me[3],
      n12 = me[4], n22 = me[5], n32 = me[6], n42 = me[7],
      n13 = me[8], n23 = me[9], n33 = me[10], n43 = me[11],
      n14 = me[12], n24 = me[13], n34 = me[14], n44 = me[15],

      t11 = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44,
      t12 = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44,
      t13 = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44,
      t14 = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;

    let det = n11 * t11 + n21 * t12 + n31 * t13 + n41 * t14;

    if ( det === 0 ) {
      let msg = 'THREE.Matrix4.getInverse(): can\'t invert matrix, determinant is 0';

      if ( throwOnDegenerate === true ) {
        throw new Error( msg );
      } else {
        console.warn( msg );
      }

      return this.identity();
    }

    let detInv = 1 / det;

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

function EventDispatcher() {}

Object.assign( EventDispatcher.prototype, {

  addEventListener: function ( type, listener ) {
    if ( this._listeners === undefined ) this._listeners = {};

    let listeners = this._listeners;

    if ( listeners[type] === undefined ) {
      listeners[type] = [];
    }

    if ( listeners[type].indexOf( listener ) === - 1 ) {
      listeners[type].push( listener );
    }
  },

  hasEventListener: function ( type, listener ) {
    if ( this._listeners === undefined ) return false;

    let listeners = this._listeners;

    return listeners[type] !== undefined && listeners[type].indexOf( listener ) !== - 1;
  },

  removeEventListener: function ( type, listener ) {
    if ( this._listeners === undefined ) return;

    let listeners = this._listeners;
    let listenerArray = listeners[type];

    if ( listenerArray !== undefined ) {
      let index = listenerArray.indexOf( listener );

      if ( index !== - 1 ) {
        listenerArray.splice( index, 1 );
      }
    }
  },

  dispatchEvent: function ( event ) {
    if ( this._listeners === undefined ) return;

    let listeners = this._listeners;
    let listenerArray = listeners[event.type];

    if ( listenerArray !== undefined ) {
      event.target = this;

      let array = [], i = 0;
      let length = listenerArray.length;

      for ( i = 0; i < length; i ++ ) {
        array[i] = listenerArray[i];
      }

      for ( i = 0; i < length; i ++ ) {
        array[i].call( this, event );
      }
    }
  },

} );

/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
* @author szimek / https://github.com/szimek/
*/

let textureId = 0;

function Texture( image, mapping, wrapS, wrapT, magFilter, minFilter, format, type, anisotropy, encoding ) {
  Object.defineProperty( this, 'id', { value: textureId ++ } );

  this.uuid = _Math.generateUUID();

  this.name = '';

  this.image = image !== undefined ? image : Texture.DEFAULT_IMAGE;
  this.mipmaps = [];

  this.mapping = mapping !== undefined ? mapping : Texture.DEFAULT_MAPPING;

  this.wrapS = wrapS !== undefined ? wrapS : ClampToEdgeWrapping;
  this.wrapT = wrapT !== undefined ? wrapT : ClampToEdgeWrapping;

  this.magFilter = magFilter !== undefined ? magFilter : LinearFilter;
  this.minFilter = minFilter !== undefined ? minFilter : LinearMipMapLinearFilter;

  this.anisotropy = anisotropy !== undefined ? anisotropy : 1;

  this.format = format !== undefined ? format : RGBAFormat;
  this.type = type !== undefined ? type : UnsignedByteType;

  this.offset = new Vector2( 0, 0 );
  this.repeat = new Vector2( 1, 1 );

  this.generateMipmaps = true;
  this.premultiplyAlpha = false;
  this.flipY = true;
  this.unpackAlignment = 4; // valid values: 1, 2, 4, 8 (see http://www.khronos.org/opengles/sdk/docs/man/xhtml/glPixelStorei.xml)

  this.version = 0;
  this.onUpdate = null;
}

Texture.DEFAULT_IMAGE = undefined;
Texture.DEFAULT_MAPPING = UVMapping;

Texture.prototype = {

  constructor: Texture,

  isTexture: true,

  set needsUpdate( value ) {
    if ( value === true ) this.version ++;
  },

  dispose: function () {
    this.dispatchEvent( { type: 'dispose' } );
  },
};

Object.assign( Texture.prototype, EventDispatcher.prototype );


/**
* @author tschw
*
* Uniforms of a program.
* Those form a tree structure with a special top-level container for the root,
* which you get by calling 'new WebGLUniforms( gl, program, renderer )'.
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
* .setValue( gl, value, [renderer] )
*
*     uploads a uniform value(s)
*   the 'renderer' parameter is needed for sampler uniforms
*
*
* Static methods of the top-level container (renderer factorizations):
*
* .upload( gl, seq, values, renderer )
*
*     sets uniforms in 'seq' to 'values[id].value'
*
* .seqWithValue( seq, values ) : filteredSeq
*
*     filters 'seq' entries with corresponding entry in values
*
*
* Methods of the top-level container (renderer factorizations):
*
* .setValue( gl, name, value )
*
*     sets uniform with  name 'name' to 'value'
*
* .set( gl, obj, prop )
*
*     sets uniform from object and property with same name than uniform
*
* .setOptional( gl, obj, prop )
*
*     like .set for an optional property of the object
*
*/

let emptyTexture = new Texture();

// --- Base for inner nodes (including the root) ---

function UniformContainer() {
  this.seq = [];
  this.map = {};
}

// --- Utilities ---

// Array Caches (provide typed arrays for temporary by size)

let arrayCacheF32 = [];
let arrayCacheI32 = [];

// --- Setters ---

// Note: Defining these methods externally, because they come in a bunch
// and this way their names minify.

// Single scalar

function setValue1f( gl, v ) { gl.uniform1f( this.addr, v ); }
function setValue1i( gl, v ) { gl.uniform1i( this.addr, v ); }

// Single float vector (from flat array or THREE.VectorN)

function setValue2fv( gl, v ) {
  if ( v.x === undefined ) gl.uniform2fv( this.addr, v );
  else gl.uniform2f( this.addr, v.x, v.y );
}

function setValue3fv( gl, v ) {
  if ( v.x !== undefined ) { gl.uniform3f( this.addr, v.x, v.y, v.z ); } else if ( v.r !== undefined ) { gl.uniform3f( this.addr, v.r, v.g, v.b ); } else { gl.uniform3fv( this.addr, v ); }
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

function SingleUniform( id, activeInfo, addr ) {
  this.id = id;
  this.addr = addr;
  this.setValue = getSingularSetter( activeInfo.type );

  // this.path = activeInfo.name; // DEBUG
}

// --- Top-level ---

// Parser - builds up the property tree from the path strings

let RePathPart = /([\w\d_]+)(\])?(\[|\.)?/g;

// extracts
//  - the identifier (member name or array index)
//  - followed by an optional right bracket (found when array index)
//  - followed by an optional left bracket or dot (type of subscript)
//
// Note: These portions can be read in a non-overlapping fashion and
// allow straightforward parsing of the hierarchy that WebGL encodes
// in the uniform names.

function addUniform( container, uniformObject ) {
  container.seq.push( uniformObject );
  container.map[uniformObject.id] = uniformObject;
}

function parseUniform( activeInfo, addr, container ) {
  let path = activeInfo.name,
    pathLength = path.length;

  // reset RegExp object, because of the early exit of a previous run
  RePathPart.lastIndex = 0;

  for (; ;) {
    let match = RePathPart.exec( path ),
      matchEnd = RePathPart.lastIndex,

      id = match[1],
      idIsIndex = match[2] === ']',
      subscript = match[3];

    if ( idIsIndex ) id = id | 0; // convert to integer

    if ( subscript === undefined ||
      subscript === '[' && matchEnd + 2 === pathLength ) {
      // bare name or "pure" bottom-level array "[0]" suffix

      addUniform( container, subscript === undefined ?
        new SingleUniform( id, activeInfo, addr ) :
        new PureArrayUniform( id, activeInfo, addr ) );

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

/**
* Uniform Utilities
*/

let UniformsUtils = {

  merge: function ( uniforms ) {
    let merged = {};

    for ( let u = 0; u < uniforms.length; u ++ ) {
      let tmp = this.clone( uniforms[u] );

      for ( let p in tmp ) {
        merged[p] = tmp[p];
      }
    }

    return merged;
  },

  clone: function ( uniforms_src ) {
    let uniforms_dst = {};

    for ( let u in uniforms_src ) {
      uniforms_dst[u] = {};

      for ( let p in uniforms_src[u] ) {
        let parameter_src = uniforms_src[u][p];

        if ( parameter_src && ( parameter_src.isColor ||
        parameter_src.isMatrix3 || parameter_src.isMatrix4 ||
        parameter_src.isVector2 || parameter_src.isVector3 || parameter_src.isVector4 ||
        parameter_src.isTexture ) ) {
          uniforms_dst[u][p] = parameter_src.clone();
        } else if ( Array.isArray( parameter_src ) ) {
          uniforms_dst[u][p] = parameter_src.slice();
        } else {
          uniforms_dst[u][p] = parameter_src;
        }
      }
    }

    return uniforms_dst;
  },

};


let begin_vertex = `
vec3 transformed = vec3( position );
`;

let color_fragment = `
#ifdef USE_COLOR
diffuseColor.rgb *= vColor;
#endif`;

let color_pars_fragment = `
#ifdef USE_COLOR
varying vec3 vColor;
#endif
`;

let color_pars_vertex = `
#ifdef USE_COLOR
varying vec3 vColor;
#endif`;

let color_vertex = `
#ifdef USE_COLOR
vColor.xyz = color.xyz;
#endif`;

let fog_fragment = `
#ifdef USE_FOG
float depth = gl_FragCoord.z / gl_FragCoord.w;
float fogFactor = smoothstep( fogNear, fogFar, depth );
gl_FragColor.rgb = mix( gl_FragColor.rgb, fogColor, fogFactor );
#endif
`;

let fog_pars_fragment = `
#ifdef USE_FOG
uniform vec3 fogColor;
uniform float fogNear;
uniform float fogFar;
#endif`;


let packing = `
vec3 packNormalToRGB( const in vec3 normal ) {
return normalize( normal ) * 0.5 + 0.5;
}

vec3 unpackRGBToNormal( const in vec3 rgb ) {
return 1.0 - 2.0 * rgb.xyz;
}

const float PackUpscale = 256. / 255.; // fraction -> 0..1 (including 1)
const float UnpackDownscale = 255. / 256.; // 0..1 -> fraction (excluding 1)

const vec3 PackFactors = vec3( 256. * 256. * 256., 256. * 256.,  256. );
const vec4 UnpackFactors = UnpackDownscale / vec4( PackFactors, 1. );

const float ShiftRight8 = 1. / 256.;

vec4 packDepthToRGBA( const in float v ) {

vec4 r = vec4( fract( v * PackFactors ), v );
r.yzw -= r.xyz * ShiftRight8; // tidy overflow
return r * PackUpscale;

}

float unpackRGBAToDepth( const in vec4 v ) {

return dot( v, UnpackFactors );

}

// NOTE: viewZ/eyeZ is < 0 when in front of the camera per OpenGL conventions

float viewZToOrthographicDepth( const in float viewZ, const in float near, const in float far ) {
return ( viewZ + near ) / ( near - far );
}
float orthographicDepthToViewZ( const in float linearClipZ, const in float near, const in float far ) {
return linearClipZ * ( near - far ) - near;
}

float viewZToPerspectiveDepth( const in float viewZ, const in float near, const in float far ) {
return (( near + viewZ ) * far ) / (( far - near ) * viewZ );
}
float perspectiveDepthToViewZ( const in float invClipZ, const in float near, const in float far ) {
return ( near * far ) / ( ( far - near ) * invClipZ - far );
}
`;

let premultiplied_alpha_fragment = `
#ifdef PREMULTIPLIED_ALPHA

// Get get normal blending with premultipled, use with CustomBlending, OneFactor, OneMinusSrcAlphaFactor, AddEquation.
gl_FragColor.rgb *= gl_FragColor.a;

#endif
`;

let project_vertex = `
vec4 mvPosition = modelViewMatrix * vec4( transformed, 1.0 );

gl_Position = projectionMatrix * mvPosition;
`;


let depth_frag = `
#if DEPTH_PACKING == 3200

uniform float opacity;

#endif

#include <packing>

void main() {

vec4 diffuseColor = vec4( 1.0 );

#if DEPTH_PACKING == 3200

  diffuseColor.a = opacity;

#endif

#if DEPTH_PACKING == 3200

  gl_FragColor = vec4( vec3( gl_FragCoord.z ), opacity );

#elif DEPTH_PACKING == 3201

  gl_FragColor = packDepthToRGBA( gl_FragCoord.z );

#endif

}
`;

let depth_vert = `

void main() {

#include <begin_vertex>
#include <project_vertex>

}
`;


let meshbasic_frag = `
uniform vec3 diffuse;
uniform float opacity;

#include <fog_pars_fragment>

void main() {

gl_FragColor = vec4( diffuse, opacity );

#include <fog_fragment>

}
`;

let meshbasic_vert = `
#include <color_pars_vertex>

void main() {

#include <color_vertex>

#include <begin_vertex>
#include <project_vertex>

}
`;

let ShaderChunk = {
  begin_vertex: begin_vertex,
  color_fragment: color_fragment,
  color_pars_fragment: color_pars_fragment,
  color_pars_vertex: color_pars_vertex,
  color_vertex: color_vertex,
  fog_fragment: fog_fragment,
  fog_pars_fragment: fog_pars_fragment,
  packing: packing,
  premultiplied_alpha_fragment: premultiplied_alpha_fragment,
  project_vertex: project_vertex,

  depth_frag: depth_frag,
  depth_vert: depth_vert,
  meshbasic_frag: meshbasic_frag,
  meshbasic_vert: meshbasic_vert,
};

/**
* @author supereggbert / http://www.paulbrunt.co.uk/
* @author philogb / http://blog.thejit.org/
* @author mikael emtinger / http://gomo.se/
* @author egraether / http://egraether.com/
* @author WestLangley / http://github.com/WestLangley
*/

function Vector4( x, y, z, w ) {
  this.x = x || 0;
  this.y = y || 0;
  this.z = z || 0;
  this.w = ( w !== undefined ) ? w : 1;
}

Vector4.prototype = {

  constructor: Vector4,

  isVector4: true,

  set: function ( x, y, z, w ) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.w = w;

    return this;
  },

  clone: function () {
    return new this.constructor( this.x, this.y, this.z, this.w );
  },

  copy: function ( v ) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;
    this.w = ( v.w !== undefined ) ? v.w : 1;

    return this;
  },

  multiplyScalar: function ( scalar ) {
    if ( isFinite( scalar ) ) {
      this.x *= scalar;
      this.y *= scalar;
      this.z *= scalar;
      this.w *= scalar;
    } else {
      this.x = 0;
      this.y = 0;
      this.z = 0;
      this.w = 0;
    }

    return this;
  },

  equals: function ( v ) {
    return ( ( v.x === this.x ) && ( v.y === this.y ) && ( v.z === this.z ) && ( v.w === this.w ) );
  },
};

/**
* @author mrdoob / http://mrdoob.com/
*/

function Color( r, g, b ) {
  if ( g === undefined && b === undefined ) {
    // r is THREE.Color, hex or string
    return this.set( r );
  }

  return this.setRGB( r, g, b );
}

Color.prototype = {

  constructor: Color,

  isColor: true,

  r: 1, g: 1, b: 1,

  set: function ( value ) {
    if ( value && value.isColor ) {
      this.copy( value );
    } else if ( typeof value === 'number' ) {
      this.setHex( value );
    }

    return this;
  },

  setHex: function ( hex ) {
    hex = Math.floor( hex );

    this.r = ( hex >> 16 & 255 ) / 255;
    this.g = ( hex >> 8 & 255 ) / 255;
    this.b = ( hex & 255 ) / 255;

    return this;
  },

  setRGB: function ( r, g, b ) {
    this.r = r;
    this.g = g;
    this.b = b;

    return this;
  },

  setHSL: function () {
    function hue2rgb( p, q, t ) {
      if ( t < 0 ) t += 1;
      if ( t > 1 ) t -= 1;
      if ( t < 1 / 6 ) return p + ( q - p ) * 6 * t;
      if ( t < 1 / 2 ) return q;
      if ( t < 2 / 3 ) return p + ( q - p ) * 6 * ( 2 / 3 - t );
      return p;
    }

    return function setHSL( h, s, l ) {
      // h,s,l ranges are in 0.0 - 1.0
      h = _Math.euclideanModulo( h, 1 );
      s = _Math.clamp( s, 0, 1 );
      l = _Math.clamp( l, 0, 1 );

      if ( s === 0 ) {
        this.r = this.g = this.b = l;
      } else {
        let p = l <= 0.5 ? l * ( 1 + s ) : l + s - ( l * s );
        let q = ( 2 * l ) - p;

        this.r = hue2rgb( q, p, h + 1 / 3 );
        this.g = hue2rgb( q, p, h );
        this.b = hue2rgb( q, p, h - 1 / 3 );
      }

      return this;
    };
  }(),

  clone: function () {
    return new this.constructor( this.r, this.g, this.b );
  },

  copy: function ( color ) {
    this.r = color.r;
    this.g = color.g;
    this.b = color.b;

    return this;
  },

  getHex: function () {
    return ( this.r * 255 ) << 16 ^ ( this.g * 255 ) << 8 ^ ( this.b * 255 ) << 0;
  },

  getHexString: function () {
    return ( '000000' + this.getHex().toString( 16 ) ).slice( - 6 );
  },
};

/**
* Uniforms library for shared webgl shaders
*/

let UniformsLib = {

  common: {

    diffuse: { value: new Color( 0xeeeeee ) },
    opacity: { value: 1.0 },

    map: { value: null },
    offsetRepeat: { value: new Vector4( 0, 0, 1, 1 ) },
  },

  fog: {

    fogNear: { value: 1 },
    fogFar: { value: 2000 },
    fogColor: { value: new Color( 0xffffff ) },

  },

  lights: {
    ambientLightColor: { value: [] },
  },

  points: {

    diffuse: { value: new Color( 0xeeeeee ) },
    opacity: { value: 1.0 },
    size: { value: 1.0 },
    scale: { value: 1.0 },
    map: { value: null },
    offsetRepeat: { value: new Vector4( 0, 0, 1, 1 ) },

  },

};

/**
* @author alteredq / http://alteredqualia.com/
* @author mrdoob / http://mrdoob.com/
* @author mikael emtinger / http://gomo.se/
*/

let ShaderLib = {

  basic: {

    uniforms: UniformsUtils.merge( [
      UniformsLib.common,
      UniformsLib.fog,
    ] ),

    vertexShader: ShaderChunk.meshbasic_vert,
    fragmentShader: ShaderChunk.meshbasic_frag,

  },

  depth: {

    uniforms: UniformsUtils.merge( [
      UniformsLib.common,
    ] ),

    vertexShader: ShaderChunk.depth_vert,
    fragmentShader: ShaderChunk.depth_frag,

  },

};


/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
*/

let materialId = 0;

function Material() {
  Object.defineProperty( this, 'id', { value: materialId ++ } );

  this.uuid = _Math.generateUUID();

  this.name = '';
  this.type = 'Material';

  this.fog = true;
  this.lights = true;

  this.blending = NormalBlending;
  this.shading = SmoothShading; // THREE.FlatShading, THREE.SmoothShading
  this.vertexColors = NoColors; // THREE.NoColors, THREE.VertexColors, THREE.FaceColors

  this.opacity = 1;
  this.transparent = false;

  this.blendSrc = SrcAlphaFactor;
  this.blendDst = OneMinusSrcAlphaFactor;
  this.blendEquation = AddEquation;
  this.blendSrcAlpha = null;
  this.blendDstAlpha = null;
  this.blendEquationAlpha = null;

  this.depthFunc = LessEqualDepth;
  this.depthTest = true;
  this.depthWrite = true;

  this.colorWrite = true;

  this.precision = null; // override the renderer's default precision for this material

  this.polygonOffset = false;
  this.polygonOffsetFactor = 0;
  this.polygonOffsetUnits = 0;

  this.premultipliedAlpha = false;

  this.overdraw = 0; // Overdrawn pixels (typically between 0 and 1) for fixing antialiasing gaps in CanvasRenderer

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
    if ( values === undefined ) return;

    for ( let key in values ) {
      let newValue = values[key];

      if ( newValue === undefined ) {
        console.warn( 'THREE.Material: \'' + key + '\' parameter is undefined.' );
        continue;
      }

      let currentValue = this[key];

      if ( currentValue === undefined ) {
        console.warn( 'THREE.' + this.type + ': \'' + key + '\' is not a property of this material.' );
        continue;
      }

      if ( currentValue && currentValue.isColor ) {
        currentValue.set( newValue );
      } else if ( (currentValue && currentValue.isVector3) && (newValue && newValue.isVector3) ) {
        currentValue.copy( newValue );
      } else if ( key === 'overdraw' ) {
        // ensure overdraw is backwards-compatible with legacy boolean type
        this[key] = Number( newValue );
      } else {
        this[key] = newValue;
      }
    }
  },

  clone: function () {
    return new this.constructor().copy( this );
  },

  copy: function ( source ) {
    this.name = source.name;

    this.fog = source.fog;
    this.lights = source.lights;

    this.blending = source.blending;
    this.shading = source.shading;
    this.vertexColors = source.vertexColors;

    this.opacity = source.opacity;
    this.transparent = source.transparent;

    this.blendSrc = source.blendSrc;
    this.blendDst = source.blendDst;
    this.blendEquation = source.blendEquation;
    this.blendSrcAlpha = source.blendSrcAlpha;
    this.blendDstAlpha = source.blendDstAlpha;
    this.blendEquationAlpha = source.blendEquationAlpha;

    this.depthFunc = source.depthFunc;
    this.depthTest = source.depthTest;
    this.depthWrite = source.depthWrite;

    this.colorWrite = source.colorWrite;

    this.precision = source.precision;

    this.polygonOffset = source.polygonOffset;
    this.polygonOffsetFactor = source.polygonOffsetFactor;
    this.polygonOffsetUnits = source.polygonOffsetUnits;

    this.premultipliedAlpha = source.premultipliedAlpha;

    this.overdraw = source.overdraw;

    this.visible = source.visible;

    return this;
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
*  defines: { "label" : "value" },
*  uniforms: { "parameter1": { value: 1.0 }, "parameter2": { value2: 2 } },
*
*  fragmentShader: <string>,
*  vertexShader: <string>,
*
*  wireframe: <boolean>,
*  wireframeLinewidth: <float>,
*
*  lights: <bool>,
* }
*/

function ShaderMaterial( parameters ) {
  Material.call( this );

  this.type = 'ShaderMaterial';

  this.defines = {};
  this.uniforms = {};

  this.vertexShader = 'void main() {\n\tgl_Position = projectionMatrix * modelViewMatrix * vec4( position, 1.0 );\n}';
  this.fragmentShader = 'void main() {\n\tgl_FragColor = vec4( 1.0, 0.0, 0.0, 1.0 );\n}';

  this.linewidth = 1;

  this.wireframe = false;
  this.wireframeLinewidth = 1;

  this.fog = false; // set to use scene fog
  this.lights = false; // set to use scene lights

  this.extensions = {
    fragDepth: false, // set to use fragment depth values
    drawBuffers: false, // set to use draw buffers
  };

  // When rendered geometry doesn't include these attributes but the material does,
  // use these default values in WebGL. This avoids errors when buffer data is missing.
  this.defaultAttributeValues = {
    'color': [1, 1, 1],
    'uv': [0, 0],
    'uv2': [0, 0],
  };

  this.index0AttributeName = undefined;

  if ( parameters !== undefined ) {
    if ( parameters.attributes !== undefined ) {
      console.error( 'THREE.ShaderMaterial: attributes should now be defined in THREE.BufferGeometry instead.' );
    }

    this.setValues( parameters );
  }
}

ShaderMaterial.prototype = Object.create( Material.prototype );
ShaderMaterial.prototype.constructor = ShaderMaterial;

ShaderMaterial.prototype.isShaderMaterial = true;

ShaderMaterial.prototype.copy = function ( source ) {
  Material.prototype.copy.call( this, source );

  this.fragmentShader = source.fragmentShader;
  this.vertexShader = source.vertexShader;

  this.uniforms = UniformsUtils.clone( source.uniforms );

  this.defines = source.defines;

  this.wireframe = source.wireframe;
  this.wireframeLinewidth = source.wireframeLinewidth;

  this.lights = source.lights;

  this.extensions = source.extensions;

  return this;
};


/**
* @author bhouston / http://clara.io
*/

function Plane( normal, constant ) {
  this.normal = ( normal !== undefined ) ? normal : new Vector3( 1, 0, 0 );
  this.constant = ( constant !== undefined ) ? constant : 0;
}

Plane.prototype = {

  constructor: Plane,

  setComponents: function ( x, y, z, w ) {
    this.normal.set( x, y, z );
    this.constant = w;

    return this;
  },

  normalize: function () {
    // Note: will lead to a divide by zero if the plane is invalid.

    let inverseNormalLength = 1.0 / this.normal.length();
    this.normal.multiplyScalar( inverseNormalLength );
    this.constant *= inverseNormalLength;

    return this;
  },

  distanceToPoint: function ( point ) {
    return this.normal.dot( point ) + this.constant;
  },
};

/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
* @author bhouston / http://clara.io
*/

function Frustum( p0, p1, p2, p3, p4, p5 ) {
  this.planes = [

  ( p0 !== undefined ) ? p0 : new Plane(),
  ( p1 !== undefined ) ? p1 : new Plane(),
  ( p2 !== undefined ) ? p2 : new Plane(),
  ( p3 !== undefined ) ? p3 : new Plane(),
  ( p4 !== undefined ) ? p4 : new Plane(),
  ( p5 !== undefined ) ? p5 : new Plane(),

  ];
}

Frustum.prototype = {

  constructor: Frustum,

  setFromMatrix: function ( m ) {
    let planes = this.planes;
    let me = m.elements;
    let me0 = me[0], me1 = me[1], me2 = me[2], me3 = me[3];
    let me4 = me[4], me5 = me[5], me6 = me[6], me7 = me[7];
    let me8 = me[8], me9 = me[9], me10 = me[10], me11 = me[11];
    let me12 = me[12], me13 = me[13], me14 = me[14], me15 = me[15];

    planes[0].setComponents( me3 - me0, me7 - me4, me11 - me8, me15 - me12 ).normalize();
    planes[1].setComponents( me3 + me0, me7 + me4, me11 + me8, me15 + me12 ).normalize();
    planes[2].setComponents( me3 + me1, me7 + me5, me11 + me9, me15 + me13 ).normalize();
    planes[3].setComponents( me3 - me1, me7 - me5, me11 - me9, me15 - me13 ).normalize();
    planes[4].setComponents( me3 - me2, me7 - me6, me11 - me10, me15 - me14 ).normalize();
    planes[5].setComponents( me3 + me2, me7 + me6, me11 + me10, me15 + me14 ).normalize();

    return this;
  },

};


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

  this.uuid = _Math.generateUUID();

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
    normalMatrix: {
      value: new Matrix3(),
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

  this.onBeforeRender = function () {};
  this.onAfterRender = function () {};
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

    if ( object === this ) {
      console.error( 'THREE.Object3D.add: object can\'t be added as a child of itself.', object );
      return this;
    }

    if ( ( object && object.isObject3D ) ) {
      if ( object.parent !== null ) {
        object.parent.remove( object );
      }

      object.parent = this;
      object.dispatchEvent( { type: 'added' } );

      this.children.push( object );
    } else {
      console.error( 'THREE.Object3D.add: object not an instance of THREE.Object3D.', object );
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

  raycast: function () {},

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
    throw new TypeError( 'THREE.BufferAttribute: array should be a Typed Array.' );
  }

  this.uuid = _Math.generateUUID();

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

  getX: function ( index ) {
    return this.array[index * this.itemSize];
  },

  getY: function ( index ) {
    return this.array[index * this.itemSize + 1];
  },

  getZ: function ( index ) {
    return this.array[index * this.itemSize + 2];
  },
};


let count = 0;
function GeometryIdCount() { return count++; }

/**
* @author alteredq / http://alteredqualia.com/
* @author mrdoob / http://mrdoob.com/
*/

function BufferGeometry() {
  Object.defineProperty( this, 'id', { value: GeometryIdCount() } );

  this.uuid = _Math.generateUUID();

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

  addAttribute: function ( name, attribute ) {
    if ( ( attribute && attribute.isBufferAttribute ) === false && ( attribute && attribute.isInterleavedBufferAttribute ) === false ) {
      console.warn( 'THREE.BufferGeometry: .addAttribute() now expects ( name, attribute ).' );

      this.addAttribute( name, new BufferAttribute( arguments[1], arguments[2] ) );

      return;
    }

    if ( name === 'index' ) {
      console.warn( 'THREE.BufferGeometry.addAttribute: Use .setIndex() for index attribute.' );
      this.setIndex( attribute );

      return;
    }

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

  isOrthographicCamera: true,

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
    console.error( 'THREE.WebGLShader: Shader couldn\'t compile.' );
  }

  if ( gl.getShaderInfoLog( shader ) !== '' ) {
    let info = gl.getShaderInfoLog( shader );
    // workaround for https://github.com/mrdoob/three.js/issues/9716
    if (info.indexOf('GL_ARB_gpu_shader5') === -1) {
      console.warn( 'THREE.WebGLShader: gl.getShaderInfoLog()', type === gl.VERTEX_SHADER ? 'vertex' : 'fragment', info, string );
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
  ( extensions.drawBuffers ) && rendererExtensions.get( 'WEBGL_draw_buffers' ) ? '#extension GL_EXT_draw_buffers : require' : '',
  ];

  return chunks.filter( filterEmptyLine ).join( '\n' );
}

function generateDefines( defines ) {
  let chunks = [];

  for ( let name in defines ) {
    let value = defines[name];

    if ( value === false ) continue;

    chunks.push( '#define ' + name + ' ' + value );
  }

  return chunks.join( '\n' );
}

function fetchAttributeLocations( gl, program, identifiers ) {
  let attributes = {};

  let n = gl.getProgramParameter( program, gl.ACTIVE_ATTRIBUTES );

  for ( let i = 0; i < n; i ++ ) {
    let info = gl.getActiveAttrib( program, i );
    let name = info.name;

    // console.log("THREE.WebGLProgram: ACTIVE VERTEX ATTRIBUTE:", name, i );

    attributes[name] = gl.getAttribLocation( program, name );
  }

  return attributes;
}

function filterEmptyLine( string ) {
  return string !== '';
}

function parseIncludes( string ) {
  let pattern = /#include +<([\w\d.]+)>/g;

  function replace( match, include ) {
    let replace = ShaderChunk[include];

    if ( replace === undefined ) {
      throw new Error( 'Can not resolve #include <' + include + '>' );
    }

    return parseIncludes( replace );
  }

  return string.replace( pattern, replace );
}

function WebGLProgram( renderer, code, material, parameters ) {
  let gl = renderer.context;

  let extensions = material.extensions;
  let defines = material.defines;

  let vertexShader = material.__webglShader.vertexShader;
  let fragmentShader = material.__webglShader.fragmentShader;

  let gammaFactorDefine = ( renderer.gammaFactor > 0 ) ? renderer.gammaFactor : 1.0;

  // console.log( 'building new program ' );

  //

  let customExtensions = generateExtensions( extensions, parameters, renderer.extensions );

  let customDefines = generateDefines( defines );

  //

  let program = gl.createProgram();

  let prefixVertex, prefixFragment;

  {
    prefixVertex = [

      'precision ' + parameters.precision + ' float;',
      'precision ' + parameters.precision + ' int;',

      '#define SHADER_NAME ' + material.__webglShader.name,

      customDefines,

      '#define GAMMA_FACTOR ' + gammaFactorDefine,

    parameters.vertexColors ? '#define USE_COLOR' : '',

      'uniform mat4 modelMatrix;',
      'uniform mat4 modelViewMatrix;',
      'uniform mat4 projectionMatrix;',
      'uniform mat4 viewMatrix;',
      'uniform mat3 normalMatrix;',
      'uniform vec3 cameraPosition;',

      'attribute vec3 position;',
      'attribute vec3 normal;',
      'attribute vec2 uv;',

      '#ifdef USE_COLOR',

      ' attribute vec3 color;',

      '#endif',
      '\n',
    ].filter( filterEmptyLine ).join( '\n' );

    prefixFragment = [

      customExtensions,

      'precision ' + parameters.precision + ' float;',
      'precision ' + parameters.precision + ' int;',

      '#define SHADER_NAME ' + material.__webglShader.name,

      customDefines,

      '#define GAMMA_FACTOR ' + gammaFactorDefine,

    ( parameters.useFog && parameters.fog ) ? '#define USE_FOG' : '',

    parameters.vertexColors ? '#define USE_COLOR' : '',

    parameters.premultipliedAlpha ? '#define PREMULTIPLIED_ALPHA' : '',

      'uniform mat4 viewMatrix;',
      'uniform vec3 cameraPosition;',

    parameters.depthPacking ? '#define DEPTH_PACKING ' + material.depthPacking : '',

      '\n',

    ].filter( filterEmptyLine ).join( '\n' );
  }

  vertexShader = parseIncludes( vertexShader, parameters );

  fragmentShader = parseIncludes( fragmentShader, parameters );

  let vertexGlsl = prefixVertex + vertexShader;
  let fragmentGlsl = prefixFragment + fragmentShader;

  // console.log( '*VERTEX*', vertexGlsl );
  // console.log( '*FRAGMENT*', fragmentGlsl );

  let glVertexShader = WebGLShader( gl, gl.VERTEX_SHADER, vertexGlsl );
  let glFragmentShader = WebGLShader( gl, gl.FRAGMENT_SHADER, fragmentGlsl );

  gl.attachShader( program, glVertexShader );
  gl.attachShader( program, glFragmentShader );

  // Force a particular attribute to index 0.

  if ( material.index0AttributeName !== undefined ) {
    gl.bindAttribLocation( program, 0, material.index0AttributeName );
  }

  gl.linkProgram( program );

  let programLog = gl.getProgramInfoLog( program );
  let vertexLog = gl.getShaderInfoLog( glVertexShader );
  let fragmentLog = gl.getShaderInfoLog( glFragmentShader );

  let runnable = true;
  let haveDiagnostics = true;

  // console.log( '**VERTEX**', gl.getExtension( 'WEBGL_debug_shaders' ).getTranslatedShaderSource( glVertexShader ) );
  // console.log( '**FRAGMENT**', gl.getExtension( 'WEBGL_debug_shaders' ).getTranslatedShaderSource( glFragmentShader ) );

  if ( gl.getProgramParameter( program, gl.LINK_STATUS ) === false ) {
    runnable = false;

    console.error( 'THREE.WebGLProgram: shader error: ', gl.getError(), 'gl.VALIDATE_STATUS', gl.getProgramParameter( program, gl.VALIDATE_STATUS ), 'gl.getProgramInfoLog', programLog, vertexLog, fragmentLog );
  } else if ( programLog !== '' ) {
    console.warn( 'THREE.WebGLProgram: gl.getProgramInfoLog()', programLog );
  } else if ( vertexLog === '' || fragmentLog === '' ) {
    haveDiagnostics = false;
  }

  if ( haveDiagnostics ) {
    this.diagnostics = {

      runnable: runnable,
      material: material,

      programLog: programLog,

      vertexShader: {

        log: vertexLog,
        prefix: prefixVertex,

      },

      fragmentShader: {

        log: fragmentLog,
        prefix: prefixFragment,

      },

    };
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

  let shaderIDs = {
    LineBasicMaterial: 'basic',
  };

  let parameterNames = [
    'precision',
    'combine', 'vertexColors', 'fog', 'useFog',
    'premultipliedAlpha',
    'depthPacking',
  ];

  this.getParameters = function ( material, fog, object ) {
    let shaderID = shaderIDs[material.type];

    // heuristics to create shader parameters according to lights in the scene
    // (not to blow over maxLights budget)

    let precision = renderer.getPrecision();

    if ( material.precision !== null ) {
      precision = capabilities.getMaxPrecision( material.precision );

      if ( precision !== material.precision ) {
        console.warn( 'THREE.WebGLProgram.getParameters:', material.precision, 'not supported, using', precision, 'instead.' );
      }
    }

    let currentRenderTarget = renderer.getCurrentRenderTarget();

    let parameters = {

      shaderID: shaderID,

      precision: precision,

      combine: material.combine,

      vertexColors: material.vertexColors,

      fog: !! fog,
      useFog: material.fog,

      premultipliedAlpha: material.premultipliedAlpha,

      depthPacking: ( material.depthPacking !== undefined ) ? material.depthPacking : false,

    };

    return parameters;
  };

  this.getProgramCode = function ( material, parameters ) {
    let array = [];

    if ( parameters.shaderID ) {
      array.push( parameters.shaderID );
    } else {
      array.push( material.fragmentShader );
      array.push( material.vertexShader );
    }

    if ( material.defines !== undefined ) {
      for ( let name in material.defines ) {
        array.push( name );
        array.push( material.defines[name] );
      }
    }

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

function WebGLGeometries( gl, properties, info ) {
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

    // TODO

    let property = properties.get( geometry );

    if ( property.wireframe ) {
      deleteAttribute( property.wireframe );
    }

    properties.delete( geometry );

    let bufferproperty = properties.get( buffergeometry );

    if ( bufferproperty.wireframe ) {
      deleteAttribute( bufferproperty.wireframe );
    }

    properties.delete( buffergeometry );

    //

    info.memory.geometries --;
  }

  function getAttributeBuffer( attribute ) {
    if ( attribute.isInterleavedBufferAttribute ) {
      return properties.get( attribute.data ).__webglBuffer;
    }

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
    if ( attribute.isInterleavedBufferAttribute ) {
      properties.delete( attribute.data );
    } else {
      properties.delete( attribute );
    }
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
      } else if ( geometry.isGeometry ) {
        if ( geometry._bufferGeometry === undefined ) {
          geometry._bufferGeometry = new BufferGeometry().setFromObject( object );
        }

        buffergeometry = geometry._bufferGeometry;
      }

      geometries[geometry.id] = buffergeometry;

      info.memory.geometries ++;

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

    if ( object.geometry.isGeometry ) {
      geometry.updateFromObject( object );
    }

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
    let data = ( attribute.isInterleavedBufferAttribute ) ? attribute.data : attribute;

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
      console.error( 'THREE.WebGLObjects.updateBuffer: dynamic THREE.BufferAttribute marked as needsUpdate but updateRange.count is 0, ensure you are using set methods or updating manually.' );
    } else {
      gl.bufferSubData( bufferType, data.updateRange.offset * data.array.BYTES_PER_ELEMENT,
                        data.array.subarray( data.updateRange.offset, data.updateRange.offset + data.updateRange.count ) );

      data.updateRange.count = 0; // reset range
    }

    attributeProperties.version = data.version;
  }

  function getAttributeBuffer( attribute ) {
    if ( attribute.isInterleavedBufferAttribute ) {
      return properties.get( attribute.data ).__webglBuffer;
    }

    return properties.get( attribute ).__webglBuffer;
  }

  function getAttributeProperties( attribute ) {
    if ( attribute.isInterleavedBufferAttribute ) {
      return properties.get( attribute.data );
    }

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

function WebGLTextures( _gl, extensions, state, properties, capabilities, paramThreeToGL, info ) {
  let _infoMemory = info.memory;
  let _isWebGL2 = ( typeof WebGL2RenderingContext !== 'undefined' && _gl instanceof WebGL2RenderingContext );

  //

  function clampToMaxSize( image, maxSize ) {
    if ( image.width > maxSize || image.height > maxSize ) {
      // Warning: Scaling through the canvas will only work with images that use
      // premultiplied alpha.

      let scale = maxSize / Math.max( image.width, image.height );

      let canvas = document.createElementNS( 'http://www.w3.org/1999/xhtml', 'canvas' );
      canvas.width = Math.floor( image.width * scale );
      canvas.height = Math.floor( image.height * scale );

      let context = canvas.getContext( '2d' );
      context.drawImage( image, 0, 0, image.width, image.height, 0, 0, canvas.width, canvas.height );

      console.warn( 'THREE.WebGLRenderer: image is too big (' + image.width + 'x' + image.height + '). Resized to ' + canvas.width + 'x' + canvas.height, image );

      return canvas;
    }

    return image;
  }

  function isPowerOfTwo( image ) {
    return _Math.isPowerOfTwo( image.width ) && _Math.isPowerOfTwo( image.height );
  }

  function makePowerOfTwo( image ) {
    if ( image instanceof HTMLImageElement || image instanceof HTMLCanvasElement ) {
      let canvas = document.createElementNS( 'http://www.w3.org/1999/xhtml', 'canvas' );
      canvas.width = _Math.nearestPowerOfTwo( image.width );
      canvas.height = _Math.nearestPowerOfTwo( image.height );

      let context = canvas.getContext( '2d' );
      context.drawImage( image, 0, 0, canvas.width, canvas.height );

      console.warn( 'THREE.WebGLRenderer: image is not power of two (' + image.width + 'x' + image.height + '). Resized to ' + canvas.width + 'x' + canvas.height, image );

      return canvas;
    }

    return image;
  }

  function textureNeedsPowerOfTwo( texture ) {
    return ( texture.wrapS !== ClampToEdgeWrapping || texture.wrapT !== ClampToEdgeWrapping ) ||
    ( texture.minFilter !== NearestFilter && texture.minFilter !== LinearFilter );
  }

  // Fallback filters for non-power-of-2 textures

  function filterFallback( f ) {
    if ( f === NearestFilter || f === NearestMipMapNearestFilter || f === NearestMipMapLinearFilter ) {
      return _gl.NEAREST;
    }

    return _gl.LINEAR;
  }

  //

  function onTextureDispose( event ) {
    let texture = event.target;

    texture.removeEventListener( 'dispose', onTextureDispose );

    deallocateTexture( texture );

    _infoMemory.textures --;
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
        console.warn( 'THREE.WebGLRenderer: Texture marked for update but image is undefined', texture );
      } else if ( image.complete === false ) {
        console.warn( 'THREE.WebGLRenderer: Texture marked for update but image is incomplete', texture );
      } else {
        uploadTexture( textureProperties, texture, slot );
        return;
      }
    }

    state.activeTexture( _gl.TEXTURE0 + slot );
    state.bindTexture( _gl.TEXTURE_2D, textureProperties.__webglTexture );
  }

  function setTextureParameters( textureType, texture, isPowerOfTwoImage ) {
    let extension;

    if ( isPowerOfTwoImage ) {
      _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_S, paramThreeToGL( texture.wrapS ) );
      _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_T, paramThreeToGL( texture.wrapT ) );

      _gl.texParameteri( textureType, _gl.TEXTURE_MAG_FILTER, paramThreeToGL( texture.magFilter ) );
      _gl.texParameteri( textureType, _gl.TEXTURE_MIN_FILTER, paramThreeToGL( texture.minFilter ) );
    } else {
      _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_S, _gl.CLAMP_TO_EDGE );
      _gl.texParameteri( textureType, _gl.TEXTURE_WRAP_T, _gl.CLAMP_TO_EDGE );

      if ( texture.wrapS !== ClampToEdgeWrapping || texture.wrapT !== ClampToEdgeWrapping ) {
        console.warn( 'THREE.WebGLRenderer: Texture is not power of two. Texture.wrapS and Texture.wrapT should be set to THREE.ClampToEdgeWrapping.', texture );
      }

      _gl.texParameteri( textureType, _gl.TEXTURE_MAG_FILTER, filterFallback( texture.magFilter ) );
      _gl.texParameteri( textureType, _gl.TEXTURE_MIN_FILTER, filterFallback( texture.minFilter ) );

      if ( texture.minFilter !== NearestFilter && texture.minFilter !== LinearFilter ) {
        console.warn( 'THREE.WebGLRenderer: Texture is not power of two. Texture.minFilter should be set to THREE.NearestFilter or THREE.LinearFilter.', texture );
      }
    }

    extension = extensions.get( 'EXT_texture_filter_anisotropic' );

    if ( extension ) {
      if ( texture.type === FloatType && extensions.get( 'OES_texture_float_linear' ) === null ) return;
      if ( texture.type === HalfFloatType && extensions.get( 'OES_texture_half_float_linear' ) === null ) return;

      if ( texture.anisotropy > 1 || properties.get( texture ).__currentAnisotropy ) {
        _gl.texParameterf( textureType, extension.TEXTURE_MAX_ANISOTROPY_EXT, Math.min( texture.anisotropy, capabilities.getMaxAnisotropy() ) );
        properties.get( texture ).__currentAnisotropy = texture.anisotropy;
      }
    }
  }

  function uploadTexture( textureProperties, texture, slot ) {
    if ( textureProperties.__webglInit === undefined ) {
      textureProperties.__webglInit = true;

      texture.addEventListener( 'dispose', onTextureDispose );

      textureProperties.__webglTexture = _gl.createTexture();

      _infoMemory.textures ++;
    }

    state.activeTexture( _gl.TEXTURE0 + slot );
    state.bindTexture( _gl.TEXTURE_2D, textureProperties.__webglTexture );

    _gl.pixelStorei( _gl.UNPACK_FLIP_Y_WEBGL, texture.flipY );
    _gl.pixelStorei( _gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, texture.premultiplyAlpha );
    _gl.pixelStorei( _gl.UNPACK_ALIGNMENT, texture.unpackAlignment );

    let image = clampToMaxSize( texture.image, capabilities.maxTextureSize );

    if ( textureNeedsPowerOfTwo( texture ) && isPowerOfTwo( image ) === false ) {
      image = makePowerOfTwo( image );
    }

    let isPowerOfTwoImage = isPowerOfTwo( image ),
      glFormat = paramThreeToGL( texture.format ),
      glType = paramThreeToGL( texture.type );

    setTextureParameters( _gl.TEXTURE_2D, texture, isPowerOfTwoImage );

    let mipmap, mipmaps = texture.mipmaps;

    {
      // regular Texture (image, video, canvas)

      // use manually created mipmaps if available
      // if there are no manual mipmaps
      // set 0 level mipmap and then use GL to generate other mipmap levels

      if ( mipmaps.length > 0 && isPowerOfTwoImage ) {
        for ( let i = 0, il = mipmaps.length; i < il; i ++ ) {
          mipmap = mipmaps[i];
          state.texImage2D( _gl.TEXTURE_2D, i, glFormat, glFormat, glType, mipmap );
        }

        texture.generateMipmaps = false;
      } else {
        state.texImage2D( _gl.TEXTURE_2D, 0, glFormat, glFormat, glType, image );
      }
    }

    if ( texture.generateMipmaps && isPowerOfTwoImage ) _gl.generateMipmap( _gl.TEXTURE_2D );

    textureProperties.__version = texture.version;

    if ( texture.onUpdate ) texture.onUpdate( texture );
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

function WebGLState( gl, extensions, paramThreeToGL ) {
  function ColorBuffer() {
    let locked = false;

    let color = new Vector4();
    let currentColorMask = null;
    let currentColorClear = new Vector4();

    return {

      setMask: function ( colorMask ) {
        if ( currentColorMask !== colorMask && ! locked ) {
          gl.colorMask( colorMask, colorMask, colorMask, colorMask );
          currentColorMask = colorMask;
        }
      },

      setLocked: function ( lock ) {
        locked = lock;
      },

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
    let locked = false;

    let currentDepthMask = null;
    let currentDepthFunc = null;
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
        if ( currentDepthMask !== depthMask && ! locked ) {
          gl.depthMask( depthMask );
          currentDepthMask = depthMask;
        }
      },

      setFunc: function ( depthFunc ) {
        if ( currentDepthFunc !== depthFunc ) {
          if ( depthFunc ) {
            switch ( depthFunc ) {
              case NeverDepth:

                gl.depthFunc( gl.NEVER );
                break;

              case AlwaysDepth:

                gl.depthFunc( gl.ALWAYS );
                break;

              case LessDepth:

                gl.depthFunc( gl.LESS );
                break;

              case LessEqualDepth:

                gl.depthFunc( gl.LEQUAL );
                break;

              case EqualDepth:

                gl.depthFunc( gl.EQUAL );
                break;

              case GreaterEqualDepth:

                gl.depthFunc( gl.GEQUAL );
                break;

              case GreaterDepth:

                gl.depthFunc( gl.GREATER );
                break;

              case NotEqualDepth:

                gl.depthFunc( gl.NOTEQUAL );
                break;

              default:

                gl.depthFunc( gl.LEQUAL );
            }
          } else {
            gl.depthFunc( gl.LEQUAL );
          }

          currentDepthFunc = depthFunc;
        }
      },

      setLocked: function ( lock ) {
        locked = lock;
      },

      setClear: function ( depth ) {
        if ( currentDepthClear !== depth ) {
          gl.clearDepth( depth );
          currentDepthClear = depth;
        }
      },

      reset: function () {
        locked = false;

        currentDepthMask = null;
        currentDepthFunc = null;
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
  let attributeDivisors = new Uint8Array( maxVertexAttributes );

  let capabilities = {};

  let compressedTextureFormats = null;

  let currentBlending = null;
  let currentBlendEquation = null;
  let currentBlendSrc = null;
  let currentBlendDst = null;
  let currentBlendEquationAlpha = null;
  let currentBlendSrcAlpha = null;
  let currentBlendDstAlpha = null;
  let currentPremultipledAlpha = false;

  let currentLineWidth = null;

  let currentPolygonOffsetFactor = null;
  let currentPolygonOffsetUnits = null;

  let currentScissorTest = null;

  let maxTextures = gl.getParameter( gl.MAX_TEXTURE_IMAGE_UNITS );

  let version = parseFloat( /^WebGL\ ([0-9])/.exec( gl.getParameter( gl.VERSION ) )[1] );
  let lineWidthAvailable = parseFloat( version ) >= 1.0;

  let currentTextureSlot = null;
  let currentBoundTextures = {};

  let currentScissor = new Vector4();
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
    setDepthFunc( LessEqualDepth );

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

    if ( attributeDivisors[attribute] !== 0 ) {
      let extension = extensions.get( 'ANGLE_instanced_arrays' );

      extension.vertexAttribDivisorANGLE( attribute, 0 );
      attributeDivisors[attribute] = 0;
    }
  }

  function enableAttributeAndDivisor( attribute, meshPerAttribute, extension ) {
    newAttributes[attribute] = 1;

    if ( enabledAttributes[attribute] === 0 ) {
      gl.enableVertexAttribArray( attribute );
      enabledAttributes[attribute] = 1;
    }

    if ( attributeDivisors[attribute] !== meshPerAttribute ) {
      extension.vertexAttribDivisorANGLE( attribute, meshPerAttribute );
      attributeDivisors[attribute] = meshPerAttribute;
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

  function setBlending( blending, blendEquation, blendSrc, blendDst, blendEquationAlpha, blendSrcAlpha, blendDstAlpha, premultipliedAlpha ) {
    if ( blending !== NoBlending ) {
      enable( gl.BLEND );
    } else {
      disable( gl.BLEND );
    }

    if ( blending !== currentBlending || premultipliedAlpha !== currentPremultipledAlpha ) {
      if ( blending === AdditiveBlending ) {
        if ( premultipliedAlpha ) {
          gl.blendEquationSeparate( gl.FUNC_ADD, gl.FUNC_ADD );
          gl.blendFuncSeparate( gl.ONE, gl.ONE, gl.ONE, gl.ONE );
        } else {
          gl.blendEquation( gl.FUNC_ADD );
          gl.blendFunc( gl.SRC_ALPHA, gl.ONE );
        }
      } else if ( blending === SubtractiveBlending ) {
        if ( premultipliedAlpha ) {
          gl.blendEquationSeparate( gl.FUNC_ADD, gl.FUNC_ADD );
          gl.blendFuncSeparate( gl.ZERO, gl.ZERO, gl.ONE_MINUS_SRC_COLOR, gl.ONE_MINUS_SRC_ALPHA );
        } else {
          gl.blendEquation( gl.FUNC_ADD );
          gl.blendFunc( gl.ZERO, gl.ONE_MINUS_SRC_COLOR );
        }
      } else if ( blending === MultiplyBlending ) {
        if ( premultipliedAlpha ) {
          gl.blendEquationSeparate( gl.FUNC_ADD, gl.FUNC_ADD );
          gl.blendFuncSeparate( gl.ZERO, gl.SRC_COLOR, gl.ZERO, gl.SRC_ALPHA );
        } else {
          gl.blendEquation( gl.FUNC_ADD );
          gl.blendFunc( gl.ZERO, gl.SRC_COLOR );
        }
      } else {
        if ( premultipliedAlpha ) {
          gl.blendEquationSeparate( gl.FUNC_ADD, gl.FUNC_ADD );
          gl.blendFuncSeparate( gl.ONE, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA );
        } else {
          gl.blendEquationSeparate( gl.FUNC_ADD, gl.FUNC_ADD );
          gl.blendFuncSeparate( gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA );
        }
      }

      currentBlending = blending;
      currentPremultipledAlpha = premultipliedAlpha;
    }

    if ( blending === CustomBlending ) {
      blendEquationAlpha = blendEquationAlpha || blendEquation;
      blendSrcAlpha = blendSrcAlpha || blendSrc;
      blendDstAlpha = blendDstAlpha || blendDst;

      if ( blendEquation !== currentBlendEquation || blendEquationAlpha !== currentBlendEquationAlpha ) {
        gl.blendEquationSeparate( paramThreeToGL( blendEquation ), paramThreeToGL( blendEquationAlpha ) );

        currentBlendEquation = blendEquation;
        currentBlendEquationAlpha = blendEquationAlpha;
      }

      if ( blendSrc !== currentBlendSrc || blendDst !== currentBlendDst || blendSrcAlpha !== currentBlendSrcAlpha || blendDstAlpha !== currentBlendDstAlpha ) {
        gl.blendFuncSeparate( paramThreeToGL( blendSrc ), paramThreeToGL( blendDst ), paramThreeToGL( blendSrcAlpha ), paramThreeToGL( blendDstAlpha ) );

        currentBlendSrc = blendSrc;
        currentBlendDst = blendDst;
        currentBlendSrcAlpha = blendSrcAlpha;
        currentBlendDstAlpha = blendDstAlpha;
      }
    } else {
      currentBlendEquation = null;
      currentBlendSrc = null;
      currentBlendDst = null;
      currentBlendEquationAlpha = null;
      currentBlendSrcAlpha = null;
      currentBlendDstAlpha = null;
    }
  }

  // TODO Deprecate

  function setColorWrite( colorWrite ) {
    colorBuffer.setMask( colorWrite );
  }

  function setDepthTest( depthTest ) {
    depthBuffer.setTest( depthTest );
  }

  function setDepthWrite( depthWrite ) {
    depthBuffer.setMask( depthWrite );
  }

  function setDepthFunc( depthFunc ) {
    depthBuffer.setFunc( depthFunc );
  }

  //

  function setLineWidth( width ) {
    if ( width !== currentLineWidth ) {
      if ( lineWidthAvailable ) gl.lineWidth( width );

      currentLineWidth = width;
    }
  }

  function setPolygonOffset( polygonOffset, factor, units ) {
    if ( polygonOffset ) {
      enable( gl.POLYGON_OFFSET_FILL );

      if ( currentPolygonOffsetFactor !== factor || currentPolygonOffsetUnits !== units ) {
        gl.polygonOffset( factor, units );

        currentPolygonOffsetFactor = factor;
        currentPolygonOffsetUnits = units;
      }
    } else {
      disable( gl.POLYGON_OFFSET_FILL );
    }
  }

  function getScissorTest() {
    return currentScissorTest;
  }

  function setScissorTest( scissorTest ) {
    currentScissorTest = scissorTest;

    if ( scissorTest ) {
      enable( gl.SCISSOR_TEST );
    } else {
      disable( gl.SCISSOR_TEST );
    }
  }

  // texture

  function activeTexture( webglSlot ) {
    if ( webglSlot === undefined ) webglSlot = gl.TEXTURE0 + maxTextures - 1;

    if ( currentTextureSlot !== webglSlot ) {
      gl.activeTexture( webglSlot );
      currentTextureSlot = webglSlot;
    }
  }

  function bindTexture( webglType, webglTexture ) {
    if ( currentTextureSlot === null ) {
      activeTexture();
    }

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

  function scissor( scissor ) {
    if ( currentScissor.equals( scissor ) === false ) {
      gl.scissor( scissor.x, scissor.y, scissor.z, scissor.w );
      currentScissor.copy( scissor );
    }
  }

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
    enableAttributeAndDivisor: enableAttributeAndDivisor,
    disableUnusedAttributes: disableUnusedAttributes,
    enable: enable,
    disable: disable,

    setBlending: setBlending,

    setColorWrite: setColorWrite,
    setDepthTest: setDepthTest,
    setDepthWrite: setDepthWrite,
    setDepthFunc: setDepthFunc,

    setLineWidth: setLineWidth,
    setPolygonOffset: setPolygonOffset,

    getScissorTest: getScissorTest,
    setScissorTest: setScissorTest,

    activeTexture: activeTexture,
    bindTexture: bindTexture,
    texImage2D: texImage2D,

    scissor: scissor,
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
    console.warn( 'THREE.WebGLRenderer:', precision, 'not supported, using', maxPrecision, 'instead.' );
    precision = maxPrecision;
  }

  let maxTextures = gl.getParameter( gl.MAX_TEXTURE_IMAGE_UNITS );
  let maxVertexTextures = gl.getParameter( gl.MAX_VERTEX_TEXTURE_IMAGE_UNITS );
  let maxTextureSize = gl.getParameter( gl.MAX_TEXTURE_SIZE );

  let maxAttributes = gl.getParameter( gl.MAX_VERTEX_ATTRIBS );
  let maxVertexUniforms = gl.getParameter( gl.MAX_VERTEX_UNIFORM_VECTORS );
  let maxVaryings = gl.getParameter( gl.MAX_VARYING_VECTORS );
  let maxFragmentUniforms = gl.getParameter( gl.MAX_FRAGMENT_UNIFORM_VECTORS );

  let vertexTextures = maxVertexTextures > 0;
  let floatFragmentTextures = !! extensions.get( 'OES_texture_float' );
  let floatVertexTextures = vertexTextures && floatFragmentTextures;

  return {

    getMaxPrecision: getMaxPrecision,

    precision: precision,

    maxTextures: maxTextures,
    maxVertexTextures: maxVertexTextures,
    maxTextureSize: maxTextureSize,

    maxAttributes: maxAttributes,
    maxVertexUniforms: maxVertexUniforms,
    maxVaryings: maxVaryings,
    maxFragmentUniforms: maxFragmentUniforms,

    vertexTextures: vertexTextures,
    floatFragmentTextures: floatFragmentTextures,
    floatVertexTextures: floatVertexTextures,

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

        case 'EXT_texture_filter_anisotropic':
          extension = gl.getExtension( 'EXT_texture_filter_anisotropic' ) || gl.getExtension( 'MOZ_EXT_texture_filter_anisotropic' ) || gl.getExtension( 'WEBKIT_EXT_texture_filter_anisotropic' );
          break;

        default:
          extension = gl.getExtension( name );
      }

      if ( extension === null ) {
        console.warn( 'THREE.WebGLRenderer: ' + name + ' extension not supported.' );
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
  console.log( 'THREE.WebGLRenderer', REVISION );

  parameters = parameters || {};

  let _canvas = parameters.canvas !== undefined ? parameters.canvas : document.createElementNS( 'http://www.w3.org/1999/xhtml', 'canvas' ),
    _context = parameters.context !== undefined ? parameters.context : null,

    _alpha = parameters.alpha !== undefined ? parameters.alpha : false,
    _depth = parameters.depth !== undefined ? parameters.depth : true,
    _antialias = parameters.antialias !== undefined ? parameters.antialias : false,
    _premultipliedAlpha = parameters.premultipliedAlpha !== undefined ? parameters.premultipliedAlpha : true,
    _preserveDrawingBuffer = parameters.preserveDrawingBuffer !== undefined ? parameters.preserveDrawingBuffer : false;

  let lights = [];

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

  // physically based shading

  this.gammaFactor = 2.0; // for backwards compatibility
  this.gammaInput = false;
  this.gammaOutput = false;

  // internal properties

  let _this = this,

    // internal state cache

    _currentProgram = null,
    _currentRenderTarget = null,
    _currentFramebuffer = null,
    _currentMaterialId = - 1,
    _currentGeometryProgram = '',
    _currentCamera = null,

    _currentScissor = new Vector4(),
    _currentScissorTest = null,

    _currentViewport = new Vector4(),

    //

    _usedTextureUnits = 0,

    //

    _clearColor = new Color( 0x000000 ),
    _clearAlpha = 0,

    _width = _canvas.width,
    _height = _canvas.height,

    _pixelRatio = 1,

    _scissor = new Vector4( 0, 0, _width, _height ),
    _scissorTest = false,

    _viewport = new Vector4( 0, 0, _width, _height ),

    // camera matrices cache

    _projScreenMatrix = new Matrix4(),

    _vector3 = new Vector3(),

    // light arrays cache

    _ambient_light = [0, 0, 0],

    // info

    _infoRender = {

      calls: 0,
      vertices: 0,
      faces: 0,
      points: 0,

    };

  this.info = {

    render: _infoRender,
    memory: {

      geometries: 0,
      textures: 0,

    },
    programs: null,

  };


  // initialize

  let _gl;

  try {
    let attributes = {
      alpha: _alpha,
      depth: _depth,
      antialias: _antialias,
      premultipliedAlpha: _premultipliedAlpha,
      preserveDrawingBuffer: _preserveDrawingBuffer,
    };

    _gl = _context || _canvas.getContext( 'webgl', attributes ) || _canvas.getContext( 'experimental-webgl', attributes );

    if ( _gl === null ) {
      if ( _canvas.getContext( 'webgl' ) !== null ) {
        throw 'Error creating WebGL context with your selected attributes.';
      } else {
        throw 'Error creating WebGL context.';
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
    console.error( 'THREE.WebGLRenderer: ' + error );
  }

  let extensions = new WebGLExtensions( _gl );

  extensions.get( 'WEBGL_depth_texture' );
  extensions.get( 'OES_texture_float' );
  extensions.get( 'OES_texture_float_linear' );
  extensions.get( 'OES_texture_half_float' );
  extensions.get( 'OES_texture_half_float_linear' );
  extensions.get( 'OES_standard_derivatives' );
  extensions.get( 'ANGLE_instanced_arrays' );

  if ( extensions.get( 'OES_element_index_uint' ) ) {
    BufferGeometry.MaxIndex = 4294967296;
  }

  let capabilities = new WebGLCapabilities( _gl, extensions, parameters );

  let state = new WebGLState( _gl, extensions, paramThreeToGL );
  let properties = new WebGLProperties();
  let textures = new WebGLTextures( _gl, extensions, state, properties, capabilities, paramThreeToGL, this.info );
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

    state.scissor( _currentScissor.copy( _scissor ).multiplyScalar( _pixelRatio ) );
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

  this.getContextAttributes = function () {
    return _gl.getContextAttributes();
  };

  this.forceContextLoss = function () {
    extensions.get( 'WEBGL_lose_context' ).loseContext();
  };

  this.getMaxAnisotropy = function () {
    return capabilities.getMaxAnisotropy();
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

  this.setScissor = function ( x, y, width, height ) {
    state.scissor( _scissor.set( x, y, width, height ) );
  };

  this.setScissorTest = function ( boolean ) {
    state.setScissorTest( _scissorTest = boolean );
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

  this.getClearAlpha = function () {
    return _clearAlpha;
  };

  this.setClearAlpha = function ( alpha ) {
    _clearAlpha = alpha;

    state.buffers.color.setClear( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha, _premultipliedAlpha );
  };

  this.clear = function ( color, depth ) {
    let bits = 0;

    if ( color === undefined || color ) bits |= _gl.COLOR_BUFFER_BIT;
    if ( depth === undefined || depth ) bits |= _gl.DEPTH_BUFFER_BIT;

    _gl.clear( bits );
  };

  this.clearColor = function () {
    this.clear( true, false, false );
  };

  this.clearDepth = function () {
    this.clear( false, true, false );
  };

  this.clearTarget = function ( renderTarget, color, depth ) {
    this.setRenderTarget( renderTarget );
    this.clear( color, depth );
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
    let geometryProgram = geometry.id + '_' + program.id + '_' + material.wireframe;

    if ( geometryProgram !== _currentGeometryProgram ) {
      _currentGeometryProgram = geometryProgram;
      updateBuffers = true;
    }

    //

    let index = geometry.index;
    let position = geometry.attributes.position;
    let rangeFactor = 1;

    if ( material.wireframe === true ) {
      index = objects.getWireframeAttribute( geometry );
      rangeFactor = 2;
    }

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
      if ( material.wireframe === true ) {
        state.setLineWidth( material.wireframeLinewidth * getTargetPixelRatio() );
        renderer.setMode( _gl.LINES );
      } else {
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
    let extension;

    if ( startIndex === undefined ) startIndex = 0;

    state.initAttributes();

    let geometryAttributes = geometry.attributes;

    let programAttributes = program.getAttributes();

    let materialDefaultAttributeValues = material.defaultAttributeValues;

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

          if ( geometryAttribute.isInterleavedBufferAttribute ) {
            let data = geometryAttribute.data;
            let stride = data.stride;
            let offset = geometryAttribute.offset;

            if ( data && data.isInstancedInterleavedBuffer ) {
              state.enableAttributeAndDivisor( programAttribute, data.meshPerAttribute, extension );

              if ( geometry.maxInstancedCount === undefined ) {
                geometry.maxInstancedCount = data.meshPerAttribute * data.count;
              }
            } else {
              state.enableAttribute( programAttribute );
            }

            _gl.bindBuffer( _gl.ARRAY_BUFFER, buffer );
            _gl.vertexAttribPointer( programAttribute, size, type, normalized, stride * bytesPerElement, ( startIndex * stride + offset ) * bytesPerElement );
          } else {
            if ( geometryAttribute.isInstancedBufferAttribute ) {
              state.enableAttributeAndDivisor( programAttribute, geometryAttribute.meshPerAttribute, extension );

              if ( geometry.maxInstancedCount === undefined ) {
                geometry.maxInstancedCount = geometryAttribute.meshPerAttribute * geometryAttribute.count;
              }
            } else {
              state.enableAttribute( programAttribute );
            }

            _gl.bindBuffer( _gl.ARRAY_BUFFER, buffer );
            _gl.vertexAttribPointer( programAttribute, size, type, normalized, 0, startIndex * size * bytesPerElement );
          }
        } else if ( materialDefaultAttributeValues !== undefined ) {
          let value = materialDefaultAttributeValues[name];

          if ( value !== undefined ) {
            switch ( value.length ) {
              case 2:
                _gl.vertexAttrib2fv( programAttribute, value );
                break;

              case 3:
                _gl.vertexAttrib3fv( programAttribute, value );
                break;

              case 4:
                _gl.vertexAttrib4fv( programAttribute, value );
                break;

              default:
                _gl.vertexAttrib1fv( programAttribute, value );
            }
          }
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
      console.error( 'THREE.WebGLRenderer.render: camera is not an instance of THREE.Camera.' );
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

    camera.matrixWorldInverse.getInverse( camera.matrixWorld );

    _projScreenMatrix.multiplyMatrices( camera.projectionMatrix, camera.matrixWorldInverse );

    lights.length = 0;

    opaqueObjectsLastIndex = - 1;
    transparentObjectsLastIndex = - 1;

    projectObject( scene, camera );

    opaqueObjects.length = opaqueObjectsLastIndex + 1;
    transparentObjects.length = transparentObjectsLastIndex + 1;

    if ( _this.sortObjects === true ) {
      opaqueObjects.sort( painterSortStable );
      transparentObjects.sort( reversePainterSortStable );
    }

    //
    if ( lights[0] && lights[0].isAmbientLight ) {
      _ambient_light[0] = lights[0].color.r;
      _ambient_light[1] = lights[0].color.g;
      _ambient_light[2] = lights[0].color.b;
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

    //

    let background = scene.background;

    if ( background === null ) {
      state.buffers.color.setClear( _clearColor.r, _clearColor.g, _clearColor.b, _clearAlpha, _premultipliedAlpha );
    } else if ( background && background.isColor ) {
      state.buffers.color.setClear( background.r, background.g, background.b, 1, _premultipliedAlpha );
      forceClear = true;
    }

    if ( this.autoClear || forceClear ) {
      this.clear( this.autoClearColor, this.autoClearDepth );
    }

    if ( scene.overrideMaterial ) {
      let overrideMaterial = scene.overrideMaterial;

      renderObjects( opaqueObjects, scene, camera, overrideMaterial );
      renderObjects( transparentObjects, scene, camera, overrideMaterial );
    } else {
      // opaque pass (front-to-back order)

      state.setBlending( NoBlending );
      renderObjects( opaqueObjects, scene, camera );

      // transparent pass (back-to-front order)

      renderObjects( transparentObjects, scene, camera );
    }

    // Generate mipmap if we're using any kind of mipmap filtering

    if ( renderTarget ) {
      textures.updateRenderTargetMipmap( renderTarget );
    }

    // Ensure depth buffer writing is enabled so it can be cleared on next render

    state.setDepthTest( true );
    state.setDepthWrite( true );
    state.setColorWrite( true );

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

  function projectObject( object, camera ) {
    if ( object.visible === false ) return;

    if ( object.isAmbientLight ) {
      lights.push( object );
    } else if ( object.isMesh || object.isLine || object.isPoints ) {
      let material = object.material;

      if ( material.visible === true ) {
        if ( _this.sortObjects === true ) {
          _vector3.setFromMatrixPosition( object.matrixWorld );
          _vector3.applyProjection( _projScreenMatrix );
        }

        let geometry = objects.update( object );

        pushRenderItem( object, geometry, material, _vector3.z, null );
      }
    }

    let children = object.children;

    for ( let i = 0, l = children.length; i < l; i ++ ) {
      projectObject( children[i], camera );
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
      object.normalMatrix.getNormalMatrix( object.modelViewMatrix );

      object.onBeforeRender( _this, scene, camera, geometry, material, group );

      _this.renderBufferDirect( camera, scene.fog, geometry, material, object, group );

      object.onAfterRender( _this, scene, camera, geometry, material, group );
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
    } else if ( parameters.shaderID !== undefined ) {
      // same glsl and uniform list
      return;
    } else {
      // only rebuild uniform list
      programChange = false;
    }

    if ( programChange ) {
      if ( parameters.shaderID ) {
        let shader = ShaderLib[parameters.shaderID];

        materialProperties.__webglShader = {
          name: material.type,
          uniforms: UniformsUtils.clone( shader.uniforms ),
          vertexShader: shader.vertexShader,
          fragmentShader: shader.fragmentShader,
        };
      } else {
        materialProperties.__webglShader = {
          name: material.type,
          uniforms: material.uniforms,
          vertexShader: material.vertexShader,
          fragmentShader: material.fragmentShader,
        };
      }

      material.__webglShader = materialProperties.__webglShader;

      program = programCache.acquireProgram( material, parameters, code );

      materialProperties.program = program;
      material.program = program;
    }

    let attributes = program.getAttributes();

    let uniforms = materialProperties.__webglShader.uniforms;

    materialProperties.fog = fog;

    if ( material.lights ) {
      // wire up the material to this renderer's lighting state

      uniforms.ambientLightColor.value = _ambient_light;
    }

    let progUniforms = materialProperties.program.getUniforms(),
      uniformsList =
      WebGLUniforms.seqWithValue( progUniforms.seq, uniforms );

    materialProperties.uniformsList = uniformsList;
  }

  function setMaterial( material ) {
  material.transparent === true ?
    state.setBlending( material.blending, material.blendEquation, material.blendSrc, material.blendDst, material.blendEquationAlpha, material.blendSrcAlpha, material.blendDstAlpha, material.premultipliedAlpha )
    : state.setBlending( NoBlending );

  state.setDepthFunc( material.depthFunc );
  state.setDepthTest( material.depthTest );
  state.setDepthWrite( material.depthWrite );
  state.setColorWrite( material.colorWrite );
  state.setPolygonOffset( material.polygonOffset, material.polygonOffsetFactor, material.polygonOffsetUnits );
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
        let uCamPos = p_uniforms.map.cameraPosition;

        if ( uCamPos !== undefined ) {
          uCamPos.setValue( _gl,
                            _vector3.setFromMatrixPosition( camera.matrixWorld ) );
        }
      }

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

      if ( material.isLineBasicMaterial ) {
        refreshUniformsLine( m_uniforms, material );
      }

      WebGLUniforms.upload(
        _gl, materialProperties.uniformsList, m_uniforms, _this );
    }


    // common matrices

    p_uniforms.set( _gl, object, 'modelViewMatrix' );
    p_uniforms.set( _gl, object, 'normalMatrix' );
    p_uniforms.setValue( _gl, 'modelMatrix', object.matrixWorld );

    return program;
  }

  // Uniforms (refresh uniforms objects)

  function refreshUniformsCommon( uniforms, material ) {
    uniforms.opacity.value = material.opacity;
    uniforms.diffuse.value = material.color;
  }

  function refreshUniformsLine( uniforms, material ) {
    uniforms.diffuse.value = material.color;
    uniforms.opacity.value = material.opacity;
  }

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

    if ( textureUnit >= capabilities.maxTextures ) {
      console.warn( 'WebGLRenderer: trying to use ' + textureUnit + ' texture units while this GPU supports only ' + capabilities.maxTextures );
    }

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
          console.warn( 'THREE.WebGLRenderer.setTexture2D: don\'t use render targets as textures. Use their .texture property instead.' );
          warned = true;
        }

        texture = texture.texture;
      }

      textures.setTexture2D( texture, slot );
    };
  }() );

  this.setTexture = ( function () {
    let warned = false;

    return function setTexture( texture, slot ) {
      if ( ! warned ) {
        console.warn( 'THREE.WebGLRenderer: .setTexture is deprecated, use setTexture2D instead.' );
        warned = true;
      }

      textures.setTexture2D( texture, slot );
    };
  }() );

  this.getCurrentRenderTarget = function () {
    return _currentRenderTarget;
  };

  this.setRenderTarget = function ( renderTarget ) {
    _currentRenderTarget = renderTarget;

    if ( renderTarget && properties.get( renderTarget ).__webglFramebuffer === undefined ) {
      textures.setupRenderTarget( renderTarget );
    }

    let framebuffer;

    if ( renderTarget ) {
      let renderTargetProperties = properties.get( renderTarget );

      framebuffer = renderTargetProperties.__webglFramebuffer;

      _currentScissor.copy( renderTarget.scissor );
      _currentScissorTest = renderTarget.scissorTest;

      _currentViewport.copy( renderTarget.viewport );
    } else {
      framebuffer = null;

      _currentScissor.copy( _scissor ).multiplyScalar( _pixelRatio );
      _currentScissorTest = _scissorTest;

      _currentViewport.copy( _viewport ).multiplyScalar( _pixelRatio );
    }

    if ( _currentFramebuffer !== framebuffer ) {
      _gl.bindFramebuffer( _gl.FRAMEBUFFER, framebuffer );
      _currentFramebuffer = framebuffer;
    }

    state.scissor( _currentScissor );
    state.setScissorTest( _currentScissorTest );

    state.viewport( _currentViewport );
  };


  // Map three.js constants to WebGL constants

  function paramThreeToGL( p ) {
    let extension;

    if ( p === RepeatWrapping ) return _gl.REPEAT;
    if ( p === ClampToEdgeWrapping ) return _gl.CLAMP_TO_EDGE;
    if ( p === MirroredRepeatWrapping ) return _gl.MIRRORED_REPEAT;

    if ( p === NearestFilter ) return _gl.NEAREST;
    if ( p === NearestMipMapNearestFilter ) return _gl.NEAREST_MIPMAP_NEAREST;
    if ( p === NearestMipMapLinearFilter ) return _gl.NEAREST_MIPMAP_LINEAR;

    if ( p === LinearFilter ) return _gl.LINEAR;
    if ( p === LinearMipMapNearestFilter ) return _gl.LINEAR_MIPMAP_NEAREST;
    if ( p === LinearMipMapLinearFilter ) return _gl.LINEAR_MIPMAP_LINEAR;

    if ( p === UnsignedByteType ) return _gl.UNSIGNED_BYTE;
    if ( p === UnsignedShort4444Type ) return _gl.UNSIGNED_SHORT_4_4_4_4;
    if ( p === UnsignedShort5551Type ) return _gl.UNSIGNED_SHORT_5_5_5_1;
    if ( p === UnsignedShort565Type ) return _gl.UNSIGNED_SHORT_5_6_5;

    if ( p === ByteType ) return _gl.BYTE;
    if ( p === ShortType ) return _gl.SHORT;
    if ( p === UnsignedShortType ) return _gl.UNSIGNED_SHORT;
    if ( p === IntType ) return _gl.INT;
    if ( p === UnsignedIntType ) return _gl.UNSIGNED_INT;
    if ( p === FloatType ) return _gl.FLOAT;

    if ( p === HalfFloatType ) {
      extension = extensions.get( 'OES_texture_half_float' );

      if ( extension !== null ) return extension.HALF_FLOAT_OES;
    }

    if ( p === AlphaFormat ) return _gl.ALPHA;
    if ( p === RGBFormat ) return _gl.RGB;
    if ( p === RGBAFormat ) return _gl.RGBA;
    if ( p === LuminanceFormat ) return _gl.LUMINANCE;
    if ( p === LuminanceAlphaFormat ) return _gl.LUMINANCE_ALPHA;
    if ( p === DepthFormat ) return _gl.DEPTH_COMPONENT;
    if ( p === DepthStencilFormat ) return _gl.DEPTH_STENCIL;

    if ( p === AddEquation ) return _gl.FUNC_ADD;
    if ( p === SubtractEquation ) return _gl.FUNC_SUBTRACT;
    if ( p === ReverseSubtractEquation ) return _gl.FUNC_REVERSE_SUBTRACT;

    if ( p === ZeroFactor ) return _gl.ZERO;
    if ( p === OneFactor ) return _gl.ONE;
    if ( p === SrcColorFactor ) return _gl.SRC_COLOR;
    if ( p === OneMinusSrcColorFactor ) return _gl.ONE_MINUS_SRC_COLOR;
    if ( p === SrcAlphaFactor ) return _gl.SRC_ALPHA;
    if ( p === OneMinusSrcAlphaFactor ) return _gl.ONE_MINUS_SRC_ALPHA;
    if ( p === DstAlphaFactor ) return _gl.DST_ALPHA;
    if ( p === OneMinusDstAlphaFactor ) return _gl.ONE_MINUS_DST_ALPHA;

    if ( p === DstColorFactor ) return _gl.DST_COLOR;
    if ( p === OneMinusDstColorFactor ) return _gl.ONE_MINUS_DST_COLOR;
    if ( p === SrcAlphaSaturateFactor ) return _gl.SRC_ALPHA_SATURATE;

    if ( p === RGB_S3TC_DXT1_Format || p === RGBA_S3TC_DXT1_Format ||
    p === RGBA_S3TC_DXT3_Format || p === RGBA_S3TC_DXT5_Format ) {
      extension = extensions.get( 'WEBGL_compressed_texture_s3tc' );

      if ( extension !== null ) {
        if ( p === RGB_S3TC_DXT1_Format ) return extension.COMPRESSED_RGB_S3TC_DXT1_EXT;
        if ( p === RGBA_S3TC_DXT1_Format ) return extension.COMPRESSED_RGBA_S3TC_DXT1_EXT;
        if ( p === RGBA_S3TC_DXT3_Format ) return extension.COMPRESSED_RGBA_S3TC_DXT3_EXT;
        if ( p === RGBA_S3TC_DXT5_Format ) return extension.COMPRESSED_RGBA_S3TC_DXT5_EXT;
      }
    }

    if ( p === RGB_PVRTC_4BPPV1_Format || p === RGB_PVRTC_2BPPV1_Format ||
    p === RGBA_PVRTC_4BPPV1_Format || p === RGBA_PVRTC_2BPPV1_Format ) {
      extension = extensions.get( 'WEBGL_compressed_texture_pvrtc' );

      if ( extension !== null ) {
        if ( p === RGB_PVRTC_4BPPV1_Format ) return extension.COMPRESSED_RGB_PVRTC_4BPPV1_IMG;
        if ( p === RGB_PVRTC_2BPPV1_Format ) return extension.COMPRESSED_RGB_PVRTC_2BPPV1_IMG;
        if ( p === RGBA_PVRTC_4BPPV1_Format ) return extension.COMPRESSED_RGBA_PVRTC_4BPPV1_IMG;
        if ( p === RGBA_PVRTC_2BPPV1_Format ) return extension.COMPRESSED_RGBA_PVRTC_2BPPV1_IMG;
      }
    }

    if ( p === RGB_ETC1_Format ) {
      extension = extensions.get( 'WEBGL_compressed_texture_etc1' );

      if ( extension !== null ) return extension.COMPRESSED_RGB_ETC1_WEBGL;
    }

    if ( p === MinEquation || p === MaxEquation ) {
      extension = extensions.get( 'EXT_blend_minmax' );

      if ( extension !== null ) {
        if ( p === MinEquation ) return extension.MIN_EXT;
        if ( p === MaxEquation ) return extension.MAX_EXT;
      }
    }

    if ( p === UnsignedInt248Type ) {
      extension = extensions.get( 'WEBGL_depth_texture' );

      if ( extension !== null ) return extension.UNSIGNED_INT_24_8_WEBGL;
    }

    return 0;
  }
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

  this.background = null;
  this.fog = null;
  this.overrideMaterial = null;

  this.autoUpdate = true; // checked by the renderer
}

Scene.prototype = Object.create( Object3D.prototype );

Scene.prototype.constructor = Scene;


/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
*
* parameters = {
*  color: <hex>,
*  opacity: <float>,
*
*  linewidth: <float>,
*  linecap: "round",
*  linejoin: "round"
* }
*/

function LineBasicMaterial( parameters ) {
  Material.call( this );

  this.type = 'LineBasicMaterial';

  this.color = new Color( 0xffffff );

  this.linewidth = 1;

  this.lights = false;

  this.setValues( parameters );
}

LineBasicMaterial.prototype = Object.create( Material.prototype );
LineBasicMaterial.prototype.constructor = LineBasicMaterial;

LineBasicMaterial.prototype.isLineBasicMaterial = true;


/**
* @author mrdoob / http://mrdoob.com/
*/

function Line( geometry, material, mode ) {
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

/**
* @author mrdoob / http://mrdoob.com/
* @author alteredq / http://alteredqualia.com/
*/

function AmbientLight( color ) {
  Object3D.call( this );
  this.type = 'AmbientLight';
  this.color = new Color( color );
}

AmbientLight.prototype = Object.assign( Object.create( Object3D.prototype ), {
  constructor: AmbientLight,
  isAmbientLight: true,
} );

/**
* @author mrdoob / http://mrdoob.com/
* @author bhouston / http://clara.io/
* @author stephomi / http://stephaneginier.com/
*/

function Raycaster( origin, direction, near, far ) {
  this.ray = new Ray( origin, direction );
  // direction is assumed to be normalized (for accurate distance calculations)

  this.near = near || 0;
  this.far = far || Infinity;
}

function ascSort( a, b ) {
  return a.distance - b.distance;
}

function intersectObject( object, raycaster, intersects, recursive ) {
  if ( object.visible === false ) return;

  object.raycast( raycaster, intersects );

  if ( recursive === true ) {
    let children = object.children;

    for ( let i = 0, l = children.length; i < l; i ++ ) {
      intersectObject( children[i], raycaster, intersects, true );
    }
  }
}

//

Raycaster.prototype = {

  constructor: Raycaster,

  linePrecision: 1,

  setFromCamera: function ( coords, camera ) {
    if ( (camera && camera.isOrthographicCamera) ) {
      this.ray.origin.set( coords.x, coords.y, ( camera.near + camera.far ) / ( camera.near - camera.far ) ).unproject( camera ); // set origin in plane of camera
      this.ray.direction.set( 0, 0, - 1 ).transformDirection( camera.matrixWorld );
    } else {
      console.error( 'THREE.Raycaster: Unsupported camera type.' );
    }
  },

  intersectObjects: function ( objects, recursive ) {
    let intersects = [];
    for ( let i = 0, l = objects.length; i < l; i ++ ) {
      intersectObject( objects[i], this, intersects, recursive );
    }
    intersects.sort( ascSort );
    return intersects;
  },
};

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
  LineBasicMaterial,
  AmbientLight,
  OrthographicCamera,
  BufferGeometry,
  BufferAttribute,
  Object3D,
  Raycaster,
  Ray,
  Matrix4,
  Vector3,
  Vector2,
  Quaternion,
  Color,
  CatmullRomCurve3,
  TriangleStripDrawMode,
  VertexColors,
  Texture,
};