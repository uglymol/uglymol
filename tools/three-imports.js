/*
 * This is import "proxy" for the subset of Three.js used in UglyMol.
 * It is used when making single bundle (no dependencies).
 * We avoid importing Three.js to make Rollup tree-shaking more effective.
 */
import '../node_modules/three/src/polyfills.js';
export { WebGLRenderer } from '../node_modules/three/src/renderers/WebGLRenderer.js';
export { Fog } from '../node_modules/three/src/scenes/Fog.js';
export { Scene } from '../node_modules/three/src/scenes/Scene.js';
export { Mesh } from '../node_modules/three/src/objects/Mesh.js';
export { LineSegments } from '../node_modules/three/src/objects/LineSegments.js';
export { Line } from '../node_modules/three/src/objects/Line.js';
export { Points } from '../node_modules/three/src/objects/Points.js';
export { ShaderMaterial } from '../node_modules/three/src/materials/ShaderMaterial.js';
export { PointsMaterial } from '../node_modules/three/src/materials/PointsMaterial.js';
export { MeshBasicMaterial } from '../node_modules/three/src/materials/MeshBasicMaterial.js';
export { LineBasicMaterial } from '../node_modules/three/src/materials/LineBasicMaterial.js';
export { AmbientLight } from '../node_modules/three/src/lights/AmbientLight.js';
export { OrthographicCamera } from '../node_modules/three/src/cameras/OrthographicCamera.js';
export { BufferGeometry } from '../node_modules/three/src/core/BufferGeometry.js';
export { Geometry } from '../node_modules/three/src/core/Geometry.js';
export { BufferAttribute } from '../node_modules/three/src/core/BufferAttribute.js';
export { Object3D } from '../node_modules/three/src/core/Object3D.js';
export { Raycaster } from '../node_modules/three/src/core/Raycaster.js';
export { Ray } from '../node_modules/three/src/math/Ray.js';
export { Matrix4 } from '../node_modules/three/src/math/Matrix4.js';
export { Vector3 } from '../node_modules/three/src/math/Vector3.js';
export { Vector2 } from '../node_modules/three/src/math/Vector2.js';
export { Quaternion } from '../node_modules/three/src/math/Quaternion.js';
export { Color } from '../node_modules/three/src/math/Color.js';
export { CatmullRomCurve3 } from '../node_modules/three/src/extras/curves/CatmullRomCurve3.js';
export { TriangleStripDrawMode, VertexColors } from '../node_modules/three/src/constants.js';
