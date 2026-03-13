const js = require('@eslint/js');
const tseslint = require('typescript-eslint');

const sharedGlobals = {
  window: 'readonly',
  document: 'readonly',
  navigator: 'readonly',
  location: 'readonly',
  console: 'readonly',
  fetch: 'readonly',
  Blob: 'readonly',
  FileReader: 'readonly',
  XMLHttpRequest: 'readonly',
  performance: 'readonly',
  requestAnimationFrame: 'readonly',
  cancelAnimationFrame: 'readonly',
  setTimeout: 'readonly',
  clearTimeout: 'readonly',
  URL: 'readonly',
  Image: 'readonly',
  Audio: 'readonly',
  DOMParser: 'readonly',
  process: 'readonly',
  require: 'readonly',
  module: 'readonly',
  exports: 'readonly',
  __dirname: 'readonly',
  global: 'readonly',
  Float32Array: 'readonly',
  Float64Array: 'readonly',
  Int8Array: 'readonly',
  Int16Array: 'readonly',
  Int32Array: 'readonly',
  Uint8Array: 'readonly',
  Uint16Array: 'readonly',
  Uint32Array: 'readonly',
  ArrayBuffer: 'readonly',
  DataView: 'readonly',
  VERSION: 'readonly',
};

const styleRules = {
  indent: ['error', 2, {
    SwitchCase: 1,
    outerIIFEBody: 0,
    ArrayExpression: 'first',
    ObjectExpression: 'first',
    ImportDeclaration: 'first',
    CallExpression: {arguments: 'first'},
    FunctionDeclaration: {parameters: 'first'},
    FunctionExpression: {parameters: 'first'},
    ignoredNodes: ['ConditionalExpression'],
  }],
  'brace-style': ['error', '1tbs', {allowSingleLine: true}],
  'block-spacing': ['error', 'always'],
  'space-before-function-paren': ['error', {anonymous: 'always', named: 'never'}],
  curly: ['error', 'multi-line', 'consistent'],
  'spaced-comment': 'off',
  'object-curly-spacing': 'off',
  'operator-linebreak': ['error', 'after', {overrides: {':': 'ignore'}}],
  camelcase: 'off',
  'linebreak-style': 'off',
  'require-jsdoc': 'off',
  'no-console': 'off',
  'prefer-const': 'off',
  'prefer-spread': 'off',
  'no-multi-spaces': ['error', {ignoreEOLComments: true}],
};

const tsFiles = ['**/*.ts', '**/*.tsx', '**/*.mts', '**/*.cts'];

module.exports = [
  {
    ignores: ['coverage/**', 'node_modules/**', 'src/wasm/*', 'src/three/*', 'three.module.js'],
  },
  js.configs.recommended,
  {
    files: ['**/*.js', '**/*.cjs', '**/*.mjs'],
    languageOptions: {
      ecmaVersion: 'latest',
      sourceType: 'module',
      globals: {
        ...sharedGlobals,
        jest: 'readonly',
        test: 'readonly',
        expect: 'readonly',
        beforeAll: 'readonly',
        afterAll: 'readonly',
        describe: 'readonly',
        it: 'readonly',
        Benchmark: 'readonly',
      },
    },
    rules: styleRules,
  },
  ...tseslint.configs.recommended.map((config) => (
    config.files ? config : {...config, files: tsFiles}
  )),
  {
    files: tsFiles,
    languageOptions: {
      ecmaVersion: 'latest',
      sourceType: 'module',
      globals: sharedGlobals,
    },
    rules: {
      ...styleRules,
      '@typescript-eslint/no-explicit-any': 'off',
      '@typescript-eslint/no-this-alias': 'off',
    },
  },
];
