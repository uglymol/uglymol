'use strict';

var Benchmark = require('benchmark');
var suite = new Benchmark.Suite();

var VALUES = [21.2, 30.5, 25.9, -7.2, 0.4, -0.4];
//var VALUES = [21, 30, 26, -7, 0, -1];
var M = 3;
var result1;
var result2;

// like Math.euclideanModulo() from three.js
suite.add('%,+,%', function () {
  result1 = 0;
  for (var i = 0; i < VALUES.length; i++) {
    var n = VALUES[i];
    result1 += ((n % M) + M) % M;
  }
});

suite.add('%,if,+=', function () {
  result2 = 0;
  for (var i = 0; i < VALUES.length; i++) {
    var n = VALUES[i];
    var x = n % M;
    if (x < 0) {
      x += M;
    }
    result2 += x;
  }
});

suite.on('cycle', function (event) { console.log(' ' + event.target); });
suite.run();
console.log('  results: ' + result1 + ' and ' + result2);

