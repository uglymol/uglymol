// Check that scripts in perf/ can be run without errors.
// To make it quick use the least possible number of iterations.
var Benchmark = require('benchmark');
Benchmark.options.maxTime = -Infinity;
Benchmark.options.minTime = -Infinity;
Benchmark.options.minSamples = 1;
Benchmark.options.initCount = 0;

describe('perf', function () {
  'use strict';
  function add(name) {
    it(name, function () {
      var save_console_log = console.log;
      console.log = function () {};
      try {
        require('../perf/' + name);
      } finally {
        console.log = save_console_log;
      }
    });
  }
  add('model');
  add('elmap');
  add('isosurface');
  add('viewer');
});

