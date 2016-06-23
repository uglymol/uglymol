'use strict';
// Node.js only

var fs = require('fs');
var Benchmark = require('benchmark');

function data_path(filename) {
  return __dirname + '/../data/' + filename;
}

exports.open_as_utf8 = function (filename) {
  var path = data_path(filename);
  return fs.readFileSync(path, {encoding: 'utf8'});
};

exports.open_as_array_buffer = function (filename) {
  var path = data_path(filename);
  var buffer = fs.readFileSync(path);
  // http://stackoverflow.com/a/12101012/104453
  var ab = new ArrayBuffer(buffer.length);
  var view = new Uint8Array(ab);
  for (var i = 0; i < buffer.length; ++i) {
    view[i] = buffer[i];
  }
  return ab;
};

var bench_to_run = process.argv[2];

exports.bench = function (name, fn, options) {
  var b = new Benchmark(name, fn, options);
  //b.on('start', function () { console.log('started ' + b.name); });
  b.on('complete', function (event) {
    console.log(' ' + event.target);
  });
  b.on('error', function () {
    console.log(b.error);
  });
  if (bench_to_run === undefined || name.indexOf(bench_to_run) > -1) {
    b.run();
  } else {
    b.fn(); // run once, for possible side effects
  }
  return b;
};

