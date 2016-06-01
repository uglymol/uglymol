'use strict';

var fs = require('fs');
var Benchmark = require('benchmark');

function data_path(filename) {
  return __dirname + '/../data/' + filename;
}

function open_as_utf8(filename) {
  var path = data_path(filename);
  return fs.readFileSync(path, {encoding: 'utf8'});
}

function open_as_array_buffer(filename) {
  var path = data_path(filename);
  var buffer = fs.readFileSync(path);
  // http://stackoverflow.com/a/12101012/104453
  var ab = new ArrayBuffer(buffer.length);
  var view = new Uint8Array(ab);
  for (var i = 0; i < buffer.length; ++i) {
    view[i] = buffer[i];
  }
  return ab;
}

function new_benchmark_suite() {
  var suite = new Benchmark.Suite();
  suite.on('cycle', function (event) {
    console.log(' ' + event.target);
  });
  return suite;
}

exports.data_path = data_path;
exports.open_as_array_buffer = open_as_array_buffer;
exports.new_benchmark_suite = new_benchmark_suite;
exports.open_as_utf8 = open_as_utf8;
