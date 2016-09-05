'use strict';
// Node.js only

var child_process = require('child_process');
var fs = require('fs');
var Benchmark = require('benchmark');

function data_path(filename) {
  if (filename.charAt(0) === '/') return filename;
  var path = __dirname + '/../data/' + filename;
  try {
    fs.statSync(path);
  } catch (e) {
    console.log('we need to download ' + filename + ' (only once)');
    var url = 'http://uglymol.github.io/data/' + filename;
    var cmd = 'curl -f ' + url + ' -o ' + path;
    if (!child_process.execSync) {
      console.log('no execSync (in Node v0.11.12+), you need to:\n' + cmd);
      process.exit(1);
    }
    child_process.execSync(cmd);
  }
  return path;
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

if (bench_to_run === 'download-data') {
  (function () {
    var files = ['1mru.pdb', '1mru.map', '1mru_diff.map', '1mru_m0.map',
                 '1mru.omap', '1mru_diff.omap'];
    for (var i = 0; i < files.length; i++) {
      data_path(files[i]);
    }
  })();
}
