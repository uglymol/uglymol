
var assert = require('chai').assert;
var fs = require('fs');
var ElMap = require('../src/elmap');

// http://stackoverflow.com/a/12101012/104453
function toArrayBuffer(buffer) {
  var ab = new ArrayBuffer(buffer.length);
  var view = new Uint8Array(ab);
  for (var i = 0; i < buffer.length; ++i) {
    view[i] = buffer[i];
  }
  return ab;
}

describe('ElMap', function () {
  'use strict';
  var dsn6_file = __dirname + '/../data/1mru.omap';
  var ccp4_file = __dirname + '/../data/1mru_2mFo-DFc.ccp4';
  var dmap_buf = toArrayBuffer(fs.readFileSync(dsn6_file));
  var cmap_buf = toArrayBuffer(fs.readFileSync(ccp4_file));
  var dmap = new ElMap();
  var cmap = new ElMap();
  it('#from_dsn6', function () {
    dmap.from_dsn6(dmap_buf);
  });
  it('#from_ccp4', function () {
    cmap.from_ccp4(cmap_buf);
  });
  it('should have same the unit cell', function () {
    for (var i = 0; i < 6; i++) {
      assert.closeTo(dmap.unit_cell.parameters[i],
                     cmap.unit_cell.parameters[i], 0.02);
    }
  });
});

