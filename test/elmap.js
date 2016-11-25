
var assert = require('chai').assert;
var ElMap = require('../uglymol').ElMap;
var util = require('../perf/util');


/* Note: axis order in ccp4 maps is tricky. It was tested by re-sectioning,
 * i.e. changing axis order, of a map with CCP4 mapmask:
  mapmask mapin 1mru_2mFo-DFc.ccp4 mapout 1mru_yzx.ccp4 << eof
  AXIS Y Z X
  MODE mapin
  eof
*/

describe('ElMap', function () {
  'use strict';
  var dmap_buf = util.open_as_array_buffer('1mru.omap');
  var cmap_buf = util.open_as_array_buffer('1mru.map');
  var dmap = new ElMap();
  var cmap = new ElMap();
  it('#from_dsn6', function () {
    dmap.from_dsn6(dmap_buf);
  });
  it('#from_ccp4', function () {
    cmap.from_ccp4(cmap_buf);
  });
  it('compare unit cells', function () {
    for (var i = 0; i < 6; i++) {
      assert.closeTo(dmap.unit_cell.parameters[i],
                     cmap.unit_cell.parameters[i], 0.02);
    }
  });
});

