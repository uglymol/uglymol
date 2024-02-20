
var ElMap = require('../uglymol').ElMap;
var util = require('../perf/util');


/* Note: axis order in ccp4 maps is tricky. It was tested by re-sectioning,
 * i.e. changing axis order, of a map with CCP4 mapmask:
  mapmask mapin 1mru_2mFo-DFc.ccp4 mapout 1mru_yzx.ccp4 << eof
  AXIS Y Z X
  MODE mapin
  eof
*/

describe('ElMap', () => {
  'use strict';
  var dmap_buf = util.open_as_array_buffer('1mru.omap');
  var cmap_buf = util.open_as_array_buffer('1mru.map');
  var dmap = new ElMap();
  var cmap = new ElMap();
  it('#from_dsn6', () => {
    dmap.from_dsn6(dmap_buf);
  });
  it('#from_ccp4', () => {
    cmap.from_ccp4(cmap_buf);
  });
  it('compare unit cells', () => {
    for (var i = 0; i < 6; i++) {
      var p1 = dmap.unit_cell.parameters[i];
      var p2 = cmap.unit_cell.parameters[i];
      expect(Math.abs(p1 - p2)).toBeLessThan(0.02);
    }
  });
});

