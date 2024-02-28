
var UM = require('../uglymol');
//var util = require('../perf/util');

describe('ReciprocalViewer', () => {
  'use strict';
  var viewer = new UM.ReciprocalViewer();
  //var data = util.open_as_utf8('rlp.csv');
  //viewer.load_from_string(data, false);
  it('misc calls', () => {
    viewer.controls.update();
    viewer.update_camera();
    viewer.shift_clip();
    viewer.recenter();
    viewer.recenter([1, 2, 3]);
  });
});
