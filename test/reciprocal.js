
var UM = require('../uglymol');
//var util = require('../perf/util');

describe('ReciprocalViewer', function () {
  'use strict';
  var viewer = new UM.ReciprocalViewer('rv');
  //var data = util.open_as_utf8('rlp.csv');
  //viewer.load_from_string(data, false);
  it('misc calls', function () {
    viewer.controls.update();
    viewer.update_camera();
    viewer.shift_clip();
    viewer.recenter();
    viewer.recenter([1, 2, 3]);
  });
});
