var UnitCell = require('../uglymol').UnitCell;

function assert_equal_arrays(a, b, rtol) {
  'use strict';
  expect(a.length).toEqual(b.length);
  for (var i = 0; i < a.length; i++) {
    var eps = rtol * Math.max(Math.abs(a[i]), Math.abs(b[i]));
    expect(Math.abs(a[i] - b[i])).toBeLessThanOrEqual(eps);
  }
}

describe('UnitCell', () => {
  'use strict';
  it('#orthogonalize', () => {
    // derived from dimple/tests/utest.py
    var uc = new UnitCell(22.84, 32.84, 42.84, 80.84, 90.84, 100.84);
    expect(uc).toBeTruthy();
    var orth = uc.orthogonalize.bind(uc);
    assert_equal_arrays([22.84, 0, 0], orth([1, 0, 0]), 1e-6);
    assert_equal_arrays([-6.17612, 32.254, 0], orth([0, 1, 0]), 1e-6);
    assert_equal_arrays([-0.628045, 6.82343, 42.2884], orth([0, 0, 1]), 1e-6);
  });
  it('#orthogonalize * #fractionalize', () => {
    var uc = [new UnitCell(63.10, 50.17, 111.07, 90, 96.19, 90.00),
              new UnitCell(70, 80, 90, 48, 58, 68),
              new UnitCell(27.07, 31.25, 33.76, 87.98, 108.0, 112.11)];
    for (var i = 0; i < uc.length; ++i) {
      var frac = [Math.random(), Math.random(), Math.random()];
      var orth = uc[i].orthogonalize(frac);
      var frac2 = uc[i].fractionalize(orth);
      assert_equal_arrays(frac, frac2, 1e-12);
    }
  });
});

