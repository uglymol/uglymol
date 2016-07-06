
- **Important**:
  For now, like in Coot, bonds are drawn using wide GL_LINES,
  but some browser, in particular on Windows, show only thin lines (width 1).
  It needs to be fixed in UglyMol - i.e. we cannot rely on wide GL_LINES
  and should write shaders to more reliably render bonds.

  Temporary workaround: [switch browser settings](https://github.com/mrdoob/three.js/wiki/How-to-use-OpenGL-or-ANGLE-rendering-on-Windows)
  from the ANGLE backend to OpenGL.

- **Important**: no refreshing when nothing moves

- It would nice to calculate maps from a reflection file (mtz)
  in the browser, but currently no one is working on it. Volunteers needed.
  Requires 3D FFT.

- catch errors from reading data files

- add a test case with hydrogens and multiple conformers,
  and a keybinding to hide/show hydrogens

- move isosurface calculation to update()

- test it on a phone

- permalinks to xyz positions

- Compare performance of marching cubes with
  github.com/mikolalysenko/isosurface and @alteredq's / stemkoski's demos.
  xtal.js and uglymol implementation is based on the latter.

- atom labels

- a few more keybindings:
  L = next ligand (ctrl-L in coot),
  u = move back to old centre,
  Ctrl-G - go to residue

