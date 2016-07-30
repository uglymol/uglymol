
- **Important**:
  For now, like in Coot, bonds are drawn using wide GL_LINES,
  but some browser, in particular on Windows, show only thin lines (width 1).
  * use own shaders for trace (done)
  * use own shaders for bonds

  Temporary workaround: [switch browser settings](https://github.com/mrdoob/three.js/wiki/How-to-use-OpenGL-or-ANGLE-rendering-on-Windows)
  from the ANGLE backend to OpenGL.

- It would nice to calculate maps from a reflection file (mtz)
  in the browser, but currently no one is working on it. Volunteers needed.
  Requires 3D FFT.

- atom labels

- use UMD format for built files. Maybe use rollupjs
  if https://github.com/mrdoob/three.js/pull/9310 gets merged.

- "ribbon" representation like in xtal.js

- move isosurface calculation to update()

- test it on a phone

- a few more keybindings:
  L = next ligand (ctrl-L in coot),
  u = move back to old centre,
  Ctrl-G - go to residue

