
UglyMol is yet another web-based macromolecular viewer, with emphasis
on **electron density maps**.
It makes models and maps easy to recognize, navigate and interpret --
for [coot](http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/) users.
It looks like coot and walks like coot (mouse and some keyboard controls
are the same).

It is only a viewer and it is for situations when you just want
a quick look without downloading the data and starting coot.
For example, when screening results from
the [dimple](http://ccp4.github.io/dimple/) pipeline in a synchrotron.
Of course, for this to work, it needs to be integrated into a website
that provides the data access.

*Caveat*: for now, like in coot, bonds are drawn using wide `GL_LINES`,
which sadly doesn't work in some browser, in particular on Windows.
It will be fixed in UglyMol at some point, but for now switching
browser settings
[from the ANGLE backend to OpenGL](https://github.com/mrdoob/three.js/wiki/How-to-use-OpenGL-or-ANGLE-rendering-on-Windows)
may help.

### Demo: [1mru](https://uglymol.github.io/1mru.html)

Technically, it is a small (~2 KLOC) project forked from Nat Echols'
[xtal.js](https://github.com/natechols/xtal.js/).
The only dependency is [three.js](https://github.com/mrdoob/three.js/).
The plan is to keep it small and ugly rather than to add features.

Actually the project is an experiment and further development, if any,
will depend on received feedback:
email wojdyr@gmail.com or
use [Issues](https://github.com/uglymol/uglymol/issues).

If you needs pretty visualization,
or if you trust the model and don't need to see the electron density,
then UglyMol is not a good fit. Try one of those:
[3Dmol](https://github.com/3dmol/3Dmol.js),
[NGL](https://github.com/arose/ngl),
[jsmol](http://wiki.jmol.org/index.php/JSmol),
[pv](https://github.com/biasmv/pv),
[GLmol](https://github.com/biochem-fan/GLmol).
And if you don't mind Java applets -
[AstexViewer](http://openastexviewer.net/)
and [JMol](http://www.jmol.org/).
