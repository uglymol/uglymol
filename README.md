UglyMol is a web-based macromolecular viewer focused on electron density.

It makes models and e.den. maps easy to recognize, navigate and interpret --
for [Coot](http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/) users.
It looks like coot and walks (mouse controls) like coot.
Except when it doesn't -- see the
[FAQ](https://github.com/uglymol/uglymol/wiki).

It's only a viewer. For situations when you want
a quick look without downloading the data and starting Coot.
Say, for screening results from
the [dimple](http://ccp4.github.io/dimple/) pipeline in a synchrotron.
Of course, for this to work, it needs to be integrated into a website
that provides the data access.

Have a look:
[1MRU](https://uglymol.github.io/1mru.html) (60kDa, 3A),
xxxx (>200kDa),
dimple blob,
SAD phasing.

(if it doesn't look well and you are on Windows -- see the first
point [here](https://github.com/uglymol/uglymol/blob/master/TODO.md))

Technically, UglyMol is a small
[project](https://github.com/uglymol/uglymol) (~2 KLOC)
forked from Nat Echols' [xtal.js](https://github.com/natechols/xtal.js/).
The only dependency is [three.js](https://github.com/mrdoob/three.js/).
Built files (uglymol.js and .min.js) are
[here](https://github.com/uglymol/uglymol.github.io)
and in [npm](https://www.npmjs.com/package/uglymol)
and [bower]().

The [plan](https://github.com/uglymol/uglymol/blob/master/TODO.md)
is to keep UglyMol small and ugly rather than to add many features.
And to make it as [fast](https://uglymol.github.io/perf.html) as possible.
Actually this project is an experiment and further development, if any,
will depend on received feedback. So, what do you think should be added
or changed?
Use the [issue tracker](https://github.com/uglymol/uglymol/issues)
or gitter [chat](https://gitter.im/ccp4/dimple)
or email wojdyr@gmail.com.

Finally, if you needs pretty visualization,
or if you trust the model and don't need to see the electron density,
then UglyMol is not a good fit. Try one of many
[other viewers](https://github.com/uglymol/uglymol/wiki/MolecularViewers).
