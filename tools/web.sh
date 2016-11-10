#!/bin/sh -eu
# takes one arg - directory of a clone of uglymol.github.io
outdir="$1"

[ -e src/elmap.js ] || { echo "Run me from top-level uglymol dir"; exit 1; }

# use README.md without badges
cat >$outdir/index.md <<EOF
---
layout: default
---

$(cat README.md | grep -v '^\[!\[' | sed s,https://uglymol.github.io/,,)
EOF

cp src/* $outdir/src/
for path in three/build/three.min.js benchmark/benchmark.js \
            lodash/lodash.min.js platform/platform.js; do
    npath=node_modules/$path
    diff -q $npath $outdir/$npath || cp $npath $outdir/$npath
done

cp uglymol.js uglymol.js.map uglymol.min.js LICENSE perf.html $outdir/
cp perf/* $outdir/perf/
cp test/*.html $outdir/test/

grep -v -- -DEV- dev.html > $outdir/1mru.html
grep -v -- -DEV- dual.html > $outdir/dual.html
grep -v -- -DEV- reciprocal.html > $outdir/reciprocal.html
sed -e s/1mru/4un4_final/g -e s/\\.map/_m0.map/g <$outdir/1mru.html \
    >$outdir/4un4.html
sed -e s/1mru/dimple_thaum/g \
    < $outdir/1mru.html > $outdir/dimple_thaum.html

cd $outdir
echo "=== $(pwd) ==="
git status -s
echo
git diff --stat
echo
jekyll build
du -sh _site/
