#!/bin/sh -eu
# takes one arg - directory of a clone of uglymol.github.io
outdir="$1"

[ -e src/elmap.js ] || { echo "Run me from top-level uglymol dir"; exit 1; }

# use README.md without badges
cat >$outdir/index.md <<EOF
---
layout: default
---

$(cat README.md | grep -v '^\[!\[' | sed s,https://uglymol.github.io,,)
EOF

cp src/* $outdir/src/
cp uglymol.js uglymol.min.js LICENSE perf.html $outdir/
cp perf/* $outdir/perf/

grep -v -- -DEV- dev.html > $outdir/1mru.html
sed -e s/1mru/4un4_final/g -e s/\\.map/_m0.map/g <$outdir/1mru.html \
    >$outdir/4un4.html
sed -e s/1mru/dimple_thaum/g \
    -e 's/pdb")/pdb", {center: [13.98, 18.1, 12.26]})/' \
    < $outdir/1mru.html > $outdir/dimple_thaum.html
# V.load_pdb("data/dimple_thaum.pdb", [13.98, 18.1, 12.26]);

cd $outdir
pwd
git status -s
