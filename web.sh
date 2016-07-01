#!/bin/sh -eu
# takes one arg - directory of a clone of uglymol.github.io

outdir="$1"

cat >$outdir/index.md <<EOF
---
layout: default
---

$(cat README.md)
EOF

cp src/* $outdir/src/
cp uglymol.js uglymol.min.js LICENSE perf.html $outdir/
cp perf/* $outdir/perf/
