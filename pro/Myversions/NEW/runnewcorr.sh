#!/bin/sh

idl<<EOF
run=756
camcol=1
rerun=1
newcorr, run, camcol, rerun
EOF
idl<<EOF
run=756
camcol=2
rerun=1
newcorr, run, camcol, rerun
EOF
idl<<EOF
run=756
camcol=3
rerun=1
newcorr, run, camcol, rerun
EOF
idl<<EOF
run=756
camcol=4
rerun=1
newcorr, run, camcol, rerun
EOF
idl<<EOF
run=756
camcol=5
rerun=1
newcorr, run, camcol, rerun
EOF
idl<<EOF
run=756
camcol=6
rerun=1
newcorr, run, camcol, rerun
EOF
