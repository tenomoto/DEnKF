#!/bin/sh
basedir=$(cd $(dirname $0); pwd)
rundir=${basedir}/cycle
genloc=${basedir}/genloc

cd ${rundir}
cat << EOF | ${genloc}
&qg_grid
  imax = 129
  jmax = 129
/
&loc_param
  r0 = 15.d0
/
EOF
