#!/bin/sh
basedir=$(cd $(dirname $0); pwd)
truedir=${basedir}/cycle/true
qg=${basedir}/run_qg
if [ ! -d ${truedir} ]; then
  mkdir ${truedir}
fi
cd ${truedir}
cat << EOF | ${qg}
&mg
  mg_itermax = 1, 1, 100
  mg_tol = 1.0d-4
/
&qg_grid
  imax = 129
  jmax = 129
/
&qg_time
  init = "file" 
  nstep = 1200 
  tsave = 0
  nsave = 4
  dt = 1.25d0
/
&qg_param
  qg_beta = 1.0d0
  qg_f = 1.6d3
  qg_eps = 1.0d-5
  qg_a = 2.0d-11
  qg_tau0 = -6.283185307179586d0
/
EOF
