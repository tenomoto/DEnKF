#!/bin/sh
basedir=$(cd $(dirname $0); pwd)
rundir=${basedir}/cycle
genobs=${basedir}/genobs
obsdir=${rundir}/obs
if [ ! -d ${obsdir} ]; then
  mkdir ${obsdir}
fi
cd ${rundir}
cat << EOF | ${genobs}
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
&denkf_true
  tdir = "true"
/
&denkf_obs
  odir = "obs"
  nobs = 300
  rsig = 4.0d0
/
EOF
