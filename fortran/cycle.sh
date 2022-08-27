#!/bin/sh
basedir=$(cd $(dirname $0); pwd)
exp=cycle
expdir=${basedir}/${exp}
denkf=${basedir}/run_denkf
qg=${basedir}/run_qg
ncycle=100
nens=25
nstep=4
n=16641

cd ${expdir}
for i in $(seq 0 ${ncycle}); do
  step=$((i * nstep))
  ssssss=$(printf %06d ${step})
  ((i1 = i + 1))
  step1=$((${i1} * nstep))
  ssssss1=$(printf %06d ${step1})
  obsoff=$(awk "NR==${i1}{print}" ${expdir}/obs/obsoff.txt)
  echo "DEnKF cycle ${i}"
  cat << EOF | ${denkf}
&qg_grid
  imax = 129
  jmax = 129
/
&denkf_cycle
  inflation_factor = 1.06
  localize = .false.
  step = ${step}
  nstep = ${nstep}
  nens = ${nens}
  obsoff = ${obsoff}
  fdir = "forecast"
  adir = "analysis"
/
&denkf_obs
  odir = "obs"
  nobs = 300
  rsig = 4.0d0
/
EOF
  if [ ${i} -lt ${ncycle} ]; then
    for e in $(seq ${nens}); do
      echo "Member ${e}"
      dd if=${expdir}/analysis/q${ssssss}.dat of=q000000.dat bs=4 skip=$((n * (e - 1))) count=${n}
      dd if=${expdir}/analysis/p${ssssss}.dat of=p000000.dat bs=4 skip=$((n * (e - 1))) count=${n}
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
  nstep = ${nstep}
  tsave = 0
  nsave = ${nstep}
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
      wait
      mv q$(printf %06d ${nstep}).dat q$(printf %03d ${e}).dat
      mv p$(printf %06d ${nstep}).dat p$(printf %03d ${e}).dat
  done
  cat q???.dat > ${expdir}/forecast/q${ssssss1}.dat
  cat p???.dat > ${expdir}/forecast/p${ssssss1}.dat
  rm q???.dat p???.dat q000000.dat p000000.dat
  fi
done
ls -l forecast
ls -l analysis
ls -l
