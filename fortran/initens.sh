#!/bin/sh
basedir=$(cd $(dirname $0); pwd)
freedir=${basedir}/free
rundir=${basedir}/cycle
rfile=${basedir}/rand.txt
nens=25
for d in forecast analysis; do
  if [ -d ${rundir}/${d} ]; then
    rm -rf ${rundir}/${d}
  fi
  mkdir ${rundir}/${d}
done
t=$(head -1 ${rfile})
t=$(awk "NR==2{print}" ${rfile})
#((t = t - nens / 2 * 100))
#if [ ${t} -lt 100000 ]; then
#  t=100000
#fi
for x in q p; do
  for i in $(seq ${nens}); do
    ii=$(printf %02d ${i})
    t=$(awk "NR==${i}+1{print}" ${rfile})
#    ((t = t + 100))
    cp ${freedir}/${x}${t}.dat ${rundir}/${x}${ii}.dat
  done
  cat ${rundir}/${x}??.dat > ${rundir}/forecast/${x}000000.dat 
  rm ${rundir}/${x}??.dat
done
