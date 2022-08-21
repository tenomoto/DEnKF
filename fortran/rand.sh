#!/bin/sh
basedir=$(cd $(dirname $0); pwd)
rundir=${basedir}
nens=25
outfile=${basedir}/rand.txt
if [ -f ${outfile} ]; then
  rm -f ${outfile}
fi
rm -f 
for i in $(seq 0 ${nens}); do
  echo $(( ($(od -An -N4 -tu4 /dev/random) % 4001) * 100 + 100000 )) >> ${outfile}
done
