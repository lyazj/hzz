#!/bin/bash

exec > merge.txt
for YEAR in 2016APV 2016 2017 2018; do
    for INDIR in $(ls -vd /eos/user/l/legao/NanoAODStore/V0/${YEAR}/MC/*/*); do
        OUTDIR=${INDIR/NanoAODStore/NanoAOD}
        mkdir -p ${OUTDIR}
        for i in {0..99}; do
            echo "${OUTDIR}/${i}.root ${INDIR}/{$[i*100]..$[i*100+99]}.root"
        done
    done
done
