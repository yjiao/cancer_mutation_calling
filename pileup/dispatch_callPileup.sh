#!/bin/bash
dir=/seq/hacohenlab1/yjiao/prepost/voting/pileup

for pid in *.region; do
    pid=${pid%_pass.region}
    qsub -b y -q long -l h_vmem=5g -N mpileup_$pid -cwd "$dir/callPileup.sh $pid $dir"
done