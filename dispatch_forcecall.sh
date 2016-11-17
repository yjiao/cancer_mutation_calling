#!/bin/bash

prefix=/seq/hacohenlab1/yjiao/prepost/voting/pileup/
mutations=/seq/hacohenlab1/yjiao/prepost/voting/pass
prefix=/seq/hacohenlab1/yjiao/prepost/voting/pileup/
outpath=/seq/hacohenlab1/yjiao/prepost/voting/forcecalled
dir=/seq/hacohenlab1/yjiao/prepost/voting


for pid in $prefix/*.pileup; do
    pid=${pid%.pileup}
    pid=${pid##$prefix/}
    echo $pid
    
    qsub -b y -q long -l h_vmem=5g -N forcecall_$pid -cwd "$dir/forcecall $pid $prefix $mutations $outpath"
done