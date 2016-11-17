#!/bin/bash
# note that samtools will ignore duplicates and Q<13
# https://www.biostars.org/p/60841/
# This may cause these results to be slightly different from IGV, depending on what we use

refpath=/seq/hacohenlab1/yjiao/prepost/voting/pileup/Homo_sapiens_assembly19.fasta
pid=$1
dir=$2
reg=$dir/${pid}_pass.region
bamlist=$dir/${pid}.bamlist
out=$dir/${pid}.pileup

while read region; do
    samtools mpileup --max-depth 10000 -f $refpath --region $region `cat $bamlist` >> $out
done < $reg