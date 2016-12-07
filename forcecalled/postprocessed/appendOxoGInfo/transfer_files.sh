#!/bin/bash

prefix=/local/cga-fh/cga/Tumor_Immunity/Pair/
suffix1=/jobs/capture/mut/oxog/postreal/append/latest/
suffix2=/jobs/capture/mut/oxog/metrics/latest/
newpath=/seq/hacohenlab1/yjiao/prepost/voting/forcecalled/postprocessed/appendOxoGInfo/

while read pair; do
    cp -f $prefix$pair$suffix1$pair* $newpath
    cp -f $prefix$pair$suffix2$pair* $newpath 
done < pairs.txt