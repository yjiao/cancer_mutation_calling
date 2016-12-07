#!/bin/bash

field=maf_file_oxoGInfo_capture
out=fh_info.tsv
prefix=/seq/hacohenlab1/yjiao/prepost/voting/forcecalled/postprocessed/appendOxoGInfo

echo -e "pair_id\t$field" > $out # header, flush original file

while read pair; do
    echo -e "$pair\t$prefix/$pair.oxoGInfo.maf.annotated"
done < pairs.txt >> $out

