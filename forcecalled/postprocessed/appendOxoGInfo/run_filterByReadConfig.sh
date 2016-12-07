#!/bin/bash

while read pair; do
    python filterByReadConfig.py --onlyAddColumnsToCopy ${pair}.oxoG.metrics.txt ${pair}.picard_oxoQ.maf.annotated ${pair}.oxoGInfo.maf.annotated
done < pairs.txt
