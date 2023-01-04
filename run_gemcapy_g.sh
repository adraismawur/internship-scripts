#!/bin/bash
mkdir output
for gbk in gbk_out/* ; do
    file=${gbk##*/}
    python $1 --macse_path $2 -t $file -i $3/$file -o output/$file --output_algn -g $4/$file/$file\_GeneFamilies.tsv
done
