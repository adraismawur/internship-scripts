#!/bin/bash
mkdir $4
for gbk in $3/* ; do
    file=${gbk##*/}
    python $1 --macse_path $2 -t $file -i $3/$file -o $4/$file --output_algn -g $5/$file/$file\_GeneFamilies.tsv $6
done
