#!/bin/bash
mkdir output
for gbk in gbk_out/* ; do
    file=${gbk##*/}
    python $1 --macse_path $2 -t $file -i gbk_out/$file -o output/$file --output_algn
done
