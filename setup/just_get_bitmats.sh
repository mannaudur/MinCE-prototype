#!/bin/bash

i=0
filename='atomsleft.txt'
while read line; do
    Bifrost build -r $line -c -o ${line%%.*}
    Bifrost build -r $line -a -o ${line%%.*}
    Bifrost query -g ${line%%.*}.gfa -f ${line%%.*}.bfg_colors -q ${line%%.*}.fasta -o ${line%%.*}
    ((i=i+1))
    echo "$i bitmats done"
done < $filename
