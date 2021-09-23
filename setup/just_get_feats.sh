#!/bin/bash

i=0
filename='featsleft.txt'
while read line; do
    python3 deBruijn/extract_features.py $line ${line%%.*}.fasta -d features/
    ((i=i+1))
    echo "$i atoms done"
done < $filename