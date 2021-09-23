#!/bin/bash

for FILE in $(find atoms/ -type f -name "*.txt")
do 
    Bifrost build -r $FILE -c -o ${FILE%%.*}
    Bifrost build -r $FILE -a -o ${FILE%%.*}
    Bifrost query -g ${FILE%%.*}.gfa -f ${FILE%%.*}.bfg_colors -q ${FILE%%.*}.fasta -o ${FILE%%.*}
    python3 extract_features.py ${FILE%%.*}.tsv ${FILE%%.*}.fasta -d features/
    echo "${FILE%%.*} done"
done