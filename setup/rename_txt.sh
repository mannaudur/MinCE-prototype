#!/bin/bash

for f in $(find atoms/ -type f); do 
    mv "$f" "$f.txt"
done