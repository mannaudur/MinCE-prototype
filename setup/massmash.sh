#!/bin/bash

for file in $(find toyset -type f); do smash/bin/sketch $file -d toyset/toysketches/; done
