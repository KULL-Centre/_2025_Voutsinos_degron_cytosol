#!/bin/bash

# Merge counts files and map to the 5 libraries
Rscript merge_and_map.r samples.csv ../call_zerotol/*_counts.txt > merge_and_map.out

# Calculate degron scores
Rscript scores.r counts.rda > scores.out
