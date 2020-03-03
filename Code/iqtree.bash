#!/usr/bin/env bash

# IQ-TREE was downloaded using the conda distribution package
# see here https://anaconda.org/bioconda/iqtree

# -s - alignment file
# -bb - rapid bootstrap with 1000 replicates
# -pre - prefix of output file
# -g - add a constraint tree to constrain the tree to a known topology
# -nt - auto detect number of threads
# -redo - redo and overwrite all analyses (do not consider checkpoints)
# -m - select best fit substitution model and merge partitions to stop model overfitting.
# -spp - specify partition file allowing for different evolutionary models per substitution

iqtree -s Nucleotide\ alignment\ 7.phy -bb 1000 -pre 270220 -g constraint_tree_exp_1_3_unrooted.newick -nt AUTO -redo -m TESTNEWMERGE -spp Partition_file.nexus
