# Subtree_Sampling.Rmd

Sampling subtype reference sequences from subtrees of a phylogenetic tree


## Description

For each fasta sequence file in a given directory, perform the following: 
(1) calculate TN93 pairwise distances (PWDs);  (2) create a neighbor joining tree with the PWD matrix and root it at the midpoint; (3) traverse the tree from the root to the leaves and partition it into subtrees having a local median PWD (MPWD) less than a cutoff; (4) select one sequence per country from each subtree. 


This R script is described in Rhee, S-Y and Shafer, RW (2017), "Geographically-Stratified Representative HIV-1 Group M pol Subtype and Circulating Recombinant Form Sequences," manuscript submitted for publication. In the manuscript, the maximum MPWD for the subtrees was defined as the 47-th percentile of the ranked PWDs.




## Usage

To run the script using the default parameters

Rscript -e "rmarkdown::render('Subtree_Sampling.Rmd')"  

Rscript -e "rmarkdown::render('Subtree_Sampling.Rmd')" fastadir="in.fasta", outdir="out.trees", pwd_rank_cutff=0.47



## Input and the default values

1. fastadir ("in.fasta" by default): directory where fasta sequence files are located. The name of each fasta file consists of subtype + ".txt". ie, 01_AE.txt, B.txt. The name of each sequence in a fasta file consists of sequenceid (in this case, accessionid), country (in two letter ISO code), collection year (YY), subtype delimited by "." Notes: example fasta files are available at doc/in.fasta

2. outdir ("out.trees" by default): directory for saving output files including .tre, .pdf, and .txt of sampled sequences

3. pwd_rank_cutff(0.47 by default): the percentile rank of the pairwise distance used as the cutoff value. The clusters are found while traversing from the root to the leaves, such that the the local MPWD (MPWD of the node's children) is less than this cutoff. The sequence space becomes partitioned by the clusters.



## Output

1. Subtree_Sampling.html

2. results.summary.txt: lists number of clusters, number of sequences sampled so on.

Files generated in outdir for each subtype:

3. a tree of the entire sequences (.tre) (.pdf)

4. a tree of the entire sequences with cluster ids defined by a cutoff (.tre)

5. a fasta sequence file containing sampled sequences (.txt)



## Packages required

"pandoc" should be installed:

https://pandoc.org/installing.html

R packages:

rmarkdown, knitr, ape, igraph, reshape2 and phangorn
