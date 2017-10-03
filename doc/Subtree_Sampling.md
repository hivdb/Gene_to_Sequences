# Subtree_Sampling.Rmd

Sampling subtype reference sequences from subtrees of a phylogenetic tree


## Description
This R script has been described in the manuscript by Rhee, S-Y and Shafer, RW, "Geographically-Stratified Representative HIV-1 Group M pol Subtype and Circulating Recombinant Form Sequences". For each fasta sequence file in a given directory, this script performs the followings: (1) calculates TN93 pairwise distances (PWDs);  (2) creates a neighbor joining tree with the PWD matrix and roots it at the midpoint; (3) traverses the tree from the root to leaves and partitioned it into subtrees having a local median PWD (MPWD) less than a cutoff; (4) then selects one sequence per country from each subtree. In the manuscript, 47% (0.47) of the ranked PWDs was used to define the maximum MPWD for subtrees.



## Usage

To run the script using the default parameters

Rscript -e "rmarkdown::render('Subtree_Sampling.Rmd')"  

Rscript -e "rmarkdown::render('Subtree_Sampling.Rmd')" fastadir="in.fasta", outdir="out.trees", pwd_rank_cutff=0.47



## Input and the default values

1. fastadir ("in.fasta" by default): directory where fasta sequence files are located. The name of each fasta file consists of subtype + ".txt". ie, 01_AE.txt, B.txt. The name of each sequence in a fasta file consists of sequenceid (in this case, accessionid), country (in two letter ISO code), collection year (YY), subtype delimited by "."

2. outdir ("out.trees" by default): directory for saving output files including .tre, .pdf, and .txt of sampled sequences

3. pwd_rank_cutff(0.47 by default): a node of which a local MPWD (MPWD of the node's children) is less than the global MPWD at the rank of this cutoff defined as a cluster (subtree)




## Output

1. Subtree_Sampling.html

2. results.summary.txt: lists number of clusters, number of sequences sampled so on.

Files generated in outdir for each subtype:

3. a tree of the entire sequences (.tre) (.pdf)

4. a tree of the entire sequences with cluster ids defined by a cutoff (.tre)

5. a fasta sequence file containing sampled sequences (.txt)

## Example fasta files are available at 


## Packages required

"pandoc" should be installed:

https://pandoc.org/installing.html

R packages:

rmarkdown, knitr, ape, igraph, reshape2 and phangorn
