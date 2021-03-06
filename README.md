# Gene_to_Sequences

Gene search pipeline for retrieving gene-coding nucelotide sequences and annotations from a local GenBank database.


This package is described in Rhee, S-Y and Shafer, RW (2017), "Geographically-Stratified Representative HIV-1 Group M pol Subtype and Circulating Recombinant Form Sequences," manuscript submitted for publication.


The package includes:

1. GB_to_BLASTDB.pl (a perl script)
parses GenBank files, creates a fasta file with sequence headers containing GenBank annotations and
converts the fasta file to a BLAST searchable databae.

2. Gene_to_Sequences.pl (a perl script)
as the second part of the pipeline, it performs an amino acid to nucleotide sequence search for a gene
and generates a file containing aligned full-length nucleotide sequences with associated GenBank annotations
including AccessionID, Title, Authors, PubMedID, Country, Collection_Date and TaxonomyID

3. Subtree_Sampling.Rmd (a R markdown)
samples subtype reference sequences from subtrees of a phylogenetic tree  
[Usage guidelines for this R markdown](doc/Subtree_Sampling.md)


## Prerequisites

* [NCBI blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [BioPerl](http://bioperl.org/)
* [FASTA](http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml)

## Usage

For the usage
```
perl GB_to_BLASTDB.pl --help

perl Gene_to_Sequences.pl --help

```

For more information
```
perl GB_to_BLASTDB.pl --man

perl Gene_to_Sequences.pl --man

``` 


## Authors

* **Soo-Yon Rhee** - *Initial work* - [sooyonrhee](https://github.com/sooyonrhee)



## License

You may distribute this module under the same terms as perl itself
 




