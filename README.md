# Gene_to_Sequences

Gene search pipeline for retrieving gene-coding nucelotide sequences and annotations from a local GenBank database. 

GB_to_BLASTDB.pl 
parses GenBank files, creates a fasta file with sequence headers containing GenBank annotations and
converts the fasta file to a BLAST searchable databae.

Gene_to_Sequences.pl
as the second part of the pipeline, it performs an amino acid to nucleotide sequence search for a gene
and generates a file containing aligned full-length nucleotide sequences with associated GenBank annotations
including AccessionID, Title, Authors, PubMedID, Country, Collection_Date and TaxonomyID



## Prerequisites

* NCBI blast+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* BioPerl (http://bioperl.org/)
* FASTA package (http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml)

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
 




