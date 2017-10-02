# ScriptName: Gene_to_Sequences.pl
#
# Please direct questions and support issues to <https://hivdb.stanford.edu>
#
#
# Copyright 
#
# You may distribute this module under the same terms as perl itself
#
# _history
# September 21, 2017  Written by Soo-Yon Rhee

# POD documentation - main docs before the code

#!/usr/bin/perl -w


=head1 NAME

Gene_to_Sequences.pl - perform a local amino acid to nucelotide sequence search for a gene
                              and generate a file containing aligned full-length nucleotide sequences
                              with associated GB annotations including AccessionID, Title, Authors,
                              PubMedID, Country, Collection_Date, and TaxonomyID

=head1 SYNOPSIS

Use:

    perl Gene_to_Sequences.pl [--help] [--man] 
	                      --qAA path
	                      --qNA path
                              --M {tfastx or tblastn}
                              --blastDB path
                              --outFile path
                              --eval  real number
                              --percIden  integer 
                           

Examples:

    perl   Gene_to_Sequences.pl --help

    perl   Gene_to_Sequences.pl --man


    #    All required input parameters:

    perl   Gene_to_Sequences.pl --qAA  query.AA.fasta.txt
                                --qNA  query.NA.fasta.txt
                                --M tfastx
                                --blastDB  /data/genBankTest/ViralDB
                                --outFile results.txt


    #    All required and optional input parameters:

    perl   Gene_to_Sequences.pl --qAA  query.AA.fasta.txt
                                --qNA  query.NA.fasta.txt
                                --M tfastx
                                --blastDB  /data/genBankTest/ViralDB
                                --outFile results.txt
                                --eval    0.000001
                                --percIden  80
                                --verbose


=head1 DESCRIPTION

This script is the second part of the gene sequence search pipeline. 

The following steps are performed to generate aligned full-length nucleotide 
gene sequences with associated GB annotations including AccessionID, Title, 
Authors, PubMedID, Country, Collection_Date, and TaxonomyID:

1) search a local sequence database for gene-coding sequences using a local amino-acid
to nucleotide alignment search program, either tblastn or tfastx. The database should
be set up by the first part of the pipeine "GB_to_BLASTDB.pl" - this ensures
a fast retreival of GenBank annotations

2) for each sequence in the output from a tblastn or tfastx for which 
its % identity to a query amino acid sequence is greater than a given cutoff
and covers the gene in full-length, perform a global nucleotide-to-nucleotide
sequence alignment to generate an aligned nucleotide sequence - this ensures
any frameshifts get removed as sequences with frameshifts caused by sequencing errors 
are frequently found in a public database

3) for each of full-lengh sequences identified, extract associated GenBank annotations
and retrieve the Taxon Name by querying NCBI Taxonomy Database.


 
=head1 ARGUMENTS

perl Gene_to_Sequences.pl takes the following arguments:

=over 4

=item help

  --help

(Optional) Display the usage message.

=item man

  --man

(Optional) Displays all documentation.

=item qAA

  --qAA   path

(Required) Sets the path to the fasta file containing an amino acid sequence for search query

=item qNA

  --qNA   path

(Required) Sets the path to the fasta file containing a nucleotide sequence of the query amino acid sequence


=item M

   --M [tfastx|tblastn]

(Required) Sets the search method to use, tfastx or tblastn 


=item blastDB

  --blastDB  path

(Required) Sets the file path to a BLAST searchable sequence database that was formatted by NCBI blast+ program.

=item eval
  
  --eval  real number

(Optional) Sets the minimum e_value for the search. The default is 0.0000001.


=item  percIden

  --percIden  integer

(Optional) Sets the minimum % identity for considering an alignment.  the default is 80.

                                -
=item verbose

  --verbose

(Optional) Warning messages on


=back

=head1 AUTHOR

Soo-Yon Rhee, E<lt>syrhee@stanford.eduE<gt>.

=head1 COPYRIGHT

You may distribute this module under the same terms as perl itself

=head1 DATE

September 21, 2017


=cut

#' Let the code begin...

package main;
use strict;
use FileHandle;
use File::Temp;
use LWP::Simple;
use Getopt::Long();    
use Pod::Usage(); 
use Bio::SeqIO;
use Bio::SearchIO;

# Programs to be used. Set the proper path to each of the programs below 
my $tblastnp = "tblastn";
my $tfastxp = "tfastx36";
my $alnp = "ggsearch36";
my $blastdbcmd = "blastdbcmd";
my %searchp = ('tblastn' => $tblastnp,   'tfastx' => $tfastxp);

my $e_val= '0.0000001';
my $min_perc_similarity = 80;
my $allowedAAMismatchesAtOneEnd = 1;

# Set this number to ensure we get all hits when tblastn is used
# The default set by tblastn is only 500.
my $max_target_seqs = 1000000;       

## By default,  warning messages off
my $verbose = -1;


MAIN:
{
    my ($help, $man, $queryAAFastaFile, $queryNAFastaFile, $blastDB, $outFile, $searchM);
    Getopt::Long::GetOptions(
        'help'                =>  \$help,
        'man'                 =>  \$man,
	'qAA=s'               =>  \$queryAAFastaFile,
	'qNA=s'               =>  \$queryNAFastaFile,
	'M=s'                 =>  \$searchM,           
        'blastDB=s'           =>  \$blastDB,
        'outFile=s'           =>  \$outFile, 
	'eval=f'              =>  \$e_val,
	'percIden=i'          =>  \$min_perc_similarity,
        'verbose+'            =>  \$verbose 
    );
    
    
    # Check for requests for help or for man (full documentation):
    Pod::Usage::pod2usage( -verbose => 1 ) if ( $help );
    Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 2 ) if ( $man  );


    # Check for required variables.
    unless (defined($blastDB) && defined($outFile) && 
	    defined($queryAAFastaFile) && defined($queryNAFastaFile) && defined($searchM)){
	Pod::Usage::pod2usage("\nERROR: the required input parameters not set.\n");
    }

    
    # Check whether the search method is valid
    if(! exists $searchp{$searchM}){
	Pod::Usage::pod2usage("\nERROR: the search method is not valid.\n");
    }

    
    # Check the length of the query files
    if(! -e $queryAAFastaFile or ! -e $queryNAFastaFile){ 
	die "Can't find query files: $queryAAFastaFile and/or $queryNAFastaFile";
    }
    my $query_AA_seq_obj = Bio::SeqIO->new(-file =>$queryAAFastaFile,  -format=>'fasta')->next_seq;
    my $query_NA_seq_obj = Bio::SeqIO->new(-file =>$queryNAFastaFile,  -format=>'fasta')->next_seq;
    if(length($query_AA_seq_obj->seq) * 3 != length($query_NA_seq_obj->seq)){
	die "The length of query AA sequence and NA sequence do not match";
    }
    
    # Create some temporary files with UNLINK => 1 to ensure the temporary files get removed at the end 
    my $tempobj = File::Temp->new(TEMPLATE => 'tempXXXX', DIR => ".", SUFFIX => '.txt', UNLINK => 1);
    my $tempOutFile = $tempobj->filename;
    my $tempobj2 = File::Temp->new(TEMPLATE => 'tempXXXX', DIR => ".", SUFFIX => '.txt', UNLINK => 1);
    my $searchResultsFile = $tempobj2->filename;

    &Search($searchM, $queryAAFastaFile, $blastDB, $searchResultsFile);
    
    &Extract_gene_full_length_sequences($searchResultsFile, $queryNAFastaFile, $blastDB, $tempOutFile);
    
    &Taxon_name_search($tempOutFile, $outFile);


}





=head2 Search

  Title   : Search
  Usage   : Search($searchM, $query_AA_fasta_file, $blastDB, $searchResultsFile)
  Function: Perform a search using a given method and a query amino acid sequence
  Returns : 
  Args    : a search method - either tfastx or tblastn
	    a fast file containing a query amino acid sequence
	    a file path to a BLAST-searchable  DB
	    a path to a file for saving search output from the search program

=cut

sub Search{
    my ($searchM, $queryAAFastaFile, $blastDB, $searchResultsFile) = @_;
    
    my $cmd;
    if($searchM eq 'tfastx'){
	print "\nStarting TFASTAX search...\n";
	# To make FASTA programs recognize a BLAST-searchable database, just add the number 12
	# after the name of the BLAST DB
	# -m 8 ensures that the search results is generated in a tabular format.
    # By default, tfastx uses:
    #    -f (gapopen) -12, -g (gapextend) -2
    #
	# Explicitly show the default parameters in the command.
    #
	$cmd = "$searchp{$searchM} -E $e_val -m 8  -s BL62 -f -12  -g -2 $queryAAFastaFile  \"$blastDB 12\"";
    }
    elsif($searchM eq 'tblastn'){
	print "\nStarting TBLASTX search...\n";
	# -outfmt 6 (output format) ensures that the search results is generated in a tabuar format.
	# by default, tblastn uses BL62 for the matrix,  -gapopen 11, gapextend 1

	$cmd = "$searchp{$searchM}  -query $queryAAFastaFile -db $blastDB ".
	       "-outfmt 6 -max_target_seqs $max_target_seqs -evalue $e_val -max_hsps 1 -gapopen 11 -gapextend 1  -matrix BLOSUM62";
    }
    
    print "\t", $cmd, "\n";
    `$cmd > $searchResultsFile`;
    
    my $numHits = `wc -l $searchResultsFile`;
    $numHits =~ /^(\d+) /;   $numHits = $1;
    print "\tNumber of Hits: $numHits\n";
}



=head2 Extract_gene_full_length_sequences
    
    Title   : Extract_gene_full_length_sequences
    Usage   : Extract_gene_full_length_sequences ($searchResultsFile, $queryNAFastaFile, $blastDB, $outFile)
    Function: For each sequence in the output from a tblastn or tfastx search for which the % identity 
              is greater than a cutoff and covers the gene in full-length, generate a nucleotide sequence
              alignment by performing a global nucleotide-to-nucleotide sequence alignment program. This step
              corrects the sequence alignment in case there are one or more frameshifts in the raw nucleotide
              sequence.
    Example : 
    Returns : 
    Args    : a file containing search output in tabular format
	    a file containing a query nucleotide sequence for alignment
	    a blast-searchable database to retrieve raw nucleotide sequences
	    a file for storing full-length gene sequences with GB annotations

=cut

sub Extract_gene_full_length_sequences{
    my ($searchResultsFile, $queryNAFastaFile, $blastDB, $outFile) = @_;
    
    print "\nParsing search output and extracting complete gene sequences...\n";
    
    my $fh = new FileHandle(">$outFile") || die "cannot open the file $outFile";
    print $fh join "\t", qw(Accession eValue Score PubMedID Title Authors Journal
			     Country Note Isolate CollectionDate FirstAA LastAA AASeq NASeq TaxonID);
    print $fh "\n";
    
    # read in query nucleotide sequence for alignment
    my $query_NA_seq_obj = Bio::SeqIO->new(-file =>$queryNAFastaFile,  -format=>'fasta')->next_seq;
    my $queryAALen =  length($query_NA_seq_obj->seq) / 3;
    
    my $numCompleteGeneSequences = 0;
    my $in = new FileHandle("<$searchResultsFile") || die "cannot open the file $searchResultsFile";
    while(my $result = <$in>){
	$result =~ s/\n$//;
	my ($qname, $accession, $perc_id, $alignmentLen, $mismatches,
	    $gapopenings, $qstart, $qend, $sstart, $send, $eValue, $score) = 
		split "\t", $result;
	$accession =~ s/\w+\|\w+\|//;
	
	my ($FirstAA, $LastAA, $AASeq, $NASeq);
	if($send < $sstart){
	    print "Skipping negative strand:\n";
	    print $result, "\n";
	    next;
	}
	
	# Only consider an alignment of similarity > $min_perc_similarity
	# Local alignment programs commonly miss out consecutive mismatches at both ends
	# which are in fact mutations. Here, we allow $$allowedAAMismatchesAtOneEnd mutations
	# at each end and extend the alignment
	if( $perc_id > $min_perc_similarity &&
	    $alignmentLen >= $queryAALen - ($allowedAAMismatchesAtOneEnd *2)){
	    
	    ## retrieve the raw/original nucleotide sequence from the database
	    my $cmd = join " ", ($blastdbcmd, '-db', $blastDB, '-entry', $accession);
	    my $fastaNAseq = `$cmd`;
	    
	    my $org_sbj_na_obj = Bio::SeqIO->new(-string => "$fastaNAseq",  -format=>'fasta')->next_seq;
	    ($FirstAA, $LastAA, $NASeq, $AASeq) 
		= &get_alignment( $qstart, $qend, $sstart, $send,
				  $org_sbj_na_obj, $query_NA_seq_obj, $allowedAAMismatchesAtOneEnd);
	    if(length($AASeq) == $queryAALen){
		$numCompleteGeneSequences ++;
		## sequence header contains GenBank annotations
		my $seqHeader  = (split "\n", $fastaNAseq)[0];
		$seqHeader =~ s/^([\S]+ )//;   $seqHeader =~ s/\n$//;

		my ($pubmed, $title, $authors, $journal, $country, $note, $isolate, $collection_date, $taxonID) = 
		    split /%%/, $seqHeader;
		print $fh join "\t", ($accession, $eValue, $score, $pubmed, $title, $authors, $journal, 
				      $country, $note, $isolate, $collection_date, $FirstAA, $LastAA, $AASeq, $NASeq, $taxonID);
		print $fh "\n";
	    }
	}
    }
    
    $fh->close;
    $in->close;
    print "\t$numCompleteGeneSequences Complete gene sequences retreived\n"; 
}


=head2 get_alignment

    Title   : get_alignment
    Usage   : 
    Function: Given the starting and ending amino acid positions of a query sequence
              and the starting and ending nucleotide positions of a subject sequence 
              of HSP, it performs a global nucleotide-to-nucleotide alignment between
              the query and the subject sequence to correct the alignment in case
              one and more frameshifts are present in the subject nucleotide sequence.

    Returns : First amino acid postion of an aligned sequence relative to a query amino acid sequence
              Last amino acid postion of an aligned sequence relative to a query amino acid sequence
              Aligned nucleotide sequence
              Translation of the aligned nucleotide sequence

    Args    : First amino acid position of a query in the alignment of HSP reported by tblastn or tfastx program 
              Last amino acid position of a query in the alignment of HSP reported by tblastn or tfastx program 
              First nucleotide position of a subject in the alignment of HSP reported by tblastn or tfastx program 
              Last nucleotide position of a subject in the alignment of HSP reported by tblastn or tfastx program 
              Reference to a Bio::Seq object representing a raw-uncut-unaligned nucelotide sequence of the subject
              Reference to a Bio::Seq object representing a nucleotide sequence of the query
              Number of amino acids mismatches allowed at each end of the alignment of HSP

=cut

sub get_alignment{
    my ($FirstAA, $LastAA, $NAStart, $NAEnd, $org_sbj_na_obj, $query_NA_seq_obj, $allowedAAMismatchesAtOneEnd) = @_;
    
    my ($AASeq, $NASeq);
    my $queryAALen =  length($query_NA_seq_obj->seq) /3;
    
    # Local alignment programs could easily miss a few amino acids at the ends of their alignment
    # when those amino acids are mismatches. For a local alignment that starts or ends at
    # +/- $allowedAAMismatchesAtOneEnd, the alignment gets extended by $allowedAAMismatchesAtOneEnd at an end
    # if the subject nucleotide sequence encompasses those positons.

    if(($FirstAA == ($allowedAAMismatchesAtOneEnd + 1)) and ($NAStart  >= ($allowedAAMismatchesAtOneEnd * 3) + 1)){
	$FirstAA = 1;              $NAStart = $NAStart - 3;
    }    
    if(($LastAA == $queryAALen - $allowedAAMismatchesAtOneEnd) and ($NAEnd + 3  <= length($org_sbj_na_obj->seq))){
	$LastAA = $queryAALen;     $NAEnd = $NAEnd + 3;
    }
    
    $NASeq = &run_global_aln($query_NA_seq_obj->subseq((($FirstAA-1)*3)+1, $LastAA*3), 
			     $org_sbj_na_obj->subseq($NAStart, $NAEnd));
    
    if(($LastAA - $FirstAA + 1) *3 != length($NASeq) and length($AASeq)*3 != length($NASeq)){
	die "Wrong alignment\n";
    }
    # translate the NASeq
    $AASeq = new Bio::Seq->new(-display_id => 'subject', -seq=>$NASeq)->translate->seq;

    return ($FirstAA, $LastAA, $NASeq, $AASeq);
}


=head2 run_global_aln
    
  Title   : run_global_aln
  Usage   : 
  Function: Performs a global alignment of the two nucleotide sequences
            and removes frameshifts.

  Returns : Aligned nucleotide sequence
  Args    : a query NA sequence in a string
	    a subject NA sequence in a string

=cut

sub run_global_aln{
    my ($query_NA_seq,  $sbj_NA_seq) = @_;
    
    ## the global alignment program, ggsearch36, takes an input sequence from a file

    my $tempobj1 = File::Temp->new(TEMPLATE => 'tempXXXX', DIR => ".", SUFFIX => '.txt', UNLINK => 1);
    my $temp_q = $tempobj1->filename;
    my $tempobj2 = File::Temp->new(TEMPLATE => 'tempXXXX', DIR => ".", SUFFIX => '.txt', UNLINK => 1);
    my $temp_s = $tempobj2->filename;

    my $out_query_NA = Bio::SeqIO->new(-file=> ">$temp_q", -format=>'fasta');
    if(! $out_query_NA->write_seq(Bio::Seq->new(-display_id => "pol", -seq=> $query_NA_seq))){
	die( "Problem with creating a sequence file: $temp_q");
    }    
    my $out_subject_NA = Bio::SeqIO->new(-file=> ">$temp_s", -format=>'fasta');
    if(! $out_subject_NA->write_seq(Bio::Seq->new(-display_id => "subject", -seq=> $sbj_NA_seq))){
	die( "Problem with creating  a sequence file: $temp_s");
    }   
    
    my $cmd = "$alnp $temp_q  $temp_s";
    
    my $searchio_lalign = Bio::SearchIO->new(-format => 'fasta',
					     -file   => "$cmd |",
					     -report_type => 'lalign');
    
    my $hsp_nn = $searchio_lalign->next_result->next_hit->next_hsp;
    my $sbj_nn = $hsp_nn->seq('sbjct');
    my $query_nn = $hsp_nn->seq('query');
    my $NASeq = '';
    if($query_nn->seq !~ /\-/){
	# No insertion or frameshits
	$NASeq = $sbj_nn->seq;
    }
    else{
	# there are some insertions/frameshits
	# remove them
	my @NAs_q = split '', $query_nn->seq;
	my @NAs_s = split '', $sbj_nn->seq;
	my @NAs_new = ();
	foreach (my $i = 0; $i <= $#NAs_q; $i=$i +1){
	    if($NAs_q[$i] ne '-'){
		push @NAs_new, $NAs_s[$i];
	    }
	}
	$NASeq = join "", @NAs_new;
    }

    return $NASeq;
}


=head2 Taxon_name_search

  Title   : Taxon_name_search
  Usage   : 
  Function: Query NCBI Taxonomy Database to retrieve a corresponding
            TaxonName of a TaxonID
  Returns : 
  Args    : File containing search results containing TaxonID 
            File for saving search results with additinal column - TaxonName

=cut

sub Taxon_name_search{
    my ($tempOutFile, $outFile) = @_;
    
    print "\nRetreiving Taxon Name from NCBI Taxonomy DB...\n";
    
    my $fh = new FileHandle(">$outFile") || die "cannot open the file $outFile";

    my $in = new FileHandle("<$tempOutFile") || die "cannot open the file $tempOutFile";
    my $colmNameLine = <$in>;    $colmNameLine =~ s/\n$//;
    print $fh join "\t", ($colmNameLine, "TaxonName"); 
    print $fh "\n";
    
 
    my %taxhash = ();
    while(my $row = <$in>){
	$row =~ s/\n$//;
	my $taxonID = (split "\t", $row) [15];
	$taxhash{$taxonID} = 1;
    }
    $in->close;
    
    ## query NCBI::Taxonomy database to retreive Taxon Name
    foreach my $taxonID (keys %taxhash){
	my $url = "https://www.ncbi.nlm.nih.gov/taxonomy/?term=". $taxonID . "&report=taxon&format=text";
	my $content = get($url);
	$content =~ s/\n//g;  $content =~ /<pre>(.+)<\/pre>/; 
	$taxhash{$taxonID} = $1;
    }  
    
    
    $in = new FileHandle("<$tempOutFile") || die "cannot open the file $tempOutFile";
    <$in>;                                ## skip the column header line
    while(my $row = <$in>){
	$row =~ s/\n$//;
	my @values = split "\t", $row;
	my $taxonID = $values[15];
	print $fh join "\t", (@values, $taxhash{$taxonID});
	print $fh "\n";
	
    }
    $in->close;
    $fh->close;
    
    my $numHits = `wc -l $outFile`;
    $numHits =~ s/\D+$//;
    $numHits --;
    print "\n$numHits Complete gene sequences with annotations generated in $outFile\n\n";
}
