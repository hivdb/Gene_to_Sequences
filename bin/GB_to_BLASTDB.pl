# ScriptName: GB_to_BLASTDB.pl
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

GB_to_BLASTDB.pl - parse GenBank files, create a fasta file with sequence headers 
                   containing GenBank annotations and convert the fasta file 
                   to a BLAST searchable database

=head1 SYNOPSIS

Use:
    
    perl GB_to_BLASTDB.pl [--help] [--man] 
                          --gbDir path
                          --blastDB path
                          --verbose 

Examples:

    perl GB_to_BLASTDB.pl --help

    perl GB_to_BLASTDB.pl --man

    perl GB_to_BLASTDB.pl --gbDir /data/genBank/UnZipped
                          --blastDB  /data/genBank/ViralDB
                          --verbose


=head1 DESCRIPTION

This script is the first part of the gene search pipeline. It creates 
a BLAST-searchable database containing all GenBank virus sequences
and GenBank annotations. Next, gene-coding sequences are retrieved 
using "Gene_to_Sequences.pl".

Before running this script, download all GenBank files of virus sequences,
unzipped them and store in a directory.

Detailed information about GenBank files can be found at:
  https://www.ncbi.nlm.nih.gov/genbank/
 
GenBank files can be downloaded from ftp://ftp.ncbi.nih.gov/genbank/
To download viral sequences only, run the command below on the ftp site:
>mget gbvrl*

The steps in this script:

1) parse GenBank files, extract GenBank annotations and create
a fasta sequence file such that accession numbers become the sequence
identifiers (which can be used to retrieve sequences using NCBI blast+ 
package programs) and extracted GenBank annotations become a short description
in the sequence header.  

Note: BioPerl does an excellent job but it seems very slow. Boulder::Genbank
library is a lot faster for parsing GenBank files, but the library is rarely maintained.
When NCBI changes the format of their files, it runs into problems. However, if speed is
a concern and you are willing to update the library yourself, the Boulder::Genbank is highly
recommanded.  Boulder::Genbank can be found at:
 http://search.cpan.org/~lds/Boulder-1.30/Boulder/Genbank.pm

2) convert the fasta sequence file to a BLAST-searchable database
using the NCBI blast+ program. A database formatted by NCBI blast+
is also recognized by other search programs, ie. FASTA package


=head1 ARGUMENTS

perl GB_to_BLASTDB.pl takes the following arguments:

=over 4

=item help

  --help

(Optional) Display the usage message.

=item man

  --man

(Optional) Displays all documentation.

=item gbDir

  --gbDir path

(Required) Sets the path to the directory in which the GenBank unzipped flat files are stored.

=item blastDB

  --blastDB  path

(Required) Sets the file path to BLAST searchable database.

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
use Getopt::Long();    
use Pod::Usage();     	
use Bio::SeqIO;

# makeblastdb is the Blast program for converting a fasta sequence file
# to a BLAST searchable database. makeblastdb is included in the NCBI blast+ package. 
my $formatp =  "makeblastdb";

## By default,  warning messages off
my $verbose = -1;

MAIN:
{

    my ($help, $man, $gbRecFilesDir, $blastDB);
    Getopt::Long::GetOptions(
        'help'              =>  \$help,
        'man'               =>  \$man,
        'gbDir=s'           =>  \$gbRecFilesDir,
        'blastDB=s'         =>  \$blastDB,
	'verbose+'          =>  \$verbose);
 
        
    # Check for requests for help or for man (full documentation):
    Pod::Usage::pod2usage( -verbose => 1 ) if ( $help );
    Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 2 ) if ( $man  );

    # Check for required variables.
    unless ( defined( $gbRecFilesDir ) && defined( $blastDB ) ){
	Pod::Usage::pod2usage( -exitstatus => 2 );
    }
 
    # Generate a unique temporary file name for storing fasta sequences
    # UNLINK =>1 ensures the temporary file gets removed at the end.
    my $tempobj = File::Temp->new(TEMPLATE => 'tempXXXX', DIR => ".", SUFFIX => '.txt', UNLINK => 1);
    my $fastaFile = $tempobj->filename;

    my ($numFiles, $n) = &Create_fasta_file($gbRecFilesDir, $fastaFile);
    print "$numFiles GenBank files processed\n";
    print "$n Accessions/sequences parsed/added to a fasta file\n";

    if($n == 0 or $numFiles == 0){
	die("No GB sequences found");
    }
    
    &Format_DB($fastaFile, $blastDB);
    print "$blastDB created\n";

}





=head2 Create_fasta_file

 Title   : Create_fasta_file
 Usage   : Create_fasta_file($gbRecFilesDir, $fastaFile)
 Function: Parse all *.seq files and create a fasta file
 	   with sequence headers containing GenBank annotations
 Example : 
 Returns : Number of GenBank files and sequences processed
 Args    : Directory where unzipped GenBank files are stored
           File to store fasta formatted sequences 
=cut

sub Create_fasta_file{
    my ($gbRecFilesDir, $fastaFile) = @_;
    my $n = 0;  my $numFiles = 0;
    my @gbFiles = ();
    
    opendir(DIR, $gbRecFilesDir) or die "Can't read the directory $gbRecFilesDir: $!";
    while(defined(my $file = readdir(DIR))){
	if($file =~ /\.seq$/){push @gbFiles,$file;}
    }

    # open a output file for fasta sequences
    my $fastafh = new FileHandle(">$fastaFile")||die ("Can not open the file $fastaFile: $!");  

    foreach my $gbFile (@gbFiles){
	my $file = join "/", ($gbRecFilesDir, $gbFile);   
	print "\tProcessing $file ....\n";
	
	# read in GenBank file format
	my $io = Bio::SeqIO -> new(-file => $file, -format => "genbank", -verbose => $verbose);
	$numFiles ++;
	while(my $seqobj = $io->next_seq()){
	    my $accession = $seqobj->accession;
	    my $annoobj = $seqobj->annotation;
	    
	    # use the first reference which is the lastest published
	    my $reference = ($annoobj->get_Annotations('reference'))[0];  
	    my $title = $reference->title() || '';
	    my $authors = $reference->authors() || '';
	    my $journal = $reference->location() || '';
	    my $pubmed = $reference->pubmed() || '';
	    
	    my $featobj = ($seqobj->get_SeqFeatures)[0];
	    my $isolate = $featobj->has_tag('isolate') ? ($featobj->get_tag_values('isolate'))[0] : '';
	    my $country = $featobj->has_tag('country') ? ($featobj->get_tag_values('country'))[0] : '';
	    my $collection_date =  $featobj->has_tag('collection_date') ? ($featobj->get_tag_values('collection_date'))[0] : '';
	    my $note = $featobj->has_tag('note') ? ($featobj->get_tag_values('note'))[0] : '';
	    my $taxonID = $featobj->has_tag('db_xref') ? ($featobj->get_tag_values('db_xref'))[0] : '';
	    $taxonID =~ s/taxon://;
	    
	    # include extracted GB annotations to a sequence header (delimited by %%)
	    # and add sequence to Fasta file
	    print $fastafh ">gnl|MYDB|$accession ";
	    print $fastafh join "%%",($pubmed, $title, $authors, $journal, $country, $note, $isolate, $collection_date, $taxonID); 
	    print $fastafh "\n";
	    print $fastafh  $seqobj->seq, "\n";	    
	    $n ++;
	}
    }
    $fastafh->close;
    return ($numFiles, $n);
}


=head2 Format_DB

 Title   : Format_DB
 Usage   : Format_DB ($fastaFile, $blastDB)
 Function: Convert a fasta file to a BLAST-searchable database
 Example : 
 Returns : 
 Args    : File containing sequences in fasta format
           File path for BLAST-searchable database to be created

=cut

sub Format_DB{
    my ($fastaFile, $blastDB) = @_;
    
    my $cmd = "$formatp -in $fastaFile -dbtype nucl -parse_seqids -out $blastDB";
    print "Converting fasta file to BLAST database ....\n";
    print "\t", $cmd, "\n";
    system($cmd);
}
