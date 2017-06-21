#!env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);

use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
use Bio::SeqIO;

my $man = 0;
my $help = 0;
my ($indir,$domaindir) = qw(results domainseq);
my $out; # can be empty will print to stdout
my $ext = '.domtbl';
my $sep = '__';
my $cutoff = 1e-2;
my $seqdb;
my $speciesorder;
my $verbose;

GetOptions('help|?'               => \$help, man => \$man,
	   'i|input:s'            => \$indir,
           'o|out:s'              => \$out,
	   'v|verbose!'           => \$verbose,
	   'd|domain|domainout:s' => \$domaindir,
	   'c|cutoff|evalue:s'    => \$cutoff,
	   'ext|extension:s'      => \$ext,
	   'db|database|seqs:s' => \$seqdb,
	   'v|verbose!' => \$debug,	   
	   'species:s' => \$speciesorder,
) or pod2usage(2);
pod2usage(1) if $help;

mkdir($domaindir) unless -d $domaindir;
my %counts;
my %counts_gene;
my %taxa;
while(<>) {
    next if /^\#/;
    my @row = split(/\s+/,$_);
    my $domain = $row[0];
    my $gene = $row[3];
    my $evalue = $row[6];
    warn("domain=$domain gene=$gene evalue=$evalue\n") if $debug;
#    next if $evalue > $cutoff;
    my ($tax,$gn) = split(/\|/,$gene);
    $taxa{$tax}++;
    $counts{$domain}->{$tax}++;
    $counts_gene{$domain}->{$tax}->{$gn}++;
}

my @taxanames = sort keys %taxa;
open(my $fh => ">summary/Pfam_counts.tsv") || die $!;
open(my $fhgn => ">summary/Pfam_counts_genes.tsv") || die $!;

print $fh join("\t", qw(DOMAIN), @taxanames), "\n";
print $fhgn join("\t", qw(DOMAIN), @taxanames), "\n";

for my $s ( map { $_->[0] }
	    sort { $b->[1] <=> $a->[1] }
	    map { [$_, sum(values %{$counts{$_}})] } keys %counts ) {

    print $fh join("\t",$s,map { $counts{$s}->{$_} || 0 } @taxanames), "\n";
    print $fhgn join("\t",$s,
		     map { scalar keys %{$counts_gene{$s}->{$_} || {}} } 
		     @taxanames), "\n";
}


=head1 NAME

gather_hmmscan_domain_counts.pl - generate a table of domain counts from hmmscan of Pfam or other DB

=head1 USAGE

perl gather_hmmscan_domain_counts.pl -i domains -o table [-ext .domtbl]

=head1 DESCRIPTION

Generate table of domain counts for use in comparing, plotting as heatmap, etc

=head2 ARGUMENTS

-i / --input    Input folder  [One file per taxon/proteome]
-o / --output   Output table
-db /--database folder or file of sequences
-ext            Extension of the result files. Default is .domtbl
-c / --cutoff   evalue cutoff

=head1 AUTHOR

Jason Stajich - jasonstajich.phd[at]gmail.com
University of California-Riverside

=cut
