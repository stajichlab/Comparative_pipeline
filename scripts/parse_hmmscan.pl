#!env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);

my $cutoff = 1e-3;
my $debug = 0;
GetOptions('c|e|cutof f:s' => \$cutoff,
    'v|verbose!' => \$debug);

mkdir("summary");
my %counts;
my %counts_gene;
my %taxa;
my %domain2acc;
while(<>) {
    next if /^\#/;
    my @row = split(/\s+/,$_);
    my $domain = $row[0];
    my $acc    = $row[1];
    my $gene = $row[3];
    my $evalue = $row[6];
    warn("domain=$domain gene=$gene evalue=$evalue\n") if $debug;
#    next if $evalue > $cutoff;
    my ($tax,$gn) = split(/\|/,$gene);
    $taxa{$tax}++;
    $counts{$domain}->{$tax}++;
    $counts_gene{$domain}->{$tax}->{$gn}++;
    $domain2acc{$domain} = $acc;
}

my @taxanames = sort keys %taxa;
open(my $fh => ">summary/Pfam_counts.tsv") || die $!;
open(my $fhgn => ">summary/Pfam_counts_genes.tsv") || die $!;

print $fh join("\t", qw(DOMAIN ACCESSION), @taxanames), "\n";
print $fhgn join("\t", qw(DOMAIN ACCESSION), @taxanames), "\n";

for my $s ( map { $_->[0] }
	    sort { $b->[1] <=> $a->[1] }
	    map { [$_, sum(values %{$counts{$_}})] } keys %counts ) {

    print $fh join("\t",$s,$domain2acc{$s},map { $counts{$s}->{$_} || 0 } @taxanames), "\n";
    print $fhgn join("\t",$s,$domain2acc{$s},
		     map { scalar keys %{$counts_gene{$s}->{$_} || {}} } 
		     @taxanames), "\n";
}
