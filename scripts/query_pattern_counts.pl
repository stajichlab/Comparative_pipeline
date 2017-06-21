#!/usr/bin/perl
use strict;
use warnings;

=head1 NAME

query_pattern_counts.pl

Query a tab delimited table of gene counts from orthomcl / MCL
output, with an appended Pfam domain table summary

=head1 USAGE

perl Comparative_pipeline/scripts/query_pattern_counts.pl \
  -i Orthologs/Basidiobolus.OrthoMCL.I15.counts.tab 
  --present Basme2finSC,N161,N168
  --zero REMAINDER
  --ignore Pipcy3_1

=head1 AUTHOR

Jason Stajich @hyphaltip http//github.com/hyphaltip

=cut

use File::Spec;
use Getopt::Long;

my ($input,$ignore,$present,$zero) = (undef,'','','');

GetOptions(
    'i|input:s'         => \$input,
    'ignore:s'          => \$ignore,
    'p|present:s'       => \$present,
    'z|zero:s'          => \$zero);

$input = shift @ARGV if ! $input && @ARGV;
my $fh;
if( $input ) {
    open($fh => $input) || die $!;
} else {
    $fh = \*STDIN;
}

my $header = <$fh>;
print $header;
chomp($header);
my ($orthid,@species_pref,$domain) = split(/\t/,$header);

my %keep;
my %zeroes;
my %ignores;
my %sp2info = map { $_ => '' } @species_pref;

for my $n ( split(/,/,$present)) {
    $sp2info{$n} = 'present';
    $keep{$n}++;
}
warn(join(" ", keys %keep),"\n");
for my $n ( split(/,/,$ignore)) { 
    if( $n =~ /REMAIN/ ) {
	for my $nm ( grep { $sp2info{$_} eq '' } @species_pref ) {
	    $sp2info{$nm} = 'ignore';
	    $zeroes{$nm}++;
	}
    } else {
	$ignores{$n}++;
	$sp2info{$n} = 'ignore';
    } 
}
for my $n ( split(/,/,$zero )) {
    if ( $n =~ /REMAIN/ ) {
	for my $nm ( grep { $sp2info{$_} eq '' } @species_pref ) {
	    $sp2info{$nm} = 'zero';
	    $zeroes{$nm}++;
	}
	last;
    } else { 
	$zeroes{$n}++;
	$sp2info{$n} = 'zero'; 
    }
}
while(<$fh>) {
    chomp;
    my ($id,@cols) = split(/\t/,$_);
    my $domains = pop @cols;
    my $i = 0;
    my $ok = 0;
    for my $col ( @cols ) {
	my $sp = $species_pref[$i++];
	if( $sp2info{$sp} eq 'present' && $col > 0 ) { 
	    $ok = 1
	} elsif ( $sp2info{$sp} eq 'present' && $col == 0 ) {
	    $ok = 0;
	    last;
        } elsif ( $sp2info{$sp} eq 'zero' && $col == 0 ) { 
	    $ok = 1;
        } elsif ( $sp2info{$sp} eq 'zero' && $col != 0 ) { 
	    $ok = 0;
	    last;
	}
    }
    print join("\t", $id,@cols,$domains),"\n" if $ok;
}

