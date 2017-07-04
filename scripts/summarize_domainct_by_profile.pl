#!/usr/bin/perl

=head1 summarize_domainct_by_profile

Used to compare domain counts between phylogenomic clades first built
out of ortholog clusters.

perl summarize_domainct_by_profile.pl -i PROJECT.pergene.tab

=cut

use strict;
use warnings;

use Getopt::Long;
use List::Util qw(sum min max);
my $pergene;
my $outdir = 'domain_counts';
GetOptions(
    'i|input:s' => \$pergene,
    'o|out:s'   => \$outdir,
    );

if ( ! -d $outdir ) {
    mkdir($outdir);
}

if( ! $pergene ) { 
    $pergene = shift @ARGV if @ARGV;
}
if( ! $pergene ) {
    opendir(DIR,".") || die $!;
    for my $file ( readdir(DIR) ) {
	next unless $file =~ /(\S+)\.pergene\.tab$/;
	$pergene = $file;
	last; # stop after the first one
    }
    if( ! $pergene ) {
	die"need a pergene file with -i or on cmdline or in the current folder\n";
    }
}

open(my $in => $pergene) || die $!;
my $header = <$in>;
my (%domain_ct,%orth);
while(<$in>) {
    my ($gene,$clade,$og, $domains) = split;
    for my $d ( split(/,/,$domains) ) {
	$domain_ct{$clade}->{$d}->{$og}++;
	$orth{$clade}->{$og}++;
    }
}

for my $cl ( keys %domain_ct ) {
    open(my $out => ">$outdir/$cl.domain_count") || die $!;
    print $out join("\t", qw(COUNT NORM_COUNT DOMAIN)),"\n";
    open(my $outct => ">$outdir/$cl.ortho_count") || die $!;
    my $OG_count =  scalar keys %{$orth{$cl}};
    printf $outct "%s,%d orthogroups,%d genes\n",$cl, 
    $OG_count, sum(values %{$orth{$cl}});
    for my $dm ( 
	sort { $b->[0] <=> $a->[0] }
	map { [ scalar keys %{$domain_ct{$cl}->{$_}}, $_] }
	keys %{$domain_ct{$cl}} ) {
	print $out join("\t", $dm->[0],
			sprintf("%.2f",100 * ($dm->[0] / $OG_count)),
			$dm->[1]), "\n";
    }
}

