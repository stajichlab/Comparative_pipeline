#!/usr/bin/perl

=head1 USAGE

perl scripts/orthomcl_phylogenetic_profile.pl -p Domains  -i OrthoMCL.out --cf collapsefile.txt

=head1 DESCRIPTION

Makes an ortholog phylogenetic distribution file

=head1 EXAMPLE

collapsefile.txt
Basidiobolus	Basme2finSC,N161,N168
ENT	Conco1,Conth1

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Spec;

my $input;
my $db;

my $odir;
my $pfamdir;
my $pfamext = 'domtbl';
my $hmmerversion = 3;
my $evalue_cutoff = 0.01;
my @collapse;
my $collapse_file;
GetOptions(
	   'd|p|pfam:s'  => \$pfamdir,
	   'ext:s'     => \$pfamext,
	   'hmmer|hv:s'=> \$hmmerversion,
	   'e|evalue:s' => \$evalue_cutoff,
	   'i|input:s'  => \$input,
	   'o|output:s' => \$odir,
	   'c|collapse:s{2}' => \@collapse,
	   'cf|collfile:s' => \$collapse_file,
	   );
my %collapse = @collapse;

$input ||= shift @ARGV;
$odir ||= $input."-profile";

if( $collapse_file ) {
  open(my $cfile => $collapse_file) || die $!;
  while(<$cfile>) {
   next if /^\#/ || /^\s+$/;
   my ($clade,$species) = split;
   for my $m ( split(/,/,$species) ) {
     $collapse{$m} = $clade;
   }
  }
}
die "must provide a folder with pfam table data (either domtblout or hmmer2table output\n" unless $pfamdir && -d $pfamdir;
my %pfam;
opendir(PFAM, $pfamdir) || die "cannot open $pfamdir: $!";
for my $file ( readdir(PFAM) ) {
    next unless ( $file =~ /\.\Q$pfamext\E$/);
    open(my $fh => "$pfamdir/$file" ) || die $!;
    if( $hmmerversion == 3 ) {	#parse HMMER3 domtblout output
	while(<$fh>) {
	    next if /^\#/;
	    my ($domain,$acesssion,$tlen, $qname, $qacc,
		$qlen, $seq_evalue, $seq_score,$seq_bias,
		$n, $totalhits, $dom_cvalue, $dom_ievalue,$dom_score,
		$dom_bias) = split(/\s+/,$_);
	    $pfam{$qname}->{$domain}++ if $dom_cvalue < $evalue_cutoff;
	}
    } elsif( $hmmerversion ==2 ) { # parse HMMER2 (hmmer_to_table output
	while(<$fh>) {
	    next if /^\#/;
	    my ($qname, $qstart, $qend, $domain,
		$domstart, $domend, $score, $evalue) = split;	       
	    $pfam{$qname}->{$domain}++ if $evalue < $evalue_cutoff;
	}
    } else { 
	warn("unknown HMMER version (2 or 3 are expected)\n");
    }
}
open(my $genefh =>">$odir.pergene.tab");
open(my $statsfh =>">$odir.clade_stats.tab");
print $statsfh join("\t", qw(ORTHOGROUP_COUNT CLADES)),"\n";
open(my $fh => $input) || die "cannot open $input: $!\n";
my (%sp,%group2name,%clades);
my $groupn = 0;
while(<$fh>) {   
    my (@orthologs) = split;
    my $group = sprintf("ORTHO_%05d",$groupn++);
    my %count;
    my %domains;	
    for my $gene ( @orthologs ) {
	my ($sp)= split(/\|/,$gene);
	if( $collapse{$sp} ) {
	    $sp = $collapse{$sp};
	}
	$count{$sp}++;
	$sp{$sp} = 1; # collecting all the species names
	push @{$group2name{$group}}, $gene;
    }
    my $nm = join("-", sort keys %count);
    push @{$clades{$nm}}, $group;
}
mkdir($odir) unless -d $odir;
print $genefh join("\t",qw(GENE CLADE ORTHOGROUP DOMAINS)),"\n"; 
my %stats;

for my $clade ( 
    # schwartzian transformation
    map { $_->[1] } 
    sort { $b->[0] <=> $a->[0] } 
    map { [ scalar @{$clades{$_}}, $_] }
# schwartzian transformation here to process clades by largest 
# membership to smallest
    keys %clades ) {
    my $og = $clades{$clade};
    print $statsfh join("\t",(scalar @$og), $clade),"\n";

    open(my $ofh => ">$odir/$clade.phyloprofile.tab") || die $!;
    print $ofh join("\t", 'ORTHOGROUP','GENE', 'DOMAINS'), "\n";
    for my $ortho ( @$og ) {
	for my $gene ( sort @{$group2name{$ortho}} ) {	    
	    print $ofh join("\t", $ortho, $gene, 
			    join(",", keys %{$pfam{$gene} || {}})),"\n";
	    print $genefh join("\t", $gene, $clade, $ortho, join(",", keys %{$pfam{$gene} || {}})),"\n";
	}
    }
}
