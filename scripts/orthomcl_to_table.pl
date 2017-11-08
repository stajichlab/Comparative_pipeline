#!env perl

=head1 NAME

orthomcl_to_table.pl -

Generate a tab delimited table of gene counts from orthomcl / MCL
output, with an appended Pfam domain table summary

=head1 USAGE

perl Comparative_pipeline/scripts/orthomcl_to_table.pl
  -r Orthologs/Basidiobolus.OrthoMCL.I15.counts.tab 
  -i Orthologs/Basidiobolus.OrthoMCL.I15.out 
  -db Orthologs/goodProteins.fasta 
  --pfam Domains/Pfam 
  --cazy Domains/CAZY
  --merops Domains/MEROPS
  -s Basme2finSC,N161,N168,Conco1,Conth1,Rambr1,Linpe1,SMEG,AX774,Synfus1,Synplu1

=head1 AUTHOR

Jason Stajich @hyphaltip http//github.com/hyphaltip
Yan Wang @rockunderking https://github.com/YanWangTF

=cut

use strict;
use warnings;

use File::Spec;
use Getopt::Long;
# BioPerl modules
use Bio::SeqIO;
use Bio::DB::Fasta;

# hardcoded (for now)
my $pfamext = 'domtbl';
my $cazyext = 'tsv';
my $meropsext = 'blasttab';

# user set
my $evalue_cutoff = 1e-3;
my $input;
my $outdir;
my $db;
my $report;
my $pfamdir = 'Domains/Pfam';
my $cazydir = 'Domains/CAZY';
my $meropsdir = 'Domains/MEROPS';
my $ordering;
my $onlyshowrequested = 0;
GetOptions(
    'i|input:s'         => \$input,
    'o|outdir:s'        => \$outdir,
    'r|report:s'        => \$report,
    'p|pfam:s'          => \$pfamdir,
    'c|cazy:s'          => \$cazydir,
    'm|merops:s'          => \$meropsdir,
    'db:s'              => \$db,
    's|species|order:s' => \$ordering,
    'only|onlyreq!'     => \$onlyshowrequested,
    );

if ( ! $outdir ) {
    $outdir = $input . "-seqs";
    $outdir =~ s/\.out//;
}
my $report_fh;
if( $report ) {
    open($report_fh => ">$report") || die "Cannot open $report for writing: $!";
} else{
    $report_fh = \*STDOUT;
}
my $dbh = Bio::DB::Fasta->new($db);
my %gene2domains;
my %gene2domains_cazy;
my %gene2domains_merops;
my %pfam2desc;
if( $pfamdir && -d $pfamdir ) {
    opendir(my $dir => $pfamdir) || die "cannot open $pfamdir: $!";
    for my $file ( readdir($dir) ) {
	next unless $file =~ /(\S+)\.$pfamext$/;
	my $stem = $1;
	open(my $fh => File::Spec->catfile($pfamdir,$file) ) || 
	    die"cannot opend $dir/$file: $!";

	while(<$fh>) {
	    next if /^\#/;

	    my ($domain,$domaccno,$tlen,$gene_name,$qaccno,$qlen,
		$fullevalue,$fullscore,$fullbias,$n,$ntotal,$cvalue,$ivalue,
		$score,$dombias,
		$hstart,$hend, $qstart,$qend,$envfrom,$envto,$acc,$desc) =
		    split(/\s+/,$_,23);
	    my $evalue = $ivalue;
	    if( ! exists $pfam2desc{$domain} ) {
		$pfam2desc{$domain} = [$desc, $domaccno,$tlen ];
	    }
	    # not necessary if we are using pfam Gathering cutoff??
	    # next if $evalue > $evalue_cutoff;
	    $gene2domains{$gene_name}->{$domain}++;
	}
    }    
}

my %cazy2desc;
if( $cazydir && -d $cazydir ) {
    opendir(my $dir => $cazydir) || die "cannot open $cazydir: $!";
    for my $file ( readdir($dir) ) {
	next unless $file =~ /(\S+)\.$cazyext$/;
	my $stem = $1;
	open(my $fh => File::Spec->catfile($cazydir,$file) ) || 
	    die"cannot opend $dir/$file: $!";

	while(<$fh>) {
	    next if /^\#/;

	    my ($domain_cazy,$domlen,$gene_name,$genelen,
		$ivalue,$hstart,$hend,$qstart,$qend,$coverage) =
		    split(/\t/,$_,10);
	    my $evalue = $ivalue;
	    if( ! exists $cazy2desc{$domain_cazy} ) {
		$cazy2desc{$domain_cazy} = [$domlen, $genelen,$coverage ];
	    }
	    # not necessary if we are using pfam Gathering cutoff??
	    # next if $evalue > $evalue_cutoff;
	    $gene2domains_cazy{$gene_name}->{$domain_cazy}++;
	}
    }    
}

my %merops2desc;
if( $meropsdir && -d $meropsdir ) {
    opendir(my $dir => $meropsdir) || die "cannot open $meropsdir: $!";
    for my $file ( readdir($dir) ) {
	next unless $file =~ /(\S+)\.$meropsext$/;
	my $stem = $1;
	open(my $fh => File::Spec->catfile($meropsdir,$file) ) || 
	    die"cannot opend $dir/$file: $!";

	while(<$fh>) {
	    next if /^\#/;

	    my ($gene_name,$domain_merops,$identity,$alnlen,$nmismatch,$gap,$qstart,$qend,
		$hstart,$hend,$KA_evalue,$score) =
		    split(/\t/,$_,12);
	    my $evalue = $KA_evalue;
	    if( ! exists $merops2desc{$domain_merops} ) {
		$merops2desc{$domain_merops} = [$identity, $hstart,$hend ];
	    }
	    # not necessary if we are using pfam Gathering cutoff??
	    # next if $evalue > $evalue_cutoff;
	    $gene2domains_merops{$gene_name}->{$domain_merops}++;
	}
    }    
}

mkdir($outdir);
my @groups;
my %spall;
open(my $fh => $input) || die "cannot open $input: $!";

while(<$fh>) {
    my @genes = split;
    my $count = scalar @groups;
    my $countn = sprintf("%05d",$count);
    my $seqio = Bio::SeqIO->new(-format => 'fasta', 
				-file =>">$outdir/ORTHO_$countn.fa");
    my $sp = {};
    my $domains = {};
    my $domains_cazy = {};
    my $domains_merops = {};
    for my $gene ( @genes ) {
	
	my ($spname,$genename) = split(/\|/,$gene);
	if( exists $gene2domains{$gene} ) {
	    while (my ($d,$ct) = each %{$gene2domains{$gene}} ) {
		$domains->{$d} += $ct;
	    }
	}
	$spall{$spname}++;
	$sp->{$spname}++;
	if( exists $gene2domains_cazy{$gene} ) {
	    while (my ($d,$ct) = each %{$gene2domains_cazy{$gene}} ) {
		$domains_cazy->{$d} += $ct;
	    }
	}
	if( exists $gene2domains_merops{$gene} ) {
	    while (my ($d,$ct) = each %{$gene2domains_merops{$gene}} ) {
		$domains_merops->{$d} += $ct;
	    }
	}
	if( my $seq = $dbh->get_Seq_by_acc($gene) ) {
	    $seqio->write_seq($seq);
	} else {
	    warn("cannot find $gene in db $db");
	    next;
	}
    }
    push @groups, [ $sp, &domain_string($domains), &domain_string($domains_cazy), &domain_string($domains_merops)];
    #last if @groups > 10;
}

# lots of code to determine the species ordering of the columns... 
my @species;
if ( $ordering ) {
    my %obs;
    for my $n ( split(/,/,$ordering) ) {
	if ( ! exists $spall{$n} ) {
	    warn("requesting species $n but it does not exist in input\n");
	} else {
	    push @species, $n;
	    $obs{$n}++;
	}
    }
    unless( $onlyshowrequested ) { # if set the --only flag then
	                           # don't complain that order list is
	                           # incomplete
for my $e ( keys %spall ) {
	    unless ( $obs{$e} ) {
		warn("expected to see $e in the order list");
		push @species, $e;
	    }	    
	}
    }
} else {
    @species = sort { $spall{$b} <=> $spall{$a} } keys %spall;
}

print $report_fh join("\t", qw(ORTHOMCL), @species, qw(Pfam), qw(CAZY), qw(MEROPS)),"\n";
my $i = 0;
for my $group ( @groups ) {
    print $report_fh join("\t", sprintf("ORTHO_%05d",$i++),
			  (map { $group->[0]->{$_} || 0 } @species),
			  # this turns hash of icounts to a string
			  # it may need to be written to also show counts 
			  # of domains to show their continuity

$group->[1], $group->[2], $group->[3]),
    "\n";					
}

sub domain_string {
    my $domains = shift;
    join("; ",map { sprintf("%s (%d)", $_, $domains->{$_})} 
	 sort { $domains->{$b} <=> $domains->{$a} } keys %$domains);
}
