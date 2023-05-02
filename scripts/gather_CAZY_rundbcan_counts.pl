#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);

use Getopt::Long;
use Pod::Usage;

my $man = 0;
my $help = 0;

my $sfetch = 'esl-sfetch';

my ($indir,$out,$domaindir) = qw(domains/CAZY summary domain_seq/CAZYdbcan);
my $outpref = 'CAZYdbcan';
my $ext = 'run_dbcan';
my $datfile = 'overview.txt';
my $seqdb;
my $speciesorder;
my $verbose;

GetOptions('help|?'               => \$help, 
	   'man'                  => \$man,
	   'i|input:s'            => \$indir,
	   'o|out:s'              => \$out,
	   'op|outpref:s'         => \$outpref,
	   'v|verbose!'           => \$verbose,
	   'd|domain|domainout:s' => \$domaindir,
	   'ext|extension:s'      => \$ext,
	   'fname:s'              => \$datfile,
	   'db|database|seqs:s'   => \$seqdb,
	   'species:s'            => \$speciesorder,
	   'sfetch:s'             => \$sfetch,
	   
    ) or pod2usage(2);
pod2usage(1) if $help;

if ( $ext !~ /^\./) {
  $ext = '.'. $ext;
}

mkdir($domaindir) unless -d $domaindir;

my (%table, %table_genes,%specieset,%counts);
opendir(my $ind => $indir) || die "cannot open $indir: $!";

for my $rundir ( readdir($ind) ) {
    next unless $rundir =~ /(\S+)\Q$ext\E$/;
    my $stem = $1;
    $stem =~ s/(sp|var)\./$1\_/; # get rid of sp/var
    my ($species) = split(/\./,$stem);
    warn("species=$species\n") if $verbose;

    my $filepath = File::Spec->catfile($indir,$rundir,$datfile);
    if ( ! -f $filepath ) {
	warn("no successful run to produce $datfile in $indir/$rundir");
	next;
    }
    open(my $in => $filepath ) || die "cannot open $filepath: $!";
#    my $header = <$in>;
    while(<$in>) {		# parse domtbl file
	next if /^\#/ || /^Gene ID/;		# skip comment lines
	chomp;
	my ($gene_name, $hmm, $hotpep, $diamond_name, $signalp, $toolcount) = split(/\t/,$_);
	if ( $toolcount == 3 ) {
	    $counts{$hotpep}->{$species}++;
	    if ( $hmm =~ /(\S+)\((\d+)\-(\d+)\)/ ) {
		my ($domain,$hstart,$hend) = ($1,$2);
		push @{$table{$hotpep}->{$species}}, [$gene_name,$hstart,$hend,
		$domain, $diamond_name,$signalp];
		$table_genes{$hotpep}->{$species}->{$gene_name}++;
		$specieset{$species}++;
	    }
	}
    }

    my @taxanames = sort keys %specieset;
    if ( $seqdb ) {
	my $sfetchexe = `which $sfetch`;
	chomp($sfetchexe);
	if ( ! -x $sfetchexe ) {
	    warn("cannot find esl-sfetch executable");
	    $seqdb = undef;
	}
	elsif( ! -f "$seqdb/allseq.ssi" ) {
	    `cat $seqdb/*.fasta > $seqdb/allseq`;
	    `$sfetch --index $seqdb/allseq`;
	}
    }

    mkdir($out) unless -d $out;
    open(my $fh => ">$out/$outpref\_counts.tsv") || die $!;
    open(my $fhgn => ">$out/$outpref\_counts_genes.tsv") || die $!;

    print $fh join("\t", qw(DOMAIN), @taxanames), "\n";
    print $fhgn join("\t", qw(DOMAIN), @taxanames), "\n";

    for my $s ( map { $_->[0] }
		sort { $b->[1] <=> $a->[1] }
		map { [$_, sum(values %{$counts{$_}})] } keys %table ) {

	print $fh join("\t",$s,map { scalar @{$table{$s}->{$_} || []} } @taxanames), "\n";
	print $fhgn join("\t",$s,
			 map { scalar keys %{$table_genes{$s}->{$_} || {}} }
			 @taxanames), "\n";
	if( $seqdb && -f "$seqdb/allseq.ssi" ) {
	    open(my $outseq => "| $sfetch -f $seqdb/allseq - > $domaindir/$s.fas") || die $!;
	    print $outseq join("\n", map { keys %{$table_genes{$s}->{$_} || {}} }
			       @taxanames ), "\n";
	    close($outseq);
	}
    }
}

=head1 NAME
    
gather_CAZY_domain_counts.pl - generate a table of domain counts from dbCAN results
    
=head1 USAGE
    
perl gather_CAZY_rundbcan_counts.pl -i domains -o table [-ext .run_dbcan]

=head1 DESCRIPTION

Generate table of domain counts for use in comparing, plotting as heatmap, etc from dbCAN CAZY run script (run_dbcan.py) which runs hotpep, hmmer, signalp, and diamond searches to confirm CAZY candidates

=head2 ARGUMENTS

 -i / --input    Input folder  [One file per taxon/proteome]
 -o / --output   Output table
 -db /--database folder or file of sequences
 --ext            Extension of the result folder. Default is .run_dbscan
 --fname         Name of run_dbcan outfile (default 'overiew.txt')

=head1 AUTHOR

Jason Stajich - jasonstajich.phd[at]gmail.com
University of California-Riverside

=cut
