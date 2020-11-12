#!env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);

use Getopt::Long;
use Pod::Usage;

my $man = 0;
my $help = 0;
my $sfetch = 'esl-sfetch';

my ($indir,$out,$domaindir) = qw(search/cellwall summary domain_seq/cellwall );
my $outpref = 'domains';
my $ext = '.domtbl';
my $cutoff = 1e-5;
my $seqdb;
my $speciesorder;
my $verbose;

GetOptions('help|?'               => \$help, man => \$man,
	   'i|input:s'            => \$indir,
           'o|out:s'              => \$out,
	   'op|outpref:s'         => \$outpref,
	   'v|verbose!'           => \$verbose,
	   'd|domain|domainout:s' => \$domaindir,
	   'c|cutoff|evalue:s'    => \$cutoff,
	   'ext|extension:s'      => \$ext,
	   'db|database|seqs:s'   => \$seqdb,
	   'species:s'            => \$speciesorder,
	   'sftech:s'             => \$sfetch,

) or pod2usage(2);
pod2usage(1) if $help;

if ( $ext !~ /^\./) {
  $ext = '.'. $ext;
}

mkdir($domaindir) unless -d $domaindir;

my (%table, %table_genes,%specieset,%counts,%domain2acc);
opendir(my $ind => $indir) || die "cannot open $indir: $!";
for my $file ( readdir($ind) ) {
    next unless $file =~ /(\S+)\Q$ext\E$/;
    my $stem = $1;
    $stem =~ s/(sp|var)\./$1\_/; # get rid of sp/var 
    my ($species) = split(/\./,$stem);
    warn("species=$species\n") if $verbose;
    my $filepath = File::Spec->catfile($indir,$file);
    open(my $in => $filepath ) || die "cannot open $filepath: $!";
    while(<$in>) {		# parse domtbl file
	next if /^\#/;		# skip comment lines
	my ($gene_name,$qacc,$qlen,$domain,$domainacc,$domainlen,
	    $fullevalue,$fullscore,$fullbias,$n,$ntotal,$cvalue,$ivalue,
	    $score,$dombias,
	    $hstart,$hend, $qstart,$qend,$envfrom,$envto,$acc,$desc) =
		split(/\s+/,$_,23);
	$domain2acc{$domain} = $domainacc;
	my $evalue = $ivalue;
	if( $evalue > $cutoff ) {
	    # skip
	    next
	}
	$counts{$domain}->{$species}++;
	push @{$table{$domain}->{$species}}, [$gene_name,$hstart,$hend];
	$table_genes{$domain}->{$species}->{$gene_name}++;
	$specieset{$species}++;
    }
}

my @taxanames = sort keys %specieset;

if ( $seqdb ) {
    my $sfetchexe = `which $sfetch`;
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

print $fh join("\t", qw(DOMAIN ACCESSION), @taxanames), "\n";
print $fhgn join("\t", qw(DOMAIN ACCESSION), @taxanames), "\n";

for my $s ( map { $_->[0] }
	    sort { $b->[1] <=> $a->[1] }
	    map { [$_, sum(values %{$counts{$_}})] } keys %table ) {

    print $fh join("\t",$s,$domain2acc{$s},
		   map { scalar @{$table{$s}->{$_} || []} } @taxanames), "\n";
    print $fhgn join("\t",$s,$domain2acc{$s},
		     map { scalar keys %{$table_genes{$s}->{$_} || {}} } 
		     @taxanames), "\n";

    if( $seqdb && -f "$seqdb/allseq.ssi" ) {
	open(my $outseq => "| $sfetch -f $seqdb/allseq - > $domaindir/$s.fas") || die $!;
	print $outseq join("\n", map { keys %{$table_genes{$s}->{$_} || {}} }
			   @taxanames ), "\n";
	close($outseq);
    }
}


=head1 NAME

gather_hmmsearch_domain_counts.pl - generate a table of domain counts from hmmsearch of Pfam or other DB

=head1 USAGE

perl gather_hmmsearch_domain_counts.pl -i domains/cellwall -o table [-ext .domtbl]

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
