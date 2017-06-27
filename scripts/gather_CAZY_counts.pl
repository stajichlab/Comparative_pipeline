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
my ($indir,$out,$domaindir) = qw(domains/CAZY summary domain_seq);
my $outpref = 'CAZY';
my $ext = 'tsv';
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
) or pod2usage(2);
pod2usage(1) if $help;

if ( $ext !~ /^\./) {
  $ext = '.'. $ext;
}

mkdir($domaindir) unless -d $domaindir;

my (%table, %table_genes,%specieset,%counts);
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
	chomp;
	my ($domain,$domainlen,$gene_name,$genelen,$evalue,
	    $hstart,$hend,$qstart,$qend,$coverage) = split(/\t/,$_);

	if( $evalue > $cutoff ) {
	    # skip
	    next;
	}
	$counts{$domain}->{$species}++;
	push @{$table{$domain}->{$species}}, [$gene_name,$hstart,$hend];
	$table_genes{$domain}->{$species}->{$gene_name}++;
	$specieset{$species}++;
    }
}

my @taxanames = sort keys %specieset;

my $db;
#if ( $seqdb ) {
# $db = Bio::DB::Fasta->new($seqdb);
#}
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
}


=head1 NAME

gather_CAZY_domain_counts.pl - generate a table of domain counts from dbCAN results

=head1 USAGE

perl gather_CAZY_domain_counts.pl -i domains -o table [-ext .domtbl]

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
