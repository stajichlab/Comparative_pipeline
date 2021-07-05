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
    
my ($indir,$out,$domainseq) = qw(domains/MEROPS summary domain_seq/MEROPS);
my $outpref = 'MEROPS';
my $domain2family = qw(lib/merops_lib.families.tab);
my $ext = 'blasttab';
my $cutoff = 1e-10;
my $seqdb;
my $speciesorder;
my $verbose;

GetOptions('help|?'               => \$help, man => \$man,
	   'i|input:s'            => \$indir,
           'o|out:s'              => \$out,
	   'op|outpref:s'         => \$outpref,
	   'v|verbose!'           => \$verbose,
	   'd|domain|domainout:s' => \$domainseq,
	   'f|df|domain2family'   => \$domain2family,
	   'c|cutoff|evalue:s'    => \$cutoff,
	   'ext|extension:s'      => \$ext,
	   'db|database|seqs:s'   => \$seqdb,
	   'species:s'            => \$speciesorder,
	   'sfetch:s'             => \$sfetch,

) or pod2usage(2);
pod2usage(1) if $help;

if ( $ext !~ /^\./) {
  $ext = '.'. $ext;
}

mkdir($domainseq) unless -d $domainseq;

die ("must provide a valid domain2family file - is the lib directory symlinked to this folder?") 
    if ( ! -e $domain2family || -z $domain2family );

open(my $d2f => $domain2family ) || die "cannot open $domain2family: $!";
my (%domainsFamily,%domainsSubFamily);
while(<$d2f>) {
    my ($domain,$family,$subfamily) = split;
    $domainsFamily{$domain}    = $family;
    $domainsSubFamily{$domain} = $subfamily;
}

    
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
    my %seen;
    while(<$in>) {		# parse blasttab file
	next if /^\#/;		# skip comment lines
	chomp;
	my ($gene_name,$domain,$percentid,$ident,$mismatch,$gaps,$qstart,$qend,
	    $hstart,$hend,$evalue,$bits) = split(/\t/,$_);
	if( $evalue > $cutoff ) {
	    # skip
	    next;
	}
	next if $seen{$gene_name}++;
	my $family = $domainsFamily{$domain};
	if( ! $family ) {
	    warn("cannot find a family for $domain ($gene_name, $domain, $percentid)\n");
	    $family = $domain;
	}
	$counts{$family}->{$species}++;
	push @{$table{$family}->{$species}}, [$gene_name,$hstart,$hend];
	$table_genes{$family}->{$species}->{$gene_name}++;
	$specieset{$species}++;
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
	open(my $outseq => "| $sfetch -f $seqdb/allseq - > $domainseq/$s.fas") || die $!;
	print $outseq join("\n", map { keys %{$table_genes{$s}->{$_} || {}} }
			   @taxanames ), "\n";
	close($outseq);
    }
}

=head1 NAME

gather_MEROPS_family_counts.pl - generate a table of domain counts from dbCAN results

=head1 USAGE

perl gather_MEROPS_family_counts.pl -i domains -o table [-ext .blasttab]

=head1 DESCRIPTION

Generate table of domain counts for use in comparing, plotting as heatmap, etc

=head2 ARGUMENTS

-i / --input    Input folder  [One file per taxon/proteome]
-o / --output   Output table
-db /--database folder or file of sequences
-ext            Extension of the result files. Default is .blasttab
-c / --cutoff   evalue cutoff

=head1 AUTHOR

Jason Stajich - jasonstajich.phd[at]gmail.com
University of California-Riverside

=cut
