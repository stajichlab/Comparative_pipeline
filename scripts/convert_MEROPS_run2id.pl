#!env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);

use Getopt::Long;
use Pod::Usage;


my $man = 0;
my $help = 0;
    
my ($infile);
my $outpref = 'MEROPS';
my $outext  = 'tsv';
my $domain2family = qw(lib/merops_lib.families.tab);
my $ext = 'blasttab';
my $cutoff = 1e-10;
my $seqdb;
my $speciesorder;
my $verbose;

GetOptions('help|?'               => \$help, man => \$man,
	   'i|input:s'            => \$infile,
	   'op|outpref:s'         => \$outpref,
	   'v|verbose!'           => \$verbose,
	   'f|df|domain2family'   => \$domain2family,
	   'c|cutoff|evalue:s'    => \$cutoff,
	   'ext|extension:s'      => \$ext,
	   'species:s'            => \$speciesorder,

) or pod2usage(2);
pod2usage(1) if $help;

if ( $ext !~ /^\./) {
  $ext = '.'. $ext;
}

die ("must provide a valid domain2family file - is the lib directory symlinked to this folder?") 
    if ( ! -e $domain2family || -z $domain2family );

open(my $d2f => $domain2family ) || die "cannot open $domain2family: $!";
my (%domainsFamily,%domainsSubFamily);
while(<$d2f>) {
    my ($domain,$family,$subfamily) = split;
    $domainsFamily{$domain} = $family;
    $domainsSubFamily{$domain} = $subfamily;
}

    
my (%table, %table_genes,%specieset,%counts);
my @input_files;
if ( ! -d $infile ) {
    push @input_files, $infile;
} else {
    my $indir = $infile;
    opendir(my $ind => $indir) || die "cannot open $indir: $!";
    for my $file ( readdir($ind) ) {
	next unless $file =~ /(\S+)\Q$ext\E$/;
	my $filepath = File::Spec->catfile($indir,$file);
	push @input_files, $filepath;
    }
}
for my $filepath ( @input_files ) {
    open(my $in => $filepath ) || die "cannot open $filepath: $!";
    my %seen;
    my $outname = $filepath;
    $outname =~ s/\Q$ext\E//;
    $outname = sprintf("%s.%s.%s",$outname,$outpref,$outext);
    open(my $outf => ">$outname") || die "cannot open $outname: $!";
    print $outf join("\t", qw(PROTEIN MEROPS_RUN MEROPS_ID PERCENTID QSTART QEND HSTART HEND EVALUE BITS)),"\n";
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
	    warn("cannot find a family for $domain");
	    $family = $domain;
	}
	print $outf join("\t", ($gene_name, $domain, $family,
				$qstart,$qend,$hstart,$hend, $evalue, 
				$bits)),"\n";
    }
}


=head1 NAME

convert_MEROPS_run2id.pl - generate a MEROPS table from MEROPS blast results

=head1 USAGE

perl convert_MEROPS_run2id.pl -i domains/MEROPS [-ext .blasttab]

=head1 DESCRIPTION

Generate table of domain counts for use in comparing, plotting as heatmap, etc

=head2 ARGUMENTS

-i / --input    Input folder  [One file per taxon/proteome]
-ext            Extension of the result files. Default is .blasttab
-c / --cutoff   evalue cutoff

=head1 AUTHOR

Jason Stajich - jasonstajich.phd[at]gmail.com
University of California-Riverside

=cut
