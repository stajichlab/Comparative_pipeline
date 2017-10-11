#!/usr/bin/perl -w
use strict;
use Getopt::Long;
#use Bio::DB::Fasta;
use Bio::SeqIO;

my $USAGE = <<EOF
orthomcl_extract_goodorth.pl
REQUIRED VALUES
  --db dbfile  location of fasta sequence file
  -i mcl_file     orthomcl result file
OPTIONAL OPTIONS
  --missing       number of taxa allowed to be missing in ortholog cluster (default: 2)
  --maxsp         number of gene members per species that are allowed (default: 1)

  -r report_file  write the ortholog report to a file, otherwise will be STDOUT
  -o  outputdir   where each ortholog sequence set should be written (default: orthomcl_orthologs)
  -p  prefix      prefix that ortholog sequence will be written with (default: orthomcl_)
  -s/--skip/skiptaxa skip these taxa in output
  -h/--help        display this help
EOF
    ;

my ($help,$db,$tab,$orthomcl);
my $allowed_missing  = 2;
my $outdir = 'orthomcl_good';
my $prefix = 'ORTHO_';
my $infile;
my $max_per_species = 1;
my $suffix = ".fa";
my $cdbyank = 'cdbyank';
my $cdbfasta = 'cdbfasta';
my @skip; # taxa to skip outputting
GetOptions(
    'i|in:s'         => \$infile,
    'db=s'           => \$db,
    'o|out|outdir:s' => \$outdir,
    'missing:i'      => \$allowed_missing,
    'maxsp:i'        => \$max_per_species,
    'p|prefix:s'     => \$prefix,
    'r|report:s'     => \$tab,
    'y|cdbyank|yank:s' => \$cdbyank,
    'idx|cdbfasta:s'  => \$cdbfasta,
    's|skip|skiptaxa|taxa:s' => \@skip,
	   'h|help'         => \$help,
	   );
die("usage: $USAGE\n") if( $help || ! $db);

mkdir($outdir) unless -d $outdir;
my $report_fh;

if( $tab ) {
    open($report_fh, ">$tab") || die("cannot open $tab: $!\n");
} else {
    $report_fh = \*STDOUT;
}
if ( ! -f "$db.cidx" ) {
 `cdbfasta $db`;
}
#my $dbh = Bio::DB::Fasta->new($dbdir);
my %skip_taxa = map { $_ => 1 } @skip;
my $in;
open($in, $infile) || die ("cannot open $infile, provide one via -i infile: $!\n");

my (@orthologs,%taxa);
my $groupn = 0;
while(<$in>) {
 my (@genes) = split;
 my $group = sprintf("ORTHO_%05d",$groupn);
 for my $gene ( @genes ) {
     my ($sp,$gname) = split(/\|/,$gene);
     next if ( $skip_taxa{$sp} ); # is okay to not use exists since we are looking for a non-zero (e.g was listed)
     $taxa{$sp} = 1;
     push @{$orthologs[$groupn]->{$sp}}, $gene;
 }
 $groupn++;
}

warn(scalar @orthologs, " orthologs captured\n");
my $i = 0;
my @taxa = sort keys %taxa;
my $ntaxa = scalar @taxa;
warn("ntaxa = $ntaxa\n");
print $report_fh join("\t", qw(ORTHOLOG_GRP NUM_MISSING ALL_SINGLE), @taxa), "\n";

for my $orth ( @orthologs ) {
    my @seent = keys %{$orth};
    if( ($ntaxa - scalar @seent) <= $allowed_missing ) {
	my @genes;
	my $all_single = 'yes';
	my $num_missing = ($ntaxa - scalar @seent);
	for my $sp ( sort keys %$orth ) {
	    my $g = $orth->{$sp};
	    $all_single = 'no' if ( scalar @{$orth->{$sp}} != 1 );
	    if( scalar @$g > $max_per_species ) {
		@genes = ();
		last;
	    } else {
		push @genes, @$g;
	    }
	}
	if( @genes ) {
	    my $outfile = File::Spec->catfile($outdir, $prefix . $i . $suffix);
	    open(my $outseq => "| cdbyank $db.cidx > $outfile") || die $!;
	    # @genes = sort @genes;
	    #my $outseq = Bio::SeqIO->new(-format => 'fasta',
	    #				 -file   => ">$outfile");
	    print $outseq join("\n", @genes), "\n";
#	    for my $gene ( @genes) {
#		print 
#		my $geneseq = $dbh->get_Seq_by_acc($gene);
#		if( ! defined $geneseq ) {
#		    warn("cannot find $gene in db $db\n");
#		    next;
#		}
#		$outseq->write_seq($geneseq);
#	    }
	    print $report_fh join("\t", sprintf("%s%d",$prefix,$i), 
				  $num_missing, $all_single, 
				  map {scalar @{$orth->{$_} || []} } @taxa), "\n";
	}
    } else {
#	warn("saw ", scalar @seent, " taxa in orth $i\n");
    }
    $i++;
}
