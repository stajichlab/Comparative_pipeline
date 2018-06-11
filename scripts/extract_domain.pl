#!/usr/bin/env perl
use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;

# subroutines

sub extract_seq_domains {
    my ($seq,$locations) = @_;

    my @seqs;
    my $seq_id= $seq->display_id;
    my $i = 1;
    for my $loc ( @$locations ) {
	push @seqs, Bio::Seq->new(-id => sprintf("%s_%d_of_%d",
						 $seq_id,
						 $i++, scalar @$locations),
				  -seq => $seq->subseq(@$loc));
				  
    }
    return [@seqs];
}

sub process_Pfam {
    my ($folder,$ext, $cutoff, $domains, $dbh, $outdir) = @_;
    opendir(my $domaindir => $folder) || die $!;
    for my $file ( readdir($domaindir) ) {
	if ( $file =~ /(\S+)\.\Q$ext\E/) {
	    my $fh;
	    my $name = $1;
	    if( $file =~ /\.gz$/ ) {
		open($fh => sprintf("gzip -dc %s |",
				    File::Spec->catfile($folder,$file))) || die $!;
	    } else {
		open($fh => File::Spec->catfile($folder,$file)) || die $!;
	    }
	    my %seendomains;
	    while(<$fh>) {	# parse domtbl file
		next if /^\#/;	# skip comment lines
		my ($domain,$domainacc,$tlen,$gene_name,$qacc,$qlen,
		    $fullevalue,$fullscore,$fullbias,$n,$ntotal,$cvalue,$ivalue,
		    $score,$dombias,
		    $hstart,$hend, $qstart,$qend,$envfrom,$envto,$acc,$desc) =
			split(/\s+/,$_,23);

		next unless exists $domains->{$domain};
		my $evalue = $ivalue;
		if( $evalue > $cutoff ) {
		    # skip
		    next
		}
		push @{$seendomains{$domain}->{$gene_name}}, [$qstart, $qend];
	    }
	    for my $domain ( keys %seendomains ) {
		my $hits = $seendomains{$domain};
		my $outsubdir = File::Spec->catdir($outdir,
						   $domain);
		mkdir($outsubdir) unless -d $outsubdir;
		
		my $out_seqio = Bio::SeqIO->new(-format => 'fasta',
						-file   => ">$outsubdir/$name.domains.fasta");
		while (my ($seqname,$locs) = each %$hits ) { 
		    my $hitseq = $dbh->get_Seq_by_id($seqname);
		    if ( ! $hitseq ) {
			warn("cannot find ",$seqname," in seq database\n");
			next;
		    }
		    my $domain_chunks = &extract_seq_domains($hitseq, $locs);
		    $out_seqio->write_seq(@$domain_chunks);
		}
	    }
	}
    }
}



my @domains;
my $db = 'pep';
my $outdir = 'domain_aln';
my $domain_type = 'Pfam';
my $extension = 'domtbl';
my $cutoff = 1e-2;

GetOptions(
    'db:s'        => \$db,
    'i|domain:s'  => \@domains,
    'o|outdir:s'  => \$outdir,
    't|type:s'    => \$domain_type,
    );

mkdir($outdir) unless -d $outdir;

my $idx = Bio::DB::Fasta->new($db);

my %domains;
for my $d ( @domains ) {
    $domains{$d} = 1;
}
if ( $domain_type =~ /pfam/i ) {
    process_Pfam("domains/Pfam",
		  $extension,
		  $cutoff,
		  \%domains,
		  $idx,
		  $outdir);
}


1;
