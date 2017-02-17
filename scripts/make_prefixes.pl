#!env perl
use strict;

# This script will generate a table of prefix to long name for datasets
# by parsing first sequence in fasta file and looking at header
# looking for ">prefix|genename" and the table will have the form
# prefix	FILENAME
# attempts are make to remove version nummber encoding in the filename strings
my $dir= shift || "pep";
opendir(DIR, $dir) || die "dir is $dir: $!";
print join ("\t", qw(Pref Name)),"\n";
for my $file ( sort readdir(DIR) ) {
    next unless $file =~ /(\S+)\.(pep|aa|TFASTX)\.fasta$/;
    my $stem = $1;
    $stem =~ s/\.\w+\.v\d+//;
    $stem =~ s/\.v\d+//;

    open(my $in => "head -n 1 $dir/$file |") || die $!;
    if( my $line = <$in> ) {
	if ( $line =~ />(\S+)/ ) {
	    my ($pref) = split(/\|/,$1);
	    print join("\t",$pref,$stem), "\n";
	}
    }
}
