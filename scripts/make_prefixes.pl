#!env perl
use strict;
my $dir= shift || "final_combine/pep";
opendir(DIR, $dir) || die "dir is $dir: $!";
print join ("\t", qw(Pref Name)),"\n";
for my $file ( sort readdir(DIR) ) {
    next unless $file =~ /(\S+)\.(pep|aa|TFASTX)\.fasta$/;
    my $stem = $1;
    $stem =~ s/\.\w+\.v\d+//;

    open(my $in => "head -n 1 $dir/$file |") || die $!;
    if( my $line = <$in> ) {
	if ( $line =~ />(\S+)/ ) {
	    my ($pref) = split(/\|/,$1);
	    print join("\t",$pref,$stem), "\n";
	}
    }
}
