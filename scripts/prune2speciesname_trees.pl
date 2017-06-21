#!/usr/bin/perl

use strict;
use warnings;

use Bio::TreeIO;

my $in = Bio::TreeIO->new(-format => 'newick', -file => shift @ARGV);
my $out = Bio::TreeIO->new(-format => 'newick');
while( my $t = $in->next_tree ) {
 for my $node ( grep { $_->is_Leaf } $t->get_nodes() ) {
  my ($sp,$gene) = split(/[\|_]/,$node->id,2);
  $node->id($sp);
 }
 $out->write_tree($t);
 $out->_print("\n");
}
