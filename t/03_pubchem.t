#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

open (my $f, '<', 't/PubChemData') or die;

my %tests;
while (my $line = <$f>) {
  my @elems = split ' ', $line;
  $tests{$elems[2]} = $elems[1];
}

close $f;

plan tests => scalar keys %tests;

for my $case (keys %tests) {
    is( ChemOnomatopist::get_name( $case ), $tests{$case} );
}
