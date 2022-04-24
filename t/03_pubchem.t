#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

open (my $f, '<', 'PubChemData') or die;

my %tests;
while (my $line = <$f>) {
  my @elems = split ' ', $line;
  $tests{$elems[2]} = $elems[1];
}

close $f;

plan tests => scalar keys %tests;

for my $case (keys %tests) {
    # FIXME: Chemical name may have initial letter uppercased, but it may
    #        not be the right choice to lowercase it before comparison.
    #        Need to think a bit more on how to deal with it. (A.M.)
    is( ChemOnomatopist::get_name( $case ), $tests{$case} );
}
