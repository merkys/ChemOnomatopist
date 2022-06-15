#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

plan tests => scalar @ChemOnomatopist::prefixes - 5;

for (5..$#ChemOnomatopist::prefixes) {
    is( ChemOnomatopist::IUPAC_numerical_multiplier( $_ ),
        $ChemOnomatopist::prefixes[$_],
        "Number $_" );
}
