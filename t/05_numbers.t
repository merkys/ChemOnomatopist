#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

# Taken from:
# https://en.wikipedia.org/w/index.php?title=IUPAC_numerical_multiplier&oldid=1086173027
my %cases = (
     241 => 'hentetracontadict',
     411 => 'undecatetract',
     548 => 'octatetracontapentact',
    9267 => 'heptahexacontadictanonali',
);

plan tests => scalar( @ChemOnomatopist::prefixes ) - 5 +
              scalar keys %cases;

for (5..$#ChemOnomatopist::prefixes) {
    is( ChemOnomatopist::IUPAC_numerical_multiplier( $_ ),
        $ChemOnomatopist::prefixes[$_],
        "Number $_" );
}

for (sort keys %cases) {
    is( ChemOnomatopist::IUPAC_numerical_multiplier( $_ ),
        $cases{$_},
        "Number $_" );
}
