#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    'COCCOCCOCCOCC' => '2,5,8,11-tetraoxatridecane', # BBv2 P-12.1
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
