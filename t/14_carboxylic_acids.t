#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    # From BBv2 P-65.1.2
    'CCCC(=O)O' => 'butanoic acid',
    'OC(=O)CCCCCCCCCCC(=O)O' => 'dodecanedioic acid',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
