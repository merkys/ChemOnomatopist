#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    # From BBv2 P-65.1.2
    'CCCC(=O)O' => 'butanoic acid',
    'OC(=O)CCCCCCCCCCC(=O)O' => 'dodecanedioic acid',

    # From BBv2 P-65.1.2.2.1
    'OC(=O)CCC(C(=O)O)CCC(=O)O' => 'pentane-1,3,5-tricarboxylic acid',
    'C(C(=O)O)(C(=O)O)C(C(=O)O)(C(=O)O)' => 'ethane-1,1,2,2-tetracarboxylic acid',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
