#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    'CC(CN)O' => '1-aminopropan-2-ol',
    'C(CC(=O)N)C(=O)C(=O)O' => '5-amino-2,5-dioxopentanoic acid',
    'CC(CC(C(=O)O)N)N' => '2,4-diaminopentanoic acid',
    'C(CN)C=O' => '3-aminopropanal',
    'C(CCN)CCN' => 'pentane-1,5-diamine',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
