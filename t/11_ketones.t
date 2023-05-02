#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    'CCCC=O'   => 'butan-1-one',
    'CC(=O)CC' => 'butan-2-one',
    'CCCC(=O)CCCC(=O)CCCC(=O)CCCC(=O)C' => 'heptadecane-2,6,10,14-tetrone',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
