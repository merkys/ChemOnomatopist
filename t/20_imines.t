#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    'C(C(=N)C(=O)O)C(=O)O' => '2-iminobutanedioic acid',
    'CC(C)C=N' => '2-methylpropan-1-imine',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
