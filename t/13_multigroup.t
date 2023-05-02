#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    'CC(=O)CC(=O)CCCO' => '7-hydroxyheptane-2,4-dione',
    'CC(C)C(CCC(C)O)CCC(=O)C' => '8-hydroxy-5-propan-2-ylnonan-2-one',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
