#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    'CC(=O)CC(=O)CCCO' => '7-hydroxyheptane-2,4-dione',
    'CC(C)C(CCC(C)O)CCC(=O)C' => '8-hydroxy-5-propan-2-ylnonan-2-one',

    'OOCC(=O)C' => '1-hydroperoxypropan-2-one', # Unchecked
    'OOCCO' => '2-hydroperoxyethan-1-ol', # BBv2 P-63.4.2.2
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
