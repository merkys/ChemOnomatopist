#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    'CC(=O)CC(=O)CCCO' => '7-hydroxyheptane-2,4-dione',
    'CC(C)C(CCC(C)O)CCC(=O)C' => '8-hydroxy-5-propan-2-ylnonan-2-one',

    'C(O)OO' => '1-hydroperoxymethanol',
    'OOCC(=O)C' => '1-hydroperoxypropan-2-one', # Unchecked
    'OOCCO' => '2-hydroperoxyethan-1-ol', # BBv2 P-63.4.2.2

    # From BBv2 P-65.1.2.4
    'CC(=O)CCCC(=O)O' => '5-oxohexanoic acid',
    'C(=O)(O)C(O)C(C(=O)O)C(=O)C(=O)O' => '1-hydroxy-3-oxopropane-1,2,3-tricarboxylic acid',

    # 'O=CCCC(=O)O' => '4-oxobutanoic acid', # BBv2 P-66.6.1.3 - FIXME

    'CCCCCCCCC(C=O)C(CC)O' => '2-(1-hydroxypropyl)decanal',
    'S=C(CC(=O)O)C' => '3-sulfanylidenebutanoic acid',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
