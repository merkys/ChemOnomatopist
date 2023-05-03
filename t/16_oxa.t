#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    'COCCOCCOCCOCC' => '2,5,8,11-tetraoxatridecane', # BBv2 P-12.1

    # From BBv2 P-15.4.3.1
    'COCSSCCOCC[Se]C' => '2,8-dioxa-4,5-dithia-11-selenadodecane',
    '[Si]OCS[Si]' => '2-oxa-4-thia-1,5-disilapentane',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
