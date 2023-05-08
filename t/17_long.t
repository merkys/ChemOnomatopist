#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

if( !$ENV{EXTENDED_TESTING} ) {
    plan skip_all => "Skip \$ENV{EXTENDED_TESTING} is not set\n";
}

my %SMILES_cases = (
    'CCCCCC(CC(CC)CC)CC(CC(CC(CC)CC)CC(CC)CC)CC(CC(CC(CC)CC)CC(CC)CC)CC(CC(CC(CC(CC(CC)CC)CC(CC)CC)CC(CC(CC)CC)CC(CC)CC)CC(CC(CC(CC)CC)CC(CC)CC)CC(CC(CC)CC)CC(CC)CC)CC(CC(CC(CC(CC)CC)CC(CC)CC)CC(CC(CC)CC)CC(CC)CC)CC(CC(CC(CC)CC)CC(CC)CC)CC(CC(CC)CC)CC(CC)CC' => '3-ethyl-5,17-bis(2-ethylbutyl)-11-[8-ethyl-6-(2-ethylbutyl)-2-[6-ethyl-4-(2-ethylbutyl)-2-[4-ethyl-2-(2-ethylbutyl)hexyl]octyl]-4-[4-ethyl-2-(2-ethylbutyl)hexyl]decyl]-9-[6-ethyl-4-(2-ethylbutyl)-2-[4-ethyl-2-(2-ethylbutyl)hexyl]octyl]-7,13,15-tris[4-ethyl-2-(2-ethylbutyl)hexyl]docosane',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    my $name = $SMILES_cases{$case};

    is ChemOnomatopist::get_name( $case ), $name;
}
