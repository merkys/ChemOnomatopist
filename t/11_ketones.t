#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use List::Util qw( all );
use Test::More;

my %SMILES_cases = (
    'CCCC=O'   => 'butanol',
    'CC(=O)CC' => 'butan-2-anol',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
