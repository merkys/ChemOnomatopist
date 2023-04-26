#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use List::Util qw( all );
use Test::More;

my %SMILES_cases = (
    'CO'          => 'methanol',
    'C(C)(C)(C)O' => '2-methylpropan-2-ol',
    'CC(O)CCO'    => 'butane-1,3-diol',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
