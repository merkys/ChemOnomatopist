#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use List::Util qw( all );
use Test::More;

my %SMILES_cases = (
    'CO'          => 'methanol',
    'C(O)O'       => 'methanediol',
    'C(C)(C)(C)O' => '2-methylpropan-2-ol',
    'CC(O)CCO'    => 'butane-1,3-diol',
    'C(CO)C(CCO)CO' => '3-(hydroxymethyl)pentane-1,5-diol',

    # From BBv2 P-63.1.1.2
    'OCC(CO)(CO)CO' => '2,2-bis(hydroxymethyl)propane-1,3-diol',
    'CC(O)(C)C(O)(C)C' => '2,3-dimethylbutane-2,3-diol',
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
