#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    'C1CCNC1' => 'pyrrolidine',

    # From BBv2 P-14.5.1
    'C1CCCCC1(C)CC' => '1-ethyl-1-methylcyclohexane',
    'CCC1CCC(C)CC1' => '1-ethyl-4-methylcyclohexane',
    # 'C(=O)(O)CC([Br])([Br])C1CCCCC1' => '3,3-dibromo-3-cyclohexylpropanoic acid', # FIXME: Does not work

    # From BBv2 P-14.5.3
    'C(C)(C)(C)C=1C=CC=C(C(C)CC)C=1' => '1-(butan-2-yl)-3-tert-butylbenzene', # FIXME: Does not work due to limitations on double bonds
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
