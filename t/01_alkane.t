#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use List::Util qw( all );
use Test::More;

my %SMILES_cases = (
    'CCCC'   => 'butane',
    'CCCCC'  => 'pentane',
    'CC(C)C' => '2-methylpropane', # FIXME: 'methylpropane'
    'C1CCC1' => 'cyclobutane',
    'c1ccccccccccccc1' => 'cyclotetradecaheptaene',
    'CC(C)CC(CCC(C)C)C' => '2,4,7-trimethyloctane',
    'C(C)C(CCC(CCC(C)C)(C)C)C' => '2,5,5,8-tetramethyldecane',
    'C(C)C(CC(C)C)CC' => '4-ethyl-2-methylhexane',
    'C(C)C(C(CC)(C)C)CCC' => '4-ethyl-3,3-dimethylheptane',
    'C(C)C(C(CCC)C)(C(CCCC)C)C' => '5-ethyl-4,5,6-trimethyldecane',
    'C(C)C(C(CC)C)C(C(CCC)(C)C)(CC)CC' => '4,5,5-triethyl-3,6,6-trimethylnonane',
    'C(C)C(C(C(CCC)C)(CCC)CCC)(CCCC)CCC' => '6-ethyl-4-methyl-5,5,6-tripropyldecane',
    'C(C)C(C(CCC)(C)C)(C(C(CCC)(C)CC)CCC)CCC' => '5,7-diethyl-4,4,7-trimethyl-5,6-dipropyldecane',
    'CCCCCCCCCCCCCCCCCCCCCCC' => 'tricosane',
    'CC(CC(CCC)CCC)C' => '4-(2-methylpropyl)heptane', # incorrect results for old method
    'CC(CC(CC(CC(CC)CC)C)(CC(CC(CC)C)C)CC(CC(CC)C)C)CC(CC)C' => '7,7-di(2,4-dimethylhexyl)-3-ethyl-5,9,11-trimethyltridecane', # incorrect results for both methods
    'CC(C(CCC)C)C(CC(CCCC)CC)CCCCCC' => '7-(1,2-dimethylpentyl)-5-ethyltridecane',
);

plan tests => 2 * scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case}, 'old';
    is ChemOnomatopist::get_name( $case, 1 ), $SMILES_cases{$case}, 'new';
}
