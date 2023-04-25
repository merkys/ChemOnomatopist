#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Old;
use List::Util qw( all );
use Test::More;

my %SMILES_cases = (
    'CCCC'   => 'butane',
    'CCCCC'  => 'pentane',
    'CC(C)C' => 'methylpropane',

    'C1CCC1' => 'cyclobutane',
    'c1ccccccccccccc1' => 'cyclotetradecaheptaene',

    'CCCCCCCCCCC(C)CCCC' => '5-methylpentadecane',
    'CC(C)CC(CCC(C)C)C' => '2,4,7-trimethyloctane',
    'CCC(C)(C)CCCCCC(C)C' => '2,8,8-trimethyldecane',
    'C(C)C(CCC(CCC(C)C)(C)C)C' => '2,5,5,8-tetramethyldecane',
    'C(C)C(CC(C)C)CC' => '4-ethyl-2-methylhexane',
    'C(C)C(C(CC)(C)C)CCC' => '4-ethyl-3,3-dimethylheptane',

    'CCCC(CCC)(C(C)C)C(C)C' => '4,4-di(propan-2-yl)heptane',
    'CCCC(CCC)C(C)C' => '4-propan-2-ylheptane',
    'CCCCC(CCCC)C(C)C(C)C' => '5-(3-methylbutan-2-yl)nonane',
    'CCCCCC(CC(C)CC)CC(CCCCC)C(C)CCC' => '6-(2-methylbutyl)-8-pentan-2-yltridecane',

    'C(C)C(C(CCC)C)(C(CCCC)C)C' => '5-ethyl-4,5,6-trimethyldecane',
    'C(C)C(C(CC)C)C(C(CCC)(C)C)(CC)CC' => '4,5,5-triethyl-3,6,6-trimethylnonane',
    'C(C)C(C(C(CCC)C)(CCC)CCC)(CCCC)CCC' => '6-ethyl-4-methyl-5,5,6-tripropyldecane',
    'C(C)C(C(CCC)(C)C)(C(C(CCC)(C)CC)CCC)CCC' => '5,7-diethyl-4,4,7-trimethyl-5,6-dipropyldecane',
    'CCC(CC)C(C)C' => '3-ethyl-2-methylpentane',
    'CCCCCCCCCCCCCCCCCCCCCCC' => 'tricosane',
    'CC(CC(CCC)CCC)C' => '2-methyl-4-propylheptane',
    'CC(CC(CC(CC(CC)CC)C)(CC(CC(CC)C)C)CC(CC(CC)C)C)CC(CC)C' => '7,7-bis(2,4-dimethylhexyl)-3-ethyl-5,9,11-trimethyltridecane',
    'CC(C(CCC)C)C(CC(CCCC)CC)CCCCCC' => '5-ethyl-7-(3-methylhexan-2-yl)tridecane',
    # 'CCC(CC)CCC(CCC(CC)CC)CCC(CCC(CCC(CC)CC)CCC(CC)CC)CCC(CCC(CC)CC)CCC(CC)CC' => '3,15-diethyl-9-[6-ethyl-3-(3-ethylpentyl)octyl]-6,12-bis(3-ethylpentyl)heptadecane', # different order
    'CC(C)CC(CC(C)C)CC(CC(C)C)CC(C)C' => '2,8-dimethyl-4,6-bis(2-methylpropyl)nonane',

    'CCCCCCCCCC(CCCC)(CCCC)C(C)(C)C' => '5-tert-butyl-5-butyltetradecane', # PubChem has 5-butyl-5-tert-butyltetradecane (?)
    'CCCCC(CC)C(CCCC)C(C)(C)C' => '5-tert-butyl-6-ethyldecane', # fails for old method
    'CCCCC(CCCC)CCCCCCCC(C(C)CCC)C(C)(C)C' => '5-tert-butyl-13-butyl-4-methylheptadecane',

    'CCCCCC(C)(C)C(C)CC(CC)CC(CCC(CC)CC)C(CCC(CC)CCC)(CC(CC)CCCC)C(CC)(CCCC(CC)CC)C(CC(C)C(CC)CC)(CC(C)(CC)CCC)C(C)C(CC(C)CC(C)CC)(CC(C)C(CC)CCC)CC(C)(CC)CCCC' => '10,14-diethyl-11-(2-ethylhexyl)-11-(3-ethylhexyl)-10-(4-ethylhexyl)-7-(2-ethyl-2-methylhexyl)-7-(3-ethyl-2-methylhexyl)-9-(2-ethyl-2-methylpentyl)-9-(3-ethyl-2-methylpentyl)-12-(3-ethylpentyl)-3,5,8,16,17,17-hexamethyldocosane',
);

plan tests => 2 * scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case}, 'new';
    is ChemOnomatopist::Old::get_name( $case ), $SMILES_cases{$case}, 'old';
}
