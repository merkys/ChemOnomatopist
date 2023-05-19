#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Old;
use Test::More;

my @cases = (
    { smiles => 'CCCC',   iupac => 'butane' },
    { smiles => 'CCCCC',  iupac => 'pentane' },
    { smiles => 'CC(C)C', iupac => '2-methylpropane' },

    { smiles => 'C1CCC1', iupac => 'cyclobutane' },
    { smiles => 'c1ccccccccccccc1', iupac => 'cyclotetradecaheptaene' },

    { smiles => 'CCCCCCCCCCC(C)CCCC',  iupac => '5-methylpentadecane' },
    { smiles => 'CC(C)CC(CCC(C)C)C',   iupac => '2,4,7-trimethyloctane' },
    { smiles => 'CCC(C)(C)CCCCCC(C)C', iupac => '2,8,8-trimethyldecane' },
    { smiles => 'C(C)C(CCC(CCC(C)C)(C)C)C', iupac => '2,5,5,8-tetramethyldecane' },
    { smiles => 'C(C)C(CC(C)C)CC', iupac => '4-ethyl-2-methylhexane' },
    { smiles => 'C(C)C(C(CC)(C)C)CCC', iupac => '4-ethyl-3,3-dimethylheptane' },

    { smiles => 'CCCC(CCC)(C(C)C)C(C)C', iupac => '4,4-di(propan-2-yl)heptane' },
    { smiles => 'CCCC(CCC)C(C)C', iupac => '4-propan-2-ylheptane' },
    { smiles => 'CCCCC(CCCC)C(C)C(C)C', iupac => '5-(3-methylbutan-2-yl)nonane' },
    { smiles => 'CCCCCC(CC(C)CC)CC(CCCCC)C(C)CCC', iupac => '6-(2-methylbutyl)-8-pentan-2-yltridecane' },
    { smiles => 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(CCCCCCCCCC)CCCCCCCCCCCCC', iupac => '14-decyltritetracontane' },
    { smiles => 'CCCCCCCC(CCCCCCC)(CCCCCCC)C(C)CC', iupac => '8-butan-2-yl-8-heptylpentadecane' },
    { smiles => 'CCC(C)C(CC(C(C)CC)C(C)CC)C(C)CC', iupac => '4,6-di(butan-2-yl)-3,7-dimethylnonane' },
    { smiles => 'CCCCCCCCCCCCC(CCC)(CCC)C(C(C)(C)C)C(CCCC)(CCCCCC)C(CC)(CCC)C(C)CCCC', iupac => '7-butyl-8-tert-butyl-6-ethyl-7-hexyl-5-methyl-6,9,9-tripropylhenicosane', AUTHOR => 1 },

    { smiles => 'C(C)C(C(CCC)C)(C(CCCC)C)C', iupac => '5-ethyl-4,5,6-trimethyldecane' },
    { smiles => 'C(C)C(C(CC)C)C(C(CCC)(C)C)(CC)CC', iupac => '4,5,5-triethyl-3,6,6-trimethylnonane' },
    { smiles => 'C(C)C(C(C(CCC)C)(CCC)CCC)(CCCC)CCC', iupac => '6-ethyl-4-methyl-5,5,6-tripropyldecane' },
    { smiles => 'C(C)C(C(CCC)(C)C)(C(C(CCC)(C)CC)CCC)CCC', iupac => '5,7-diethyl-4,4,7-trimethyl-5,6-dipropyldecane' },
    { smiles => 'CCC(CC)C(C)C', iupac => '3-ethyl-2-methylpentane' },
    { smiles => 'CCCCCCCCCCCCCCCCCCCCCCC', iupac => 'tricosane' },
    { smiles => 'CC(CC(CCC)CCC)C', iupac => '2-methyl-4-propylheptane' },
    { smiles => 'CC(CC(CC(CC(CC)CC)C)(CC(CC(CC)C)C)CC(CC(CC)C)C)CC(CC)C', iupac => '7,7-bis(2,4-dimethylhexyl)-3-ethyl-5,9,11-trimethyltridecane' },
    { smiles => 'CC(C(CCC)C)C(CC(CCCC)CC)CCCCCC', iupac => '5-ethyl-7-(3-methylhexan-2-yl)tridecane' },
    { smiles => 'CCC(CC)CCC(CCC(CC)CC)CCC(CCC(CCC(CC)CC)CCC(CC)CC)CCC(CCC(CC)CC)CCC(CC)CC', iupac => '3,15-diethyl-9-[6-ethyl-3-(3-ethylpentyl)octyl]-6,12-bis(3-ethylpentyl)heptadecane', AUTHOR => 1 }, # different order
    { smiles => 'CC(C)CC(CC(C)C)CC(CC(C)C)CC(C)C', iupac => '2,8-dimethyl-4,6-bis(2-methylpropyl)nonane' },

    { smiles => 'CCCCCCCCCC(CCCC)(CCCC)C(C)(C)C', iupac => '5-tert-butyl-5-butyltetradecane' }, # PubChem has 5-butyl-5-tert-butyltetradecane (?)
    { smiles => 'CCCCC(CC)C(CCCC)C(C)(C)C', iupac => '5-tert-butyl-6-ethyldecane' },
    { smiles => 'CCCCC(CCCC)CCCCCCCC(C(C)CCC)C(C)(C)C', iupac => '5-tert-butyl-13-butyl-4-methylheptadecane' },

    { smiles => 'CCCCCC(C)(C)C(C)CC(CC)CC(CCC(CC)CC)C(CCC(CC)CCC)(CC(CC)CCCC)C(CC)(CCCC(CC)CC)C(CC(C)C(CC)CC)(CC(C)(CC)CCC)C(C)C(CC(C)CC(C)CC)(CC(C)C(CC)CCC)CC(C)(CC)CCCC', iupac => '10,14-diethyl-11-(2-ethylhexyl)-11-(3-ethylhexyl)-10-(4-ethylhexyl)-7-(2-ethyl-2-methylhexyl)-7-(3-ethyl-2-methylhexyl)-9-(2-ethyl-2-methylpentyl)-9-(3-ethyl-2-methylpentyl)-12-(3-ethylpentyl)-3,5,8,16,17,17-hexamethyldocosane' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => 2 * @cases;

my %old_method_mistakes = (
    'CCCCC(CC)C(CCCC)C(C)(C)C' => '6-tert-butyl-5-ethyldecane',
);

for my $case (@cases) {
    my $name = $case->{iupac};

    is ChemOnomatopist::get_name( $case->{smiles} ), $name, 'new';

    if( $old_method_mistakes{$case->{smiles}} ) {
        $name = $old_method_mistakes{$case->{smiles}};
    }

    is ChemOnomatopist::Old::get_name( $case->{smiles} ), $name, 'old';
}
