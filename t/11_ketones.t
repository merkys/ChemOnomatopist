#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    # 'CCCC=O'   => 'butan-1-one', # FIXME
    'CC(=O)CC' => 'butan-2-one',
    'CCCC(=O)CCCC(=O)CCCC(=O)CCCC(=O)C' => 'heptadecane-2,6,10,14-tetrone',
    'CC(C)C(=O)CC(=O)C' => '5-methylhexane-2,4-dione',
    'CC(C)CC(=O)C(CCC(=O)C)C(=O)CC(C)C' => '8-methyl-5-(3-methyl-1-oxobutyl)nonane-2,6-dione', # PubChem has '8-methyl-5-(3-methylbutanoyl)nonane-2,6-dione'

    # From BBv2 P-64.6.1
    'CC(=S)CC' => 'butane-2-thione',
    'CCC(CCC)=[Se]' => 'hexane-3-selone',
    'CC(CC(C)=S)=S' => 'pentane-2,4-dithione',

    'OC(C)C(C(C(C(C)O)=C)=O)=C' => '2,6-dihydroxy-3,5-dimethylideneheptan-4-one', # BBv2 P-64.7.1
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
