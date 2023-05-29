#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-15.1.7.2.3
    { smiles => 'C=CCCO', iupac => 'but-3-en-1-ol' },
    { smiles => 'C=CC(C)C', iupac => '3-methylbut-1-ene' },
    { smiles => 'CC=CC(O)CC(C)C', iupac => '6-methylhept-2-en-4-ol' },

    { smiles => 'COCCOCCOCCOCC=C', iupac => '2,5,8,11-tetraoxatetradec-13-ene' }, # BBv2 P-31.1.2.2.2

    # From BBv2 P-44.4.1.1
    { smiles => 'CC=CC#C', iupac => 'pent-3-en-1-yne' },
    { smiles => 'CCC=CC', iupac => 'pent-2-ene' },
    { smiles => 'CCOCCOCCOCCOC=C', iupac => '3,6,9,12-tetraoxatetradec-1-ene' },

    # From BBv2 P-44.4.1.2
    { smiles => 'C1C=CCCCCCCCCCCCCCCCCC1', iupac => 'cycloicosene' },
    { smiles => 'C1C#CCCCCCCCCCCCCCCCCC1', iupac => 'cycloicosyne' },
    { smiles => 'C1=CCCCCCC=CCCCCCCCCCCC1', iupac => 'cycloicosa-1,8-diene', AUTHOR => 1 },

    { smiles => 'CCC(=C)CCC', iupac => '3-methylidenehexane' }, # BBv2 P-61.2.1
    { smiles => 'C(C)=C(C(C)C1=NC=CC=C1)CCC=C(C)C', iupac => '2-(3-ethylidene-7-methyloct-6-en-2-yl)pyridine' }, # BBv2 P-61.2.4
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $name = $case->{iupac};

    is ChemOnomatopist::get_name( $case->{smiles} ), $name;
}
