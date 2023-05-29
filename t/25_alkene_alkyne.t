#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-44.4.1.1
    { smiles => 'CC=CC#C', iupac => 'pent-3-en-1-yne' },
    { smiles => 'CCC=CC', iupac => 'pent-2-ene' },
    { smiles => 'CCOCCOCCOCCOC=C', iupac => '3,6,9,12-tetraoxatetradec-1-ene' },

    # From BBv2 P-44.4.1.2
    { smiles => 'C1C=CCCCCCCCCCCCCCCCCC1', iupac => 'cycloicosene', AUTHOR => 1 },
    { smiles => 'C1C#CCCCCCCCCCCCCCCCCC1', iupac => 'cycloicosyne', AUTHOR => 1 },
    { smiles => 'C1C=CCCCCC=CCCCCCCCCCCC1', iupac => 'cycloicosa-1,8-diene', AUTHOR => 1 },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $name = $case->{iupac};

    is ChemOnomatopist::get_name( $case->{smiles} ), $name;
}
