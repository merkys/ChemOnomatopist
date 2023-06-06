#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-66.5.1.1.4
    { smiles => 'C(#N)C1=CC=C(O1)C(=O)O', iupac => '5-cyanofuran-2-carboxylic acid', AUTHOR => 1 },
    { smiles => 'C(#N)CCC(=O)O', iupac => '3-cyanopropanoic acid' },
    { smiles => 'C(#N)CC(CCC#N)CCC#N', iupac => '4-(cyanomethyl)heptanedinitrile', AUTHOR => 1 },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $name = $case->{iupac};

    is ChemOnomatopist::get_name( $case->{smiles} ), $name;
}
