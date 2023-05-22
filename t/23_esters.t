#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'CCCCCCCC(=O)OC(C)(C)C', iupac => 'tert-butyl octanoate' },
    { smiles => 'CCOC(=O)CC(=O)OC', iupac => 'ethyl methyl propanedioate', AUTHOR => 1 },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac};
}
