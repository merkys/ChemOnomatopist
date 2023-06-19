#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'CC(CN)O', iupac => '1-aminopropan-2-ol' },
    { smiles => 'C(CC(=O)N)C(=O)C(=O)O', iupac => '5-amino-2,5-dioxopentanoic acid' },
    { smiles => 'CC(CC(C(=O)O)N)N', iupac => '2,4-diaminopentanoic acid' },
    { smiles => 'C(CN)C=O', iupac => '3-aminopropanal' },
    { smiles => 'C(CCN)CCN', iupac => 'pentane-1,5-diamine' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    eval { is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
}
