#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-65.1.2
    { smiles => 'CCCC(=O)O', iupac => 'butanoic acid' },
    { smiles => 'OC(=O)CCCCCCCCCCC(=O)O', iupac => 'dodecanedioic acid' },

    # From BBv2 P-65.1.2.2.1
    { smiles => 'OC(=O)CCC(C(=O)O)CCC(=O)O', iupac => 'pentane-1,3,5-tricarboxylic acid' },
    { smiles => 'C(C(=O)O)(C(=O)O)C(C(=O)O)(C(=O)O)', iupac => 'ethane-1,1,2,2-tetracarboxylic acid' },

    { smiles => 'C(=O)O', iupac => 'formic acid' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $ok;
    eval { $ok = is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
    diag 'test supposed to fail with AUTHOR_TESTING' if $case->{AUTHOR} && $ok;
}
