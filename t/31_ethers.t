#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-63.2.4.1
    { smiles => 'COC',  iupac => 'methoxymethane', AUTHOR => 1 },
    { smiles => 'CCOC', iupac => 'methoxyethane', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)OC', iupac => 'anisole' },
    { smiles => 'COC1=CC2=CC=CC=C2C=C1', iupac => '2-methoxynaphthalene' },
    { smiles => 'ClCCOCC', iupac => '1-chloro-2-ethoxyethane', AUTHOR => 1 },
    { smiles => 'COCCOC', iupac => '1,2-dimethoxyethane', AUTHOR => 1 },
    { smiles => 'COCCOCCOC', iupac => '1-methoxy-2-(2-methoxyethoxy)ethane', AUTHOR => 1 },
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
