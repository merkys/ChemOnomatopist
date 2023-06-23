#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'CCCCCCCC(=O)OC(C)(C)C', iupac => 'tert-butyl octanoate', AUTHOR => 1 }, # BBv2 P-65.6.3.3.1
    { smiles => 'CCOC(=O)CC(=O)OC', iupac => 'ethyl methyl propanedioate', AUTHOR => 1 }, # BBv2 P-65.6.3.3.2.1
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
