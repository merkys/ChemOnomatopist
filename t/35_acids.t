#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'C(O)S(=O)O', iupac => 'hydroxymethanesulfinic acid' }, # PubChem 9000

    { smiles => 'S(=O)(O)CC(=O)O', iupac => 'sulfinoacetic acid' }, # From BBv2 P-65.3.2.1
    { smiles => 'CC(CC)S(=O)O', iupac => 'butane-2-sulfinic acid' }, # From BBv2 P-65.3.1
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
