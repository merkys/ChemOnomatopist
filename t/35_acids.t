#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'C(O)S(=O)O', iupac => 'hydroxymethanesulfinic acid' }, # PubChem 9000

    { smiles => 'S(=O)(O)CC(=O)O', iupac => 'sulfinoacetic acid' }, # From BBv2 P-65.3.2.1

    # From BBv2 P-65.3.1
    { smiles => 'C1(=CC=CC=C1)S(=O)(=O)O', iupac => 'benzenesulfonic acid' },
    { smiles => 'CC(CC)S(=O)O', iupac => 'butane-2-sulfinic acid' },
    { smiles => 'CC1=C(C=C(C=C1)S(=O)(=O)O)S(=O)(=O)O', iupac => '4-methylbenzene-1,3-disulfonic acid' },
    { smiles => 'NC1=CC=C(C=C1)S(=O)(=O)O', iupac => '4-aminobenzene-1-sulfonic acid' },
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
