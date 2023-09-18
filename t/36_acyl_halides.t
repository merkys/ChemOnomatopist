#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-65.5.1
    { smiles => 'C(CCCCC)(=O)F', iupac => 'hexanoyl fluoride' },
    { smiles => 'C(CC(=O)Cl)(=O)Cl', iupac => 'propanedioyl dichloride', AUTHOR => 1 },
    { smiles => 'C1(=CC=C(C=C1)C(=O)Cl)C(=O)Cl', iupac => 'benzene-1,4-dicarbonyl dichloride', AUTHOR => 1 },
    { smiles => 'C(CCC(=O)Cl)(=O)Br', iupac => 'butanedioyl bromide chloride', AUTHOR => 1 },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $ok;
    eval { $ok = is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac}, $case->{smiles} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
    diag 'test supposed to fail with AUTHOR_TESTING' if $case->{AUTHOR} && $ok;
}
