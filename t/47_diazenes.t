#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-68.3.1.3.1
    { smiles => 'C1(=CC=CC=C1)N=NC#N', iupac => 'phenyldiazenecarbonitrile' },
    { smiles => 'N(=N)CCC(=O)O', iupac => '3-diazenylpropanoic acid' },

    { smiles => 'ClC=1C=C(C=CC1)N=NC1=CC=C(C=C1)Cl', iupac => '(3-chlorophenyl)(4-chlorophenyl)diazene' }, # BBv3 P-68.3.1.3.2.1
    { smiles => 'C1(=CC=CC=C1)N=NC1=CC=C(C=C1)S(=O)(=O)O', iupac => '4-(phenyldiazenyl)benzene-1-sulfonic acid' }, # BBv3 P-68.3.1.3.2.2
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
