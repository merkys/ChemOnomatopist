#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'C1=CC=CC=2C3=CC=CC=C3C=CC12', iupac => 'phenanthrene' },

    # From BBv2 P-25.2.1
    { smiles => 'N1=CC=CC2=CC=C3N=CC=CC3=C12', iupac => '1,7-phenanthroline' },
    { smiles => 'C1=CC=CC2=NC=C3C=CC=CC3=C12', iupac => 'phenanthridine' },

    { smiles => 'COC=1C(=CC=2C=C(C3=CC(=CC=C3C2C1OC)O)OC)O', iupac => '3,4,9-trimethoxyphenanthrene-2,7-diol' }, # From Wikipedia Gymnopusin
    { smiles => 'C1=CC=CC=2C3=CC=CC=C3C(C(C12)=O)=O', iupac => 'phenanthrene-9,10-dione' }, # From Wikipedia Phenanthrenequinone
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
