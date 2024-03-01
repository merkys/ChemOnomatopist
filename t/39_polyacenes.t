#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-25.1.2.1
    { smiles => 'C1=CC=CC2=CC3=CC4=CC=CC=C4C=C3C=C12', iupac => 'tetracene' },
    { smiles => 'C1=CC=CC2=CC3=CC4=CC5=CC=CC=C5C=C4C=C3C=C12', iupac => 'pentacene' },

    { smiles => 'CC1=C2C=C3C=CC=CC3=CC2=CC2=CC3=CC=CC=C3C=C12', iupac => '6-methylpentacene', AUTHOR => 1 }, # From Wikipedia Pentacene
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
