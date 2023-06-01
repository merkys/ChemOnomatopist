#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-25.0
    { smiles => 'C1=CC=CC2=CC=CC=C12', iupac => 'naphtalene' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    eval { is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
}
