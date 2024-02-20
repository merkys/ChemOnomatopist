#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # BBv2 P-25.1.2.6
    { smiles => 'C1=CC=CC2=CC=C3C=CC4=CC=C5C=CC6=CC=CC=C6C5=C4C3=C12', iupac => 'hexahelicene' },

    { smiles => 'C1=CC=CC2=CC=C3C=CC4=CC=C5C=CC6=CC=C7C=CC=CC7=C6C5=C4C3=C12', iupac => 'heptahelicene' },
    { smiles => 'C1=CC=CC2=CC=C3C=CC4=CC=C5C=CC6=CC=C7C=CC8=CC=CC=C8C7=C6C5=C4C3=C12', iupac => 'octahelicene' },
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
