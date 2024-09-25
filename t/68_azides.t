#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-15.2.1.1
    { smiles => 'N(=[N+]=[N-])C1=CC=CC=C1', iupac => 'azidobenzene' },

    # From BBv3 P-61.7
    { smiles => 'N(=[N+]=[N-])CCC1=CC=CC=C1', iupac => '(2-azidoethyl)benzene' },
    { smiles => 'N(=[N+]=[N-])C=1C(=CC2=CC=CC=C2C1)S(=O)(=O)O', iupac => '3-azidonaphthalene-2-sulfonic acid' },
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
