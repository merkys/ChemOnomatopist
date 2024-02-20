#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-63.6
    { smiles => 'C(C)S(=O)CCCC', iupac => '1-(ethanesulfinyl)butane' },
    { smiles => 'C(C)[Se](=O)C1=CC=CC=C1', iupac => '(ethaneseleninyl)benzene' },

    { smiles => 'C1(CCCCC1)S(=O)N1CC(OCC1)C(=O)O', iupac => '4-(cyclohexanesulfinyl)morpholine-2-carboxylic acid' }, # BBv2 P-65.4.1
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
