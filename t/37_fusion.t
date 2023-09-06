#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-25.3.1.3
    { smiles => '[Se]1C=CC2=C1[Se]C=C2', iupac => 'selenopheno[2,3-b]selenophene', AUTHOR => 1 },
    { smiles => '[Se]1C=2C(C=C1)=C[Se]C2', iupac => 'selenopheno[3,4-b]selenophene', AUTHOR => 1 },
    { smiles => '[Se]1C2=C(C=C1)[Se]C=C2', iupac => 'selenopheno[3,2-b]selenophene' },

    # From BBv2 P-25.3.2.4
    { smiles => 'S1CC=CSC=2C1=COC2', iupac => '2H-[1,4]dithiepino[2,3-c]furan', AUTHOR => 1 },
    { smiles => 'O1CC=C2OC=CC=C21', iupac => '2H-furo[3,2-b]pyran', AUTHOR => 1 },
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
