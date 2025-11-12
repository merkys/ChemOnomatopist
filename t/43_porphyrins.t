#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'C12=CC=C(N1)C=C1C=CC(=N1)C=C1C=CC(N1)=CC=1C=CC(N1)=C2', iupac => 'porphyrin' },

    { smiles => 'C1(=CC=CC=C1)C=1C2=CC=C(N2)C=C2C=CC(C(=C3C=CC(=CC=4C=CC1N4)N3)C3=CC=CC=C3)=N2', iupac => '5,15-diphenylporphyrin' },

    { smiles => 'BrC=1C2=C(C3=CC=C(N3)C(=C3C(=C(C(C(=C4C=CC(=C(C(C1Br)=N2)C2=CC=CC=C2)N4)C4=CC=CC=C4)=N3)Br)Br)C3=CC=CC=C3)C3=CC=CC=C3',
      iupac  => '7,8,17,18-tetrabromo-5,10,15,20-tetraphenylporphyrin' }, # COD entry 2241690
    { smiles => 'C(C)C1=C2NC(=C1CC)C(=C1C(=C(C(=N1)C=C1C(=C(C(N1)=CC=1C(=C(C(N1)=C2[N+](=O)[O-])CC)CC)CC)CC)CC)CC)[N+](=O)[O-]',
      iupac  => '2,3,7,8,12,13,17,18-octaethyl-5,20-dinitroporphyrin' }, # COD entry 4032498
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
