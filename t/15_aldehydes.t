#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-66.6.1.1.1
    { smiles => 'CCCCC=O', iupac => 'pentanal' },
    { smiles => 'O=CCCCC=O', iupac => 'pentanedial' },

    { smiles => 'O=CCC(C=O)CCC=O', iupac => 'butane-1,2,4-tricarbaldehyde' }, # From BBv2 P-66.6.1.1.2

    { smiles => 'CN(N=NN(C)C)C=O', iupac => '1,4,4-trimethyltetraaz-2-ene-1-carbaldehyde', AUTHOR => 1 }, # From BBv2 P-66.6.1.1.3
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    eval { is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
}
