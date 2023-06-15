#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'NC(=N)N', iupac => 'guanidine' },

    # From BBv2 P-66.4.1.2.1.2
    { smiles => 'CN(C(=NC1=CC=CC=C1)N(C)C)C', iupac => "N,N,N',N'-tetramethyl-N''-phenylguanidine", AUTHOR => 1 },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    eval { is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
}
