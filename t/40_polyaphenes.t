#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'OC1=C2C(C=3C=CC=CC3C(C2=C2C(C3=CC=CC=C3C(C2=C1O)=O)=O)=O)=O', iupac => '6,7-dihydroxypentaphene-5,8,13,14-tetrone' }, # PubChem 5379520
    { smiles => 'CC1=C2C(C=3C=CC=CC3C(C2=C2C(C3=CC=CC=C3C(C2=C1)=O)=O)=O)=O', iupac => '6-methylpentaphene-5,8,13,14-tetrone' }, # PubChem 93949195
    { smiles => 'FC1=C(C2=C(C3=C4C(=C5C(=C6C(=C(C(=C(C6=C(C5=C(C4=C(C(=C3C(=C2C(=C1F)F)F)F)F)F)F)F)F)F)F)F)F)F)C1=CC=CC=C1', iupac => '2,3,4,5,6,7,8,9,10,11,12,13,14,15,16-pentadecafluoro-1-phenylhexaphene', AUTHOR => 'flaky' }, # PubChem 121333830
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
