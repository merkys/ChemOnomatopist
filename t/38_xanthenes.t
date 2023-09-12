#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'C1=CC=CC=2OC3=CC=CC=C3CC12', iupac => '9H-xanthene' }, # From BBv2 P-25.2.1
    { smiles => 'C1=CC=CC=2OC3=CC=CC=C3OC12', iupac => 'oxanthrene' }, # From BBv2 P-25.2.2.2

    { smiles => 'C1=CC=CC=2OC3=CC=CC=C3C(C12)=O', iupac => '9H-xanthen-9-one' }, # From Wikipedia Xanthone
    { smiles => 'OC1=CC(=CC=2OC3=CC(=CC(=C3C(C12)=O)C)OC)OC', iupac => '1-hydroxy-3,6-dimethoxy-8-methyl-9H-xanthen-9-one' }, # From Wikipedia Lichexanthone
    { smiles => 'C1=CC=CC=2SC3=CC=CC=C3C(C12)=O', iupac => '9H-thioxanthen-9-one' }, # From Wikipedia Thioxanthone
    { smiles => 'CC(C)C1=CC=CC=2C(C3=CC=CC=C3SC12)=O', iupac => '4-(propan-2-yl)-9H-thioxanthen-9-one' }, # From Wikipedia Isopropylthioxanthone
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
