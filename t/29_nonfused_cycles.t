#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'ClC1=NC=CC(=C1)OC1=NC=C(C=C1)Cl', iupac => '2-chloro-4-[(5-chloropyridin-2-yl)oxy]pyridine', AUTHOR => 1}, # BBv2 P-15.3.2.4.2
    { smiles => 'c1ncccc1[C@@H]2CCCN2C', iupac => '3-[(2S)-1-methylpyrrolidin-2-yl]pyridine', AUTHOR => 1 }, # nicotine
    { smiles => 'C1=CC(=NC(=C1N)C2=CC(=C(C=C2Cl)Cl)Cl)C(=O)O', iupac => '5-amino-6-(2,4,5-trichlorophenyl)pyridine-2-carboxylic acid', AUTHOR => 1 }, # PubChem 133086582
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
