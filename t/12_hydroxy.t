#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'OCC(CO)(CO)CO', iupac => '2,2-bis(hydroxymethyl)propane-1,3-diol' },
    { smiles => 'CC(C)(C)O', iupac => '2-methylpropan-2-ol' }, # BBv2 P-63.1.2
    { smiles => 'ON1CC(CCC1)C#N', iupac => '1-hydroxypiperidine-3-carbonitrile' }, # BBv2 P-63.1.4

    # From BBv2 P-63.1.5
    { smiles => 'CC(C)S', iupac => 'propane-2-thiol' },
    { smiles => 'C(C)[SeH]', iupac => 'ethaneselenol' },
    { smiles => 'SCCCCS', iupac => 'butane-1,4-dithiol' },
    { smiles => 'SC1=CC=CC=C1', iupac => 'benzenethiol' },
    { smiles => 'SCCC(=O)O', iupac => '3-sulfanylpropanoic acid' },
    { smiles => 'SC1=C(C=CC=C1)O', iupac => '2-sulfanylphenol' },
    { smiles => 'OC(CS)C1CCC(C(C1)O)S', iupac => '5-(1-hydroxy-2-sulfanylethyl)-2-sulfanylcyclohexan-1-ol' },
    { smiles => 'SC(CC(=O)O)CS', iupac => '3,4-bis(sulfanyl)butanoic acid', AUTHOR => 1 },
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
