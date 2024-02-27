#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'OCC(CO)(CO)CO', iupac => '2,2-bis(hydroxymethyl)propane-1,3-diol' },
    { smiles => 'C1(CCCCC1)CCO', iupac => '2-cyclohexylethan-1-ol' }, # From BBv2 P-13.5.2
    { smiles => 'C1(CCCCC1)CO', iupac => 'cyclohexylmethanol' }, # From BBv2 P-15.6.1.1
    { smiles => 'ClCCCCC(CCO)C(C)O', iupac => '3-(4-chlorobutyl)pentane-1,4-diol' }, # From BBv2 P-44.1.1
    { smiles => 'CC(C)(C)O', iupac => '2-methylpropan-2-ol' }, # From BBv2 P-63.1.2
    { smiles => 'ON1CC(CCC1)C#N', iupac => '1-hydroxypiperidine-3-carbonitrile' }, # From BBv2 P-63.1.4
    { smiles => 'ClCCC(CC(C)O)O', iupac => '6-chlorohexane-2,4-diol' }, # From BBv3 P-92.1.4.1

    # From BBv3 P-63.1.5
    { smiles => 'CC(C)S', iupac => 'propane-2-thiol' },
    { smiles => 'C(C)[SeH]', iupac => 'ethaneselenol' },
    { smiles => 'SCCCCS', iupac => 'butane-1,4-dithiol' },
    { smiles => 'SC1=CC=CC=C1', iupac => 'benzenethiol' },
    { smiles => 'SCCC(=O)O', iupac => '3-sulfanylpropanoic acid' },
    { smiles => 'SC1=C(C=CC=C1)O', iupac => '2-sulfanylphenol' },
    { smiles => 'SC1=CC=C(C=C1)SSC=1C=C(C=CC1)S', iupac => '3-[(4-sulfanylphenyl)disulfanyl]benzene-1-thiol' },
    { smiles => 'OC(CS)C1CCC(C(C1)O)S', iupac => '5-(1-hydroxy-2-sulfanylethyl)-2-sulfanylcyclohexan-1-ol' },
    { smiles => 'SC(CC(=O)O)CS', iupac => '3,4-bis(sulfanyl)butanoic acid' },
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
