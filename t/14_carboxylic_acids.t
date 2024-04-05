#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-65.1.2
    { smiles => 'CCCC(=O)O', iupac => 'butanoic acid' },
    { smiles => 'OC(=O)CCCCCCCCCCC(=O)O', iupac => 'dodecanedioic acid' },

    # From BBv2 P-65.1.2.2.1
    { smiles => 'OC(=O)CCC(C(=O)O)CCC(=O)O', iupac => 'pentane-1,3,5-tricarboxylic acid' },
    { smiles => 'C(C(=O)O)(C(=O)O)C(C(=O)O)(C(=O)O)', iupac => 'ethane-1,1,2,2-tetracarboxylic acid' },

    # From BBv2 P-65.1.5.2
    { smiles => 'C(C)(O)=S', iupac => 'ethanethioic O-acid' },
    { smiles => 'C(S)=O', iupac => 'methanethioic S-acid', AUTHOR => 1 },
    { smiles => 'O=C(CCC(=O)O)S', iupac => '4-oxo-4-sulfanylbutanoic acid', AUTHOR => 1 },
    { smiles => 'OC(C(=O)O)=S', iupac => 'hydroxy(sulfanylidene)acetic acid', AUTHOR => 1 },
    { smiles => '[SeH]C(=O)C1=CC=C(C(=O)O)C=C1', iupac => '4-(selanylcarbonyl)benzoic acid' },
    { smiles => 'C=1(C(=CC=CC1)C(=S)S)C(=S)S', iupac => 'benzene-1,2-dicarbodithioic acid' },
    { smiles => 'C(C(=S)S)(=S)S', iupac => 'ethanebis(dithioic acid)' },

    # From BBv3 P-65.1.5.3
    { smiles => 'CC(=O)OS', iupac => 'ethane(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => 'c1ccccc1C(=O)SO', iupac => 'benzenecarbo(thioperoxoic) SO-acid' },
    { smiles => 'C1=C(C=CC2=CC=CC=C12)C(OO)=S', iupac => 'naphthalene-2-carboperoxothioic acid' },

    { smiles => 'C(=O)O', iupac => 'formic acid', AUTHOR => 1 },

    { smiles => 'C1=CC2=C(C=C1C(=O)O)SC(=N2)C(F)F', iupac => '2-(difluoromethyl)-1,3-benzothiazole-6-carboxylic acid' }, # PubChem 84692459
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
