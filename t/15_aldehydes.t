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

    # From BBv2 P-66.6.3
    { smiles => 'C(C)=S', iupac => 'ethanethial' },
    { smiles => 'C1(=CC=CC=C1)C=S', iupac => 'benzenecarbothialdehyde' },
    { smiles => 'C(CCCCC)=[Se]', iupac => 'hexaneselenal' },
    { smiles => 'C(CCCC=S)=S', iupac => 'pentanedithial' },
    { smiles => 'C(=S)C1=CC=C(C(=O)O)C=C1', iupac => '4-(methanethioyl)benzoic acid', AUTHOR => 1 },
    { smiles => 'C(=[Se])C1CCC(CC1)C(=O)O', iupac => '4-(methaneselenoyl)cyclohexane-1-carboxylic acid', AUTHOR => 1 },
    { smiles => 'S=C1CCC(CC1)C=[Se]', iupac => '4-sulfanylidenecyclohexane-1-carboselenaldehyde' },

    # From BBv2 P-66.6.4
    { smiles => 'O=C(CC=O)C', iupac => '3-oxobutanal' },
    { smiles => 'C=C(C=O)CCCC', iupac => '2-methylidenehexanal' },
    { smiles => 'OC1=C(C=O)C=CC=C1', iupac => '2-hydroxybenzaldehyde' },
    { smiles => 'OCC1=CC=C(O1)C=O', iupac => '5-(hydroxymethyl)furan-2-carbaldehyde' },
    { smiles => 'O(C1=CC=CC=C1)CC=O', iupac => 'phenoxyacetaldehyde', AUTHOR => 1 },
    { smiles => 'FC=1C(=C(C=O)C=CC1)C', iupac => '3-fluoro-2-methylbenzaldehyde' },
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
