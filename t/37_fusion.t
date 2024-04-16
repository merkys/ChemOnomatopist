#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-25.3.1.3
    { smiles => '[Se]1C=CC2=C1[Se]C=C2', iupac => 'selenopheno[2,3-b]selenophene' },
    { smiles => '[Se]1C=2C(C=C1)=C[Se]C2', iupac => 'selenopheno[3,4-b]selenophene', AUTHOR => 1 },
    { smiles => '[Se]1C2=C(C=C1)[Se]C=C2', iupac => 'selenopheno[3,2-b]selenophene' },

    # From BBv2 P-25.3.2.4
    { smiles => 'S1CC=CSC=2C1=COC2', iupac => '2H-[1,4]dithiepino[2,3-c]furan', AUTHOR => 1 },
    { smiles => 'O1CC=C2OC=CC=C21', iupac => '2H-furo[3,2-b]pyran', AUTHOR => 1 },
    { smiles => 'N1=CC=CC2=C1C=NOC2', iupac => '5H-pyrido[2,3-d][1,2]oxazine' },
    { smiles => 'O1COC2=C1C=CO2', iupac => '2H-furo[2,3-d][1,3]dioxole' },
    { smiles => 'O1P=CC2=C1OCO2', iupac => '5H-[1,3]dioxolo[4,5-d][1,2]oxaphosphole' },
    { smiles => 'S1C=NC2=C1N=C[Se]2', iupac => '[1,3]selenazolo[5,4-d][1,3]thiazole' },
    { smiles => 'O1C2=C(SC=C1)[Se]C=CO2', iupac => '[1,4]oxaselenino[2,3-b][1,4]oxathiine' },
    { smiles => 'N1=CC=NC=2C1=CN=NC2', iupac => 'pyrazino[2,3-d]pyridazine' },
    { smiles => 'O1SNC2=C1ONS2', iupac => '3H,5H-[1,3,2]oxathiazolo[4,5-d][1,2,3]oxathiazole' },

    { smiles => 'S1C=2N(C=C1)C=CN2', iupac => 'imidazo[2,1-b][1,3]thiazole' }, # From BBv2 P-25.3.2.5.1

    # From BBv2 P-25.3.3.1.2
    { smiles => 'O1C=2C(=CC=C1)C=CC2', iupac => 'cyclopenta[b]pyran', AUTHOR => 1 },
    { smiles => 'S1COC=2NC=CC21', iupac => '2H,4H-[1,3]oxathiolo[5,4-b]pyrrole', AUTHOR => 1 },
    { smiles => 'O1C2=C(C=C1)C=CS2', iupac => 'thieno[2,3-b]furan' },
    { smiles => 'N1C=NC2=C1C=CS2', iupac => '1H-thieno[2,3-d]imidazole', AUTHOR => 1 },
    { smiles => 'N=1N2C(N=CC1)=NC=C2', iupac => 'imidazo[1,2-b][1,2,4]triazine' },
    { smiles => 'O1COC2=C1N=CN2', iupac => '2H,4H-[1,3]dioxolo[4,5-d]imidazole', AUTHOR => 1 },

    { smiles => 'CN(C=1C2=C(N=CN1)NC=C2)C', iupac => 'N,N-dimethyl-7H-pyrrolo[2,3-d]pyrimidin-4-amine', AUTHOR => 1 }, # PubChem 1861
    { smiles => 'C1(=CC=CC=C1)C1=CC2=C(N=CN=C2N)N1C1=CC=CC=C1', iupac => '6,7-diphenylpyrrolo[2,3-d]pyrimidin-4-amine' }, # PubChem 49845043
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
