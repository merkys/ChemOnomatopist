#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-25.0
    { smiles => 'C1=CC=CC2=CC=CC=C12', iupac => 'naphthalene' },
    { smiles => 'C=12C=CC=CC=1C=CC=C2', iupac => 'naphthalene' },

    # From BBv2 P-25.2.1
    { smiles => 'N1=CN=CC2=NC=CN=C12', iupac => 'pteridine' },
    { smiles => 'N1=CN=C2N=CNC2=C1', iupac => 'purine' },
    { smiles => 'N1N=CC2=CC=CC=C12', iupac => '1H-indazole' },
    { smiles => 'C12=C(C=CN2)C=CC=C1', iupac => '1H-indole' },
    { smiles => 'O1CC=CC2=C1C=CC=C2', iupac => '2H-1-benzopyran' },
    { smiles => 'C1OC=CC2=C1C=CC=C2', iupac => '1H-2-benzopyran' },
    { smiles => 'C1[Se]C=CC2=C1C=CC=C2', iupac => '1H-2-benzoselenopyran' },
    { smiles => 'N1=CC=CC2=CC=CC=C12', iupac => 'quinoline' },

    { smiles => 'P1C=CC2=CC=CC=C12', iupac => 'phosphindole' }, # Not sure if H prefix is not needed
    { smiles => 'C1=CC=C2C(=C1)C(=C(P2(Cl)(Cl)Cl)Cl)Cl', iupac => '1,1,1,2,3-pentachlorophosphindole' }, # PubChem 2784508

    { smiles => 'C1=CC=C2C(C(C=CC2=C1)O)O', iupac => 'naphthalene-1,2-diol' }, # PubChem 362 has 1,2-dihydronaphthalene-1,2-diol, which is clearly incorrect
    { smiles => 'C12(C(C(C(C(C1(F)F)(F)F)(F)F)(F)F)(C(C(C(C2(F)F)(F)F)(F)F)(F)F)F)F', iupac => '1,1,2,2,3,3,4,4,4a,5,5,6,6,7,7,8,8,8a-octadecafluoronaphthalene', AUTHOR => 1 }, # PubChem 9386

    { smiles => 'C1=CC=CC=CC2=CC=CC=CC=C12', iupac => 'octalene' }, # From BBv2 P-25.1.2.3

    # From BBv2 P-25.2.2.4
    { smiles => 'C1=COC=CC2=C1C=CC=C2', iupac => '3-benzoxepine', AUTHOR => 1 },
    { smiles => 'O1C=CC2=C1C=CC=C2', iupac => '1-benzofuran', AUTHOR => 1 },
    { smiles => 'C=1OC=C2C1C=CC=C2', iupac => '2-benzofuran', AUTHOR => 1 }, # FIXME: Fails to be detected as aromatic
    { smiles => 'C1=CC=COC=CC=CC=COC=CC=CC2=C1C=CC=C2', iupac => '5,12-benzodioxacyclooctadecine', AUTHOR => 1 },

    { smiles => 'C1=C(C=CC2=CC=CC=C12)CCCO', iupac => '3-(naphthalen-2-yl)propan-1-ol', AUTHOR => 1 }, # From P-15.6.1.2

    { smiles => 'C1=CC2=C(C=C1O)C(=CN2)CCN', iupac => '3-(2-aminoethyl)-1H-indol-5-ol' }, # serotonin

    # Guanine
    { smiles => 'O=C1c2ncnc2nc(N)N1', iupac => '2-amino-1,9-dihydro-6H-purin-6-one', AUTHOR => 1 }, # FIXME: fails due to aromaticity
    { smiles => 'NC=1NC(C=2N=CNC2N1)=O', iupac => '2-amino-1,9-dihydro-6H-purin-6-one', AUTHOR => 1 },
    { smiles => 'N1C(N)=NC=2N=CNC2C1=O', iupac => '2-amino-1,9-dihydro-6H-purin-6-one', AUTHOR => 1 },

    { smiles => 'CCN1C=NC2=C(N=CN=C21)N', iupac => '9-ethylpurin-6-amine', AUTHOR => 1 }, # PubChem 7
    { smiles => 'CCCCCNC1=C(C=NC2=CC=CC=C21)[N+](=O)[O-]', iupac => '3-nitro-N-pentylquinolin-4-amine', AUTHOR => 1 }, # PubChem 21355672
    { smiles => 'CC1=C(C2=C(C=C1)C=C(S2)O)CO', iupac => '7-(hydroxymethyl)-6-methyl-1-benzothiophen-2-ol' }, # PubChem 130803523
    { smiles => 'COC1=C(C2=C(C=C1)SC(=N2)N)Br', iupac => '4-bromo-5-methoxy-1,3-benzothiazol-2-amine' }, # PubChem 131985781

    { smiles => 'C1(=CC=C(C=2C(=CC=C(C12)C(=O)O)C(=O)O)C(=O)O)C(=O)O', iupac => 'naphthalene-1,4,5,8-tetracarboxylic acid' },
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
