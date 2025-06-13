#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'NC1=CC2=CC(=CC2=CC=C1)O', iupac => '5-aminoazulen-2-ol' }, # From BBv3 P-15.1.8.1

    # From BBv2 P-25.0
    { smiles => 'C1=CC=CC2=CC=CC=C12', iupac => 'naphthalene' },
    { smiles => 'C=12C=CC=CC=1C=CC=C2', iupac => 'naphthalene' },

    # From BBv2 P-25.2.1
    { smiles => 'N1=CN=CC2=NC=CN=C12', iupac => 'pteridine' },
    { smiles => 'N1=CC=CC2=NC=CC=C12', iupac => '1,5-naphthyridine' },
    { smiles => 'C1=NC=CC2=CC=NC=C12', iupac => '2,7-naphthyridine' },
    { smiles => 'N1=CN=C2N=CNC2=C1', iupac => '7H-purine' },
    { smiles => 'N1N=CC2=CC=CC=C12', iupac => '1H-indazole' },
    { smiles => 'C12=C(C=CN2)C=CC=C1', iupac => '1H-indole' },
    { smiles => 'O1CC=CC2=C1C=CC=C2', iupac => '2H-1-benzopyran' },
    { smiles => 'C1OC=CC2=C1C=CC=C2', iupac => '1H-2-benzopyran' },
    { smiles => 'C1[Se]C=CC2=C1C=CC=C2', iupac => '1H-2-benzoselenopyran' },
    { smiles => 'N1=CC=CC2=CC=CC=C12', iupac => 'quinoline' },

    # From BBv3 P-25.2.1, Table 2.9
    { smiles => 'C1=CC=CC2=NC3=CC=CC=C3C=C12', iupac => 'acridine' },
    { smiles => 'C1=CC=CC2=[As]C3=CC=CC=C3C=C12', iupac => 'acridarsine' },
    { smiles => 'C1=CC=CC2=PC3=CC=CC=C3C=C12', iupac => 'acridophosphine' },
    { smiles => 'N1C=CC2=CC=CC=C12', iupac => '1H-indole' },
    { smiles => '[AsH]1C=CC2=CC=CC=C12', iupac => 'arsindole' },
    { smiles => 'P1C=CC2=CC=CC=C12', iupac => 'phosphindole' }, # Not sure if H prefix is not needed
    { smiles => 'C=1C=CN2C=CC=CC12', iupac => 'indolizine' },
    { smiles => 'C=1C=C[As]2C=CC=CC12', iupac => 'arsindolizine' },
    { smiles => 'C=1C=CP2C=CC=CC12', iupac => 'phosphindolizine' },
    { smiles => 'C=1NC=C2C=CC=CC12', iupac => '2H-isoindole' },
    { smiles => 'C=1[AsH]C=C2C=CC=CC12', iupac => 'isoarsindole' },
    { smiles => 'C=1PC=C2C=CC=CC12', iupac => 'isophosphindole' },
    { smiles => 'C1=[As]C=CC2=CC=CC=C12', iupac => 'isoarsinoline' },
    { smiles => 'C1=PC=CC2=CC=CC=C12', iupac => 'isophosphinoline' },
    { smiles => 'C1=CC=CC2=NC=C3C=CC=CC3=C12', iupac => 'phenanthridine' },
    { smiles => 'C1=CC=CC2=[As]C=C3C=CC=CC3=C12', iupac => 'arsanthridine' },
    { smiles => 'C1=CC=CC2=PC=C3C=CC=CC3=C12', iupac => 'phosphanthridine' },
    { smiles => '[As]1=CC=CC2=CC=CC=C12', iupac => 'arsinoline' },
    { smiles => 'P1=CC=CC2=CC=CC=C12', iupac => 'phosphinoline' },
    { smiles => 'C=1C=CC[As]2C=CC=CC12', iupac => '4H-arsinolizine' },
    { smiles => 'C=1C=CCP2C=CC=CC12', iupac => '4H-phosphinolizine' },

    { smiles => 'C1=CC=C2C(=C1)C(=C(P2(Cl)(Cl)Cl)Cl)Cl', iupac => '1,1,1,2,3-pentachlorophosphindole', AUTHOR => 1 }, # PubChem 2784508, why λ5 is not present?

    { smiles => 'C1=CC=C2C(C(C=CC2=C1)O)O', iupac => '1,2-dihydronaphthalene-1,2-diol' }, # PubChem 362
    { smiles => 'C12(C(C(C(C(C1(F)F)(F)F)(F)F)(F)F)(C(C(C(C2(F)F)(F)F)(F)F)(F)F)F)F', iupac => '1,1,2,2,3,3,4,4,4a,5,5,6,6,7,7,8,8,8a-octadecafluoronaphthalene' }, # PubChem 9386

    { smiles => 'C1=CC=CC=CC2=CC=CC=CC=C12', iupac => 'octalene' }, # From BBv2 P-25.1.2.3

    # From BBv2 P-25.2.2.4
    { smiles => 'C1=COC=CC2=C1C=CC=C2', iupac => '3-benzoxepine' },
    { smiles => 'O1C=CC2=C1C=CC=C2', iupac => '1-benzofuran' },
    { smiles => 'C=1OC=C2C1C=CC=C2', iupac => '2-benzofuran' },
    { smiles => 'C1=CC=COC=CC=CC=COC=CC=CC2=C1C=CC=C2', iupac => '5,12-benzodioxacyclooctadecine', AUTHOR => 1 },

    { smiles => 'C1=C(C=CC2=CC=CC=C12)CCCO', iupac => '3-(naphthalen-2-yl)propan-1-ol' }, # From P-15.6.1.2

    # From BBv3 P-44.4.1.5
    { smiles => '[SiH2]1O[SiH]=CC2=C1C=CC=C2', iupac => '1H-2,1,3-benzoxadisiline' },
    { smiles => '[SiH2]1OC=[SiH]C2=C1C=CC=C2', iupac => '1H-2,1,4-benzoxadisiline' },

    { smiles => 'C1=CC2=C(C=C1O)C(=CN2)CCN', iupac => '3-(2-aminoethyl)-1H-indol-5-ol' }, # serotonin

    # Guanine
    { smiles => 'O=C1c2ncnc2nc(N)N1', iupac => '2-amino-1,9-dihydro-6H-purin-6-one', AUTHOR => 1 }, # FIXME: fails due to aromaticity
    { smiles => 'NC=1NC(C=2N=CNC2N1)=O', iupac => '2-amino-1,9-dihydro-6H-purin-6-one', AUTHOR => 1 },
    { smiles => 'N1C(N)=NC=2N=CNC2C1=O', iupac => '2-amino-1,9-dihydro-6H-purin-6-one', AUTHOR => 1 },

    { smiles => 'CCN1C=NC2=C(N=CN=C21)N', iupac => '9-ethylpurin-6-amine' }, # PubChem 7
    { smiles => 'CCCCCNC1=C(C=NC2=CC=CC=C21)[N+](=O)[O-]', iupac => '3-nitro-N-pentylquinolin-4-amine' }, # PubChem 21355672
    { smiles => 'CC1=C(C2=C(C=C1)C=C(S2)O)CO', iupac => '7-(hydroxymethyl)-6-methyl-1-benzothiophen-2-ol' }, # PubChem 130803523
    { smiles => 'COC1=C(C2=C(C=C1)SC(=N2)N)Br', iupac => '4-bromo-5-methoxy-1,3-benzothiazol-2-amine' }, # PubChem 131985781

    { smiles => 'C1=CC=C2C(=C1)N=CS2', iupac => '1,3-benzothiazole' }, # PubChem 7222
    { smiles => 'C1=CC=C2C(=C1)C=NS2', iupac => '1,2-benzothiazole' }, # PubChem 9225
    { smiles => 'C1=CC2=CSN=C2C=C1', iupac => '2,1-benzothiazole' }, # PubChem 638008

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
