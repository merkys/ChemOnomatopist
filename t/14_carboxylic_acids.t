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

    # From BBv3 P-65.1.2.2.3
    { smiles => 'C(=O)(O)CC(CC(=O)O)CCCC(=O)O', iupac => '3-(carboxymethyl)heptanedioic acid', AUTHOR => 1 },

    # From BBv3 P-65.1.5.1
    { smiles => 'C(CCCCC)(O)=S', iupac => 'hexanethioic O-acid' },
    { smiles => 'C(CCCCC)(=[Se])S', iupac => 'hexaneselenothioic S-acid' },
    { smiles => 'C(CCCCC(=S)S)(=S)S', iupac => 'hexanebis(dithioic acid)' },
    { smiles => 'N1(CCCCC1)C(=S)S', iupac => 'piperidine-1-carbodithioic acid' },
    { smiles => 'C1CCCCC1C(=S)[SeH]', iupac => 'cyclohexanecarboselenothioic Se-acid' },
    { smiles => 'C(C)(=S)C1=CC=C(C(=O)O)C=C1', iupac => '4-(ethanethioyl)benzoic acid', AUTHOR => 1 },
    { smiles => 'O=C(CCC(=O)O)S', iupac => '4-oxo-4-sulfanylbutanoic acid', AUTHOR => 1 },
    { smiles => 'OC(CCC(=O)O)=S', iupac => '4-hydroxy-4-sulfanylidenebutanoic acid', AUTHOR => 1 },
    { smiles => 'O=C(C(=O)O)S', iupac => 'oxo(sulfanyl)acetic acid', AUTHOR => 1 },
    { smiles => 'OC(=S)C1=CC(=NC=C1)C(=O)O', iupac => '4-(hydroxycarbonothioyl)pyridine-2-carboxylic acid' },
    { smiles => 'SC(=O)C1=CC(=NC=C1)C(=O)O', iupac => '4-(sulfanylcarbonyl)pyridine-2-carboxylic acid' },
    { smiles => 'C(CC)(=N)S', iupac => 'propanimidothioic acid', AUTHOR => 1 },
    { smiles => 'C(CCC)(=NN)[SeH]', iupac => 'butanehydrazonoselenoic acid', AUTHOR => 1 },
    { smiles => 'SN=C(O)C1CCCC1', iupac => 'N-sulfanylcyclopentanecarboximidic acid', AUTHOR => 1 },
    { smiles => 'ON=C([SeH])C1CCCCC1', iupac => 'N-hydroxycyclohexanecarboximidoselenoic acid', AUTHOR => 1 },
    { smiles => 'NC(=CC(=S)S)SCC', iupac => '3-amino-3-(ethylsulfanyl)prop-2-ene(dithioic acid)', AUTHOR => 1 },

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
    { smiles => 'S=C(CCC(=O)O)OS', iupac => '4-sulfanylidene-4-(sulfanyloxy)butanoic acid', AUTHOR => 1 },
    { smiles => 'C(=S)(OS)CCC(=O)O', iupac => '3-(dithiocarbonoperoxoyl)propanoic acid', AUTHOR => 1 },
    { smiles => 'C(=S)(OS)C(=O)O', iupac => '(dithiocarbonoperoxoyl)formic acid', AUTHOR => 1 },
    { smiles => 'OSC(C(=O)O)=O', iupac => '(hydroxysulfanyl)oxoacetic acid', AUTHOR => 1 },
    { smiles => 'OSC(=O)C1CCC(CC1)C(=O)O', iupac => '4-[(hydroxysulfanyl)carbonyl]cyclohexanecarboxylic acid', AUTHOR => 1 },

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
