#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-44.1.1
    { smiles => 'C1(CCCCC1)CCC(=O)O', iupac => '3-cyclohexylpropanoic acid' },
    { smiles => 'C(CC)C=1C=C(C(=O)O)C=CC1', iupac => '3-propylbenzoic acid' },
    { smiles => 'ClCCCCC(CCO)C(C)O', iupac => '3-(4-chlorobutyl)pentane-1,4-diol' },
    { smiles => 'N(N)C(=O)O', iupac => 'hydrazinecarboxylic acid', AUTHOR => 1 },
    { smiles => '[SiH3]CCC(=O)O', iupac => '3-silylpropanoic acid' },
    { smiles => 'C(C)[SiH2]C(=O)O', iupac => 'ethylsilanecarboxylic acid', AUTHOR => 1 },
    { smiles => 'C(CCC)OCCOCC(SCCSCC(=O)O)(CCSCCSCC(=O)O)COCCOCCCC', iupac => '7,7-bis[(2-butoxyethoxy)methyl]-3,6,10,13-tetrathiapentadecane-1,15-dioic acid', AUTHOR => 1 },
    { smiles => '[SiH2](CC[SiH3])CC[SiH3]', iupac => '[silanediyldi(ethane-2,1-diyl)]bis(silane)', AUTHOR => 1 },

    # From BBv3 P-44.1.2.1
    { smiles => 'C[Si](C)(C)C', iupac => 'tetramethylsilane', AUTHOR => 1 },
    { smiles => 'CP[SiH3]', iupac => 'methyl(silyl)phosphane', AUTHOR => 1 },
    { smiles => 'C(C)(C)(C)[Si](OCC1OC1)(C)C', iupac => 'tert-butyldi(methyl)(oxiranylmethoxy)silane', AUTHOR => 1 },
    { smiles => 'C(=O)(O)CC[SiH2][SiH2]C(=O)O', iupac => '2-(2-carboxyethyl)disilane-1-carboxylic acid', AUTHOR => 1 },
    { smiles => 'O1C(=CC2=C1C=CC=C2)P', iupac => '(1-benzofuran-2-yl)phosphane' },
    { smiles => 'C[Si](N1C=NC=C1)(C)C', iupac => '1-(trimethylsilyl)-1H-imidazole', AUTHOR => 1 },
    { smiles => 'C(#N)C1=PC=CC(=C1)C1CC(OCC1)C#N', iupac => '4-(2-cyanophosphinin-4-yl)oxane-2-carbonitrile' },
    { smiles => 'P1=C(C=CC=C1)PC=1OC=CC1', iupac => '2-[(phosphinin-2-yl)phosphanyl]furan', AUTHOR => 1 },
    { smiles => 'O1CC(=CC=C1)NNC1[SiH2]CCC1', iupac => '1-(2H-pyran-3-yl)-2-(silolan-2-yl)hydrazine' },
    { smiles => 'C(SCSCSCSC)NCOCOCOCOC', iupac => 'N-(2,4,6,8-tetrathianonan-1-yl)-2,4,6,8-tetraoxanonan-1-amine', AUTHOR => 1 },
    { smiles => 'C([SiH2]C[SiH2]C[SiH2]C[SiH2]C)C1CC(COC1)COCOCOCOC', iupac => '1-[5-(2,4,6,8-tetrasilanonan-1-yl)oxan-3-yl]-2,4,6,8-tetraoxanonane', AUTHOR => 1 },

    # From BBv3 P-44.1.2.2
    { smiles => 'C(CCCCCC)C1=CC=CC=C1', iupac => 'heptylbenzene' },
    { smiles => 'C(=C)C1CCCCC1', iupac => 'ethenylcyclohexane' },
    { smiles => 'C(C1=CC=CC=C1)C1=CC=CC=C1', iupac => '1,1\'-methylenedibenzene', AUTHOR => 1 },
    { smiles => 'C(=CC1CCCCC1)C1CCCCC1', iupac => '1,1\'-(ethene-1,2-diyl)dicyclohexane', AUTHOR => 1 },
    { smiles => 'N(N)C1=NC=CC=C1', iupac => '2-hydrazinylpyridine', AUTHOR => 1 },
    { smiles => 'N(N)C=1NCCN1', iupac => '2-hydrazinyl-4,5-dihydro-1H-imidazole' },
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
