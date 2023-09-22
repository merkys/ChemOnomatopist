#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'ClC1=NC=CC(=C1)OC1=NC=C(C=C1)Cl', iupac => '2-chloro-4-[(5-chloropyridin-2-yl)oxy]pyridine' }, # BBv2 P-15.3.2.4.2
    { smiles => 'C(C1=CC=CC=C1)C1=NC=CC=C1', iupac => '2-benzylpyridine' }, # BBv2 P-29.6.1
    { smiles => 'SC1=CC=C(C=C1)SSC=1C=C(C=CC1)S', iupac => '3-[(4-sulfanylphenyl)disulfanyl]benzene-1-thiol', AUTHOR => 1 }, # BBv2 P-63.1.5
    { smiles => 'C1(=CC=CC=C1)SC1CCNCC1', iupac => '4-(phenylsulfanyl)piperidine', AUTHOR => 1 }, # BBv2 P-63.2.5

    # From BBv2 P-63.2.4.2
    { smiles => 'C1(CCCCC1)OC1=CC=CC=C1', iupac => '(cyclohexyloxy)benzene', AUTHOR => 1 }, # FIXME: Very close
    { smiles => 'N1=CC(=CC=C1)OC1=NC=CN=C1', iupac => '2-[(pyridin-3-yl)oxy]pyrazine' },

    { smiles => 'c1ncccc1[C@@H]2CCCN2C', iupac => '3-[(2S)-1-methylpyrrolidin-2-yl]pyridine', AUTHOR => 1 }, # nicotine
    { smiles => 'CC1C(=O)NC(C(=O)N1C(C)C(C)C)C2CCCCC2', iupac => '3-cyclohexyl-6-methyl-1-(3-methylbutan-2-yl)piperazine-2,5-dione' }, # PubChem 64959818
    { smiles => 'C1=CC(=C(N=C1)C2=NC(=NS2)N)Br', iupac => '5-(3-bromopyridin-2-yl)-1,2,4-thiadiazol-3-amine', AUTHOR => 1 }, # PubChem 107526369
    { smiles => 'C1=CC(=NC(=C1N)C2=CC(=C(C=C2Cl)Cl)Cl)C(=O)O', iupac => '5-amino-6-(2,4,5-trichlorophenyl)pyridine-2-carboxylic acid', AUTHOR => 1 }, # PubChem 133086582
    { smiles => 'CC1CCC(CC(C1)C)C2CCCCC2', iupac => '1-cyclohexyl-3,5-dimethylcycloheptane' }, # PubChem 149225482
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
