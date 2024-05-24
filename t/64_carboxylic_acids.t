#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 Table 4.3. Suffixes were turned to full names by prepending them with 'benzene'
    { smiles => 'C1(=CC=CC=C1)C(=O)O', iupac => 'benzenecarboxylic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=O)OO', iupac => 'benzenecarboperoxoic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=S)OO', iupac => 'benzenecarboperoxothioic OO-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=[Se])OO', iupac => 'benzenecarboperoxoselenoic acid' },
    { smiles => 'C1(=CC=CC=C1)C(=O)SO', iupac => 'benzenecarbo(thioperoxoic) SO-acid' },
    { smiles => 'C1(=CC=CC=C1)C(=O)OS', iupac => 'benzenecarbo(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=O)[Se]O', iupac => 'benzenecarbo(selenoperoxoic) SeO-acid' },
    { smiles => 'C1(=CC=CC=C1)C(=O)O[SeH]', iupac => 'benzenecarbo(selenoperoxoic) OSe-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=S)SO', iupac => 'benzenecarbothio(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=S)OS', iupac => 'benzenecarbothio(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=[Se])OS', iupac => 'benzenecarboseleno(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=S)[Se]O', iupac => 'benzenecarbo(selenoperoxo)thioic SeO-acid' },
    { smiles => 'C1(=CC=CC=C1)C(=S)O[SeH]', iupac => 'benzenecarbo(selenoperoxo)thioic OSe-acid' },
    { smiles => 'C1(=CC=CC=C1)C(=S)SS', iupac => 'benzenecarbo(dithioperoxo)thioic acid' },
    { smiles => 'C1(=CC=CC=C1)C(=[Se])[Se][SeH]', iupac => 'benzenecarbo(diselenoperoxo)selenoic acid' },
    { smiles => 'C1(=CC=CC=C1)C(=[Te])[Te][TeH]', iupac => 'benzenecarbo(ditelluroperoxo)telluroic acid' },
    { smiles => 'C1(=CC=CC=C1)C(O)=S', iupac => 'benzenecarbothioic O-acid' },
    { smiles => 'C1(=CC=CC=C1)C(S)=O', iupac => 'benzenecarbothioic S-acid' },
    { smiles => 'C1(=CC=CC=C1)C(=S)S', iupac => 'benzenecarbodithioic acid' },
    { smiles => 'C1(=CC=CC=C1)C(=[Se])S', iupac => 'benzenecarboselenothioic S-acid' },
    { smiles => 'C1(=CC=CC=C1)C(=S)[SeH]', iupac => 'benzenecarboselenothioic Se-acid' },
    { smiles => 'C1(=CC=CC=C1)C(=[Se])[SeH]', iupac => 'benzenecarbodiselenoic acid' },
    { smiles => 'C1(=CC=CC=C1)C(=[Te])[SeH]', iupac => 'benzenecarboselenotelluroic Se-acid' },
    { smiles => 'C1(=CC=CC=C1)C(=[Te])[TeH]', iupac => 'benzenecarboditelluroic acid' },
    { smiles => 'C1(=CC=CC=C1)C(O)=N', iupac => 'benzenecarboximidic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)OO', iupac => 'benzenecarboximidoperoxoic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)SO', iupac => 'benzenecarboximido(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)OS', iupac => 'benzenecarboximido(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)SS', iupac => 'benzenecarbo(dithioperoxo)imidic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)[Se]S', iupac => 'benzenecarboximido(selenothioperoxoic) SeS-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)S', iupac => 'benzenecarboximidothioic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)[SeH]', iupac => 'benzenecarboximidoselenoic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)[TeH]', iupac => 'benzenecarboximidotelluroic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(O)=NN', iupac => 'benzenecarbohydrazonic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=NN)OO', iupac => 'benzenecarbohydrazonoperoxoic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=NN)SO', iupac => 'benzenecarbohydrazono(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=NN)OS', iupac => 'benzenecarbohydrazono(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=NN)[Te][TeH]', iupac => 'benzenecarbo(ditelluroperoxo)hydrazonic acid', AUTHOR => 1 },
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
