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
    { smiles => '', iupac => 'benzenecarbo(selenoperoxoic) SeO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbo(selenoperoxoic) OSe-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbothio(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbothio(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarboseleno(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbo(selenoperoxo)thioic SeO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbo(selenoperoxo)thioic OSe-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbo(dithioperoxo)thioic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbo(diselenoperoxo)selenoic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbo(ditelluroperoxo)telluroic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(O)=S', iupac => 'benzenecarbothioic O-acid' },
    { smiles => 'C1(=CC=CC=C1)C(S)=O', iupac => 'benzenecarbothioic S-acid' },
    { smiles => 'C1(=CC=CC=C1)C(=S)S', iupac => 'benzenecarbodithioic acid' },
    { smiles => 'C1(=CC=CC=C1)C(=[Se])S', iupac => 'benzenecarboselenothioic S-acid' },
    { smiles => '', iupac => 'benzenecarboselenothioic Se-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=[Se])[SeH]', iupac => 'benzenecarbodiselenoic acid' },
    { smiles => '', iupac => 'benzenecarboselenotelluroic Se-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=[Te])[TeH]', iupac => 'benzenecarboditelluroic acid' },
    { smiles => 'C1(=CC=CC=C1)C(O)=N', iupac => 'benzenecarboximidic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)OO', iupac => 'benzenecarboximidoperoxoic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarboximido(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarboximido(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbo(dithioperoxo)imidic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarboximido(selenothioperoxoic) SeS-acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)S', iupac => 'benzenecarboximidothioic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)[SeH]', iupac => 'benzenecarboximidoselenoic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=N)[TeH]', iupac => 'benzenecarboximidotelluroic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(O)=NN', iupac => 'benzenecarbohydrazonic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(=NN)OO', iupac => 'benzenecarbohydrazonoperoxoic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbohydrazono(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbohydrazono(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'benzenecarbo(ditelluroperoxo)hydrazonic acid', AUTHOR => 1 },
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
