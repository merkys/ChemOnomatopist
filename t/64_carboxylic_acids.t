#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 Table 4.3. Suffixes were turned to full names by prepending them with 'methane'
    { smiles => 'CC(=O)O', iupac => 'methanecarboxylic acid', AUTHOR => 1 },
    { smiles => 'CC(=O)OO', iupac => 'methanecarboperoxoic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarboperoxothioic OO-acid', AUTHOR => 1 },
    { smiles => 'CC(OO)=[Se]', iupac => 'methanecarboperoxoselenoic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(selenoperoxoic) SeO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(selenoperoxoic) OSe-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbothio(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbothio(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarboseleno(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(selenoperoxo)thioic SeO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(selenoperoxo)thioic OSe-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(dithioperoxo)thioic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(diselenoperoxo)selenoic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(ditelluroperoxo)telluroic acid', AUTHOR => 1 },
    { smiles => 'CC(O)=S', iupac => 'methanecarbothioic O-acid', AUTHOR => 1 },
    { smiles => 'CC(S)=O', iupac => 'methanecarbothioic S-acid', AUTHOR => 1 },
    { smiles => 'CC(=S)S', iupac => 'methanecarbodithioic acid', AUTHOR => 1 },
    { smiles => 'CC(=[Se])S', iupac => 'methanecarboselenothioic S-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarboselenothioic Se-acid', AUTHOR => 1 },
    { smiles => 'CC(=[Se])[SeH]', iupac => 'methanecarbodiselenoic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarboselenotelluroic Se-acid', AUTHOR => 1 },
    { smiles => 'CC(=[Te])[TeH]', iupac => 'methanecarboditelluroic acid', AUTHOR => 1 },
    { smiles => 'CC(O)=N', iupac => 'methanecarboximidic acid', AUTHOR => 1 },
    { smiles => 'CC(=N)OO', iupac => 'methanecarboximidoperoxoic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarboximido(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarboximido(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(dithioperoxo)imidic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarboximido(selenothioperoxoic) SeS-acid', AUTHOR => 1 },
    { smiles => 'CC(=N)S', iupac => 'methanecarboximidothioic acid', AUTHOR => 1 },
    { smiles => 'CC(=N)[SeH]', iupac => 'methanecarboximidoselenoic acid', AUTHOR => 1 },
    { smiles => 'CC(=N)[TeH]', iupac => 'methanecarboximidotelluroic acid', AUTHOR => 1 },
    { smiles => 'CC(O)=NN', iupac => 'methanecarbohydrazonic acid', AUTHOR => 1 },
    { smiles => 'CC(=NN)OO', iupac => 'methanecarbohydrazonoperoxoic acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbohydrazono(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbohydrazono(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => '', iupac => 'methanecarbo(ditelluroperoxo)hydrazonic acid', AUTHOR => 1 },
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
