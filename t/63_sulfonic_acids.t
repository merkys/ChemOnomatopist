#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 Table 4.3. Suffixes were turned to full names by prepending them with 'methane'
    { smiles => 'CS(=O)(=O)O', iupac => 'methanesulfonic acid' },
    { smiles => 'CS(=O)(=O)OO', iupac => 'methanesulfonoperoxoic acid' },
    { smiles => 'CS(=O)(=S)OO', iupac => 'methanesulfonoperoxothioic OO-acid' },
    { smiles => 'CS(=O)(=[Se])OO', iupac => 'methanesulfonoperoxoselenoic OO-acid' },
    { smiles => 'CS(=O)(=O)SO', iupac => 'methanesulfono(thioperoxoic) SO-acid' },
    { smiles => 'CS(=O)(=O)OS', iupac => 'methanesulfono(thioperoxoic) OS-acid' },
    { smiles => 'CS(=S)(=S)OO', iupac => 'methanesulfonoperoxodithioic OO-acid' },
    { smiles => 'CS(=O)(=S)SO', iupac => 'methanesulfonothio(thioperoxoic) SO-acid' },
    { smiles => 'CS(=S)(=[Se])OO', iupac => 'methanesulfonoperoxoselenothioic OO-acid' },
    { smiles => 'CS(=[Se])(=[Se])SS', iupac => 'methanesulfono(dithioperoxo)diselenoic acid' },
    { smiles => 'CS(=S)(=S)[Se][SeH]', iupac => 'methanesulfono(diselenoperoxo)dithioic acid' },
    { smiles => 'CS(=[Te])(=[Te])[Te][TeH]', iupac => 'methanesulfono(ditelluroperoxo)ditelluroic acid' },
    { smiles => 'CS(=O)(O)=S', iupac => 'methanesulfonothioic O-acid' },
    { smiles => 'CS(=O)(S)=O', iupac => 'methanesulfonothioic S-acid' },
    { smiles => 'CS(=O)([SeH])=O', iupac => 'methanesulfonoselenoic Se-acid' },
    { smiles => 'CS(O)(=S)=S', iupac => 'methanesulfonodithioic O-acid' },
    { smiles => 'CS(S)(=S)=O', iupac => 'methanesulfonodithioic S-acid' },
    { smiles => 'CS(O)(=[Se])=[Te]', iupac => 'methanesulfonoselenotelluroic O-acid' },
    { smiles => 'CS([SeH])(=O)=[Te]', iupac => 'methanesulfonoselenotelluroic Se-acid' },
    { smiles => 'CS([TeH])(=[Se])=O', iupac => 'methanesulfonoselenotelluroic Te-acid' },
    { smiles => 'CS(=S)(=S)S', iupac => 'methanesulfonotrithioic acid' },
    { smiles => 'CS(=O)(O)=N', iupac => 'methanesulfonimidic acid' },
    { smiles => 'CS(=O)(=N)OO', iupac => 'methanesulfonimidoperoxoic acid' },
    { smiles => 'CS(=S)(=N)OO', iupac => 'methanesulfonimidoperoxothioic OO-acid', AUTHOR => 1 },
    { smiles => 'CS(=O)(=N)SO', iupac => 'methanesulfonimido(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => 'CS(=O)(=N)OS', iupac => 'methanesulfonimido(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => 'CS(O)(=N)=S', iupac => 'methanesulfonimidothioic O-acid' },
    { smiles => 'CS(S)(=N)=O', iupac => 'methanesulfonimidothioic S-acid' },
    { smiles => 'CS(=N)(=S)S', iupac => 'methanesulfonimidodithioic acid' },
    { smiles => 'CS(=N)(=[Se])S', iupac => 'methanesulfonimidoselenothioic S-acid' },
    { smiles => 'CS(=N)(=S)[SeH]', iupac => 'methanesulfonimidoselenothioic Se-acid', AUTHOR => 1 },
    { smiles => 'CS(=N)(=[Te])[TeH]', iupac => 'methanesulfonimidoditelluroic acid' },
    { smiles => 'CS(O)(=N)=N', iupac => 'methanesulfonodiimidic acid' },
    { smiles => 'CS(=N)(=N)OO', iupac => 'methanesulfonodiimidoperoxoic acid' },
    { smiles => 'CS(=N)(=N)SO', iupac => 'methanesulfonodiimido(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => 'CS(=N)(=N)OS', iupac => 'methanesulfonodiimido(thioperoxoic) OS-acid', AUTHOR => 1 },
    { smiles => 'CS(=N)(=N)[SeH]', iupac => 'methanesulfonodiimidoselenoic acid' },
    { smiles => 'CS(=N)(=N)[TeH]', iupac => 'methanesulfonodiimidotelluroic acid' },
    { smiles => 'CS(=O)(O)=NN', iupac => 'methanesulfonohydrazonic acid' },
    { smiles => 'CS(=O)(=NN)OO', iupac => 'methanesulfonohydrazonoperoxoic acid' },
    { smiles => 'CS(=NN)(OO)=S', iupac => 'methanesulfonohydrazonoperoxothioic acid' },
    { smiles => 'CS(O)(=NN)=S', iupac => 'methanesulfonohydrazonothioic O-acid' },
    { smiles => 'CS(S)(=NN)=O', iupac => 'methanesulfonohydrazonothioic S-acid' },
    { smiles => 'CS(=NN)(=NN)O', iupac => 'methanesulfonodihyrazonic acid', AUTHOR => 1 },
    { smiles => 'CS(=NN)(=NN)OO', iupac => 'methanesulfonodihydrazonoperoxoic acid' },
    { smiles => 'CS(=NN)(=NN)SO', iupac => 'methanesulfonodihydrazono(thioperoxoic) SO-acid', AUTHOR => 1 },
    { smiles => 'CS(=NN)(=NN)S', iupac => 'methanesulfonodihydrazonothioic acid' },
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
