#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-68.4.2.2
    { smiles => 'COOS', iupac => 'methyldioxidanethiol', AUTHOR => 1 },
    { smiles => 'COOOS', iupac => 'methyltrioxidanethiol' },
    { smiles => '', iupac => 'methyloxidane-SO-thioperoxol', AUTHOR => 1 },
    { smiles => 'CSSO', iupac => 'methyldisulfanol', AUTHOR => 1 },
    { smiles => 'CSSS[SeH]', iupac => 'methyltrisulfaneselenol', AUTHOR => 1 },
    { smiles => 'COOSS', iupac => 'methyldioxidanedithioperoxol', AUTHOR => 1 },
    { smiles => 'COS[Se][SeH]', iupac => 'methoxysulfanediselenoperoxol', AUTHOR => 1 },
    { smiles => 'CSSOO', iupac => 'methyldisulfaneperoxol', AUTHOR => 1 },
    { smiles => 'CSOSO', iupac => '[(methylsulfanyl)oxy]sulfanol', AUTHOR => 1 },
    { smiles => '[Se](O[Se]C1=CC=C(C=C1)O)C1=CC=C(C=C1)O', iupac => '4,4\'-diselenoxanediyldiphenol', AUTHOR => 1 },
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
