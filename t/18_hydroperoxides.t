#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'CCOO', iupac => 'ethaneperoxol' },

    # From BBv2 P-56.2
    { smiles => 'CSO', iupac => 'methane-SO-thioperoxol' },
    { smiles => 'C1(=CC=CC=C1)[Se][SeH]', iupac => 'benzenediselenoperoxol' },

    # From BBv2 P-63.4.2.1
    { smiles => 'CCCOS', iupac => 'propane-1-OS-thioperoxol' },
    { smiles => 'CCSS', iupac => 'ethanedithioperoxol' },
    { smiles => 'CS[Se][H]', iupac => 'methane-SSe-selenothioperoxol' },

    # From BBv2 P-63.4.2.2
    { smiles => 'OOCCO', iupac => '2-hydroperoxyethan-1-ol' },
    { smiles => 'SSCC(=O)O', iupac => 'disulfanylacetic acid', AUTHOR => 1 },
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
