#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-68.3.1.2.1
    { smiles => 'CN(N)C', iupac => '1,1-dimethylhydrazine' },
    { smiles => 'C1(=CC=CC=C1)NN', iupac => 'phenylhydrazine' },
    { smiles => 'N(N)CN' => iupac => '1-hydrazinylmethanamine', AUTHOR => 1 },
    { smiles => 'N(N)C(=O)O', iupac => 'hydrazinecarboxylic acid', AUTHOR => 1 },
    { smiles => 'FN(N(F)F)F', iupac => 'tetrafluorohydrazine', AUTHOR => 1 },
    { smiles => 'N(N)CC#N', iupac => 'hydrazinylacetonitrile', AUTHOR => 1 },

    # From BBv2 P-68.3.1.2.2
    { smiles => 'C(CC)=NN', iupac => 'propylidenehydrazine' },
    { smiles => 'CN(N=C(C)C)C', iupac => '1,1-dimethyl-2-(propan-2-ylidene)hydrazine', AUTHOR => 1 },
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
