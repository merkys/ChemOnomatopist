#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-66.3.1.1
    { smiles => 'C(CCCC)(=O)NN', iupac => 'pentanehydrazide', AUTHOR => 1 },
    { smiles => 'C(CCC(=O)NN)(=O)NN', iupac => 'butanedihydrazide', AUTHOR => 1 },
    { smiles => 'C1(CCCCC1)C(=O)NN', iupac => 'cyclohexanecarbohydrazide', AUTHOR => 1 },
    { smiles => 'N1(CCCCC1)C(=O)NN', iupac => 'piperidine-1-carbohydrazide', AUTHOR => 1 },
    { smiles => 'CS(=O)(=O)NN', iupac => 'methanesulfonohydrazide', AUTHOR => 1 },
    { smiles => 'C(C)(C(=O)NN)(C(=O)NN)C(=O)NN', iupac => 'ethane-1,1,1-tricarbohydrazide', AUTHOR => 1 },
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
