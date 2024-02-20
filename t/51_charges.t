#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-73.2.2.1.1
    { smiles => '[CH3+]', iupac => 'methylium' },
    { smiles => 'C1(=CC=CC=C1)[Si+](C1=CC=CC=C1)C1=CC=CC=C1', iupac => 'triphenylsilylium', AUTHOR => 1 },
    { smiles => '[CH2+]CC', iupac => 'propylium' },
    { smiles => '[CH+]1CCC1', iupac => 'cyclobutylium' },
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
