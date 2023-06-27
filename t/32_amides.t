#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-66.1.1.1.1.1
    { smiles => 'C(CCCCC)(=O)N', iupac => 'hexanamide' },
    { smiles => 'C(CCCC(=O)N)(=O)N', iupac => 'pentanediamide' },

    { smiles => 'C(C(CC(=O)N)C(=O)N)C(=O)N', iupac => 'propane-1,2,3-tricarboxamide', AUTHOR => 1 }, # BBv2 P-66.1.1.1.1.2

    # From BBv2 P-66.1.1.1.2.4
    { smiles => 'C(C=C)(=O)N', iupac => 'prop-2-enamide' },
    { smiles => 'OC(C(=O)N)C', iupac => '2-hydroxypropanamide' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $ok;
    eval { $ok = is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
    diag 'test supposed to fail with AUTHOR_TESTING' if $case->{AUTHOR} && $ok;
}
