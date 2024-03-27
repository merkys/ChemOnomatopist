#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-45.1.1
    { smiles => 'O(C1=CC(=C(C(=O)O)C=C1)Cl)C1=CC(=C(C(=O)O)C=C1)Cl', iupac => '4,4\'-oxybis(2-chlorobenzoic acid)', AUTHOR => 1 },
    { smiles => 'C(=O)(O)C1=CC(=C(OC2=CC(=C(C(=O)O)C=C2)Cl)C=C1)Cl', iupac => '4-(4-carboxy-2-chlorophenoxy)-2-chlorobenzoic acid', AUTHOR => 1 },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $ok;
    eval { $ok = is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac}, $case->{smiles} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
    if( $case->{AUTHOR} && $ok ) {
        diag 'test supposed to fail with AUTHOR_TESTING' .
             ( $case->{AUTHOR} !~ /^1$/ ? ': ' . $case->{AUTHOR} : '' );
    }
}
