#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'CCOO', iupac => 'ethaneperoxol' },

    # From BBv2 P-63.7
    { smiles => '[SeH]OCCOO', iupac => '2-(selanyloxy)ethane-1-peroxol', AUTHOR => 1 },
    { smiles => 'NCC(C)(OO)C', iupac => '1-amino-2-methylpropane-2-peroxol' },
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
