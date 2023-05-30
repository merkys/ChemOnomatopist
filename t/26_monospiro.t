#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-24.2.1
    { smiles => 'C1CCCC12CCCCC2', iupac => 'spiro[4.5]decane' },
    { smiles => 'C1CCCC12CCCC2',  iupac => 'spiro[4.4]nonane' },

    # From BBv2 P-24.2.4.1.1
    { smiles => 'C1CCCC12OCCCC2', iupac => '6-oxaspiro[4.5]decane' },
    { smiles => 'C1CCCC12NCCOC2', iupac => '9-oxa-6-azaspiro[4.5]decane' },
    { smiles => 'C1CCCC12CSCNC2', iupac => '7-thia-9-azaspiro[4.5]decane' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    eval { is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
}
