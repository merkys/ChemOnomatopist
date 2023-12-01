#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2
    { smiles => 'COCCOCCOCC[SH2]C', iupac => '2,5,8-trioxa-11λ4-thiadodecane' },
    { smiles => 'O1CCOP12OCCCO2', iupac => '1,4,6,10-tetraoxa-5λ5-phosphaspiro[4.5]decane' },
    { smiles => 'C1CCS12CCCCC2', iupac => '4λ4-thiaspiro[3.5]nonane' },
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
