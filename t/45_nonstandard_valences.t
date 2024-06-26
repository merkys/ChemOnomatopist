#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'C[SH5]', iupac => 'methyl-λ6-sulfane' }, # From BBv3 P-14.1.3

    # From BBv3 P-15.4.3.2.5
    { smiles => 'COCCOCCOCC[SH2]C', iupac => '2,5,8-trioxa-11λ4-thiadodecane' },
    { smiles => 'C[SH2]CCSCCSCCSC', iupac => '2λ4,5,8,11-tetrathiadodecane' },

    { smiles => 'O1CCOP12OCCCO2', iupac => '1,4,6,10-tetraoxa-5λ5-phosphaspiro[4.5]decane' }, # From BBv3 P-44.4.1.5

    # From BBv3 P-24.8.1.1
    { smiles => 'C1CCC12CC[PH3]CC2', iupac => '7λ5-phosphaspiro[3.5]nonane' },
    { smiles => 'C1[SiH2]CC12CC[PH3]CC2', iupac => '7λ5-phospha-2-silaspiro[3.5]nonane' },
    { smiles => 'C1CCS12CCCCC2', iupac => '4λ4-thiaspiro[3.5]nonane' },

    { smiles => 'n1(nccc1c1ccc(F)cc1)c1ccccc1', iupac => '5-(4-fluorophenyl)-1-phenylpyrazole' }, # COD 2201522
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
