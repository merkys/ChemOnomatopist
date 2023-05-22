#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'COCCOCCOCCOCC', iupac => '2,5,8,11-tetraoxatridecane' }, # BBv2 P-12.1
    { smiles => 'C(F)(F)(F)C(F)(F)CO', iupac => '2,2,3,3,3-pentafluoropropan-1-ol' }, # BBv2 P-14.3.4.5

    # From BBv2 P-15.4.3.1
    { smiles => 'COCSSCCOCC[Se]C', iupac => '2,8-dioxa-4,5-dithia-11-selenadodecane' },
    { smiles => '[Si]OCS[Si]', iupac => '2-oxa-4-thia-1,5-disilapentane', AUTHOR => 1 },

    # From BBv2 P-15.4.3.2.1
    { smiles => 'C[Si]C[Si]C[Si]CSCC', iupac => '8-thia-2,4,6-trisiladecane' },
    { smiles => 'C[Si]C[Si]C[Si]COC',  iupac => '2-oxa-4,6,8-trisilanonane', AUTHOR => 1 }, # Fails nondeterministically on some machines

    # From BBv2 P-15.4.3.2.3
    { smiles => 'C[Si]C[Si]C[Si]C[Si]C(=O)O', iupac => '2,4,6,8-tetrasilanonan-1-oic acid', AUTHOR => 1 }, # FIXME
    { smiles => 'C[Si]C[Si]C[Si]C[Si]CCO', iupac => '2,4,6,8-tetrasiladecan-10-ol' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac};
}
