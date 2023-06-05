#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'CC(=O)CC(=O)CCCO', iupac => '7-hydroxyheptane-2,4-dione' },
    { smiles => 'CC(C)C(CCC(C)O)CCC(=O)C', iupac => '8-hydroxy-5-propan-2-ylnonan-2-one' },

    { smiles => 'C(O)OO', iupac => '1-hydroperoxymethanol' },
    { smiles => 'OOCC(=O)C', iupac => '1-hydroperoxypropan-2-one' }, # Unchecked
    { smiles => 'OOCCO', iupac => '2-hydroperoxyethan-1-ol' }, # BBv2 P-63.4.2.2

    # From BBv2 P-65.1.2.4
    { smiles => 'CC(=O)CCCC(=O)O', iupac => '5-oxohexanoic acid' },
    { smiles => 'C(=O)(O)C(O)C(C(=O)O)C(=O)C(=O)O', iupac => '1-hydroxy-3-oxopropane-1,2,3-tricarboxylic acid' },

    { smiles => 'O=CCCC(=O)O', iupac => '4-oxobutanoic acid', AUTHOR => 1 }, # BBv2 P-66.6.1.3 - FIXME

    { smiles => 'CCCCCCCCC(C=O)C(CC)O', iupac => '2-(1-hydroxypropyl)decanal' },
    { smiles => 'S=C(CC(=O)O)C', iupac => '3-sulfanylidenebutanoic acid' },

    { smiles => '[N+](=O)([O-])C', iupac => 'nitromethane' }, # BBv2 P-61.5.1
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    eval { is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
}
