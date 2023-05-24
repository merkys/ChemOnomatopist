#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'C1CCNC1', iupac => 'pyrrolidine' },

    { smiles => 'C1=CC=CC=C1',  iupac => 'benzene' },
    { smiles => 'C=1C=CC=CC=1', iupac => 'benzene' },
    { smiles => 'c1ccccc1',     iupac => 'benzene' },

    { smiles => 'C1=CC=CC=C1O', iupac => 'phenol' },
    { smiles => 'C(=O)(O)C1=CC=CC=C1', iupac => 'benzoic acid' },

    # From BBv2 P-14.5.1
    { smiles => 'C1CCCCC1(C)CC', iupac => '1-ethyl-1-methylcyclohexane' },
    { smiles => 'CCC1CCC(C)CC1', iupac => '1-ethyl-4-methylcyclohexane' },

    { smiles => 'C1CCCCC1(C(C)(C)C)(CCCC)', iupac => '1-butyl-1-tert-butylcyclohexane' }, # Simplified version of example from P-14.5.1

    { smiles => 'C(=O)(O)CC([Br])([Br])C1CCCCC1', iupac => '3,3-dibromo-3-cyclohexylpropanoic acid' },
    { smiles => 'ClC=1C=CC=CC=1C(F)(F)C(F)(F)F',  iupac => '1-chloro-2-(pentafluoroethyl)benzene', AUTHOR => 1 }, # From BBv2 P-14.3.4.5

    # From BBv2 P-14.5.3
    { smiles => 'C(C)(C)(C)C=1C=CC=C(C(C)CC)C=1', iupac => '1-(butan-2-yl)-3-tert-butylbenzene', AUTHOR => 1 },

    { smiles => 'O=C1NC(=O)NC=C1C', iupac => '5-methylpyrimidine-2,4(1H,3H)-dione', AUTHOR => 1 }, # thymine
    { smiles => 'c1cc(oc1)C=O', iupac => 'furan-2-carbaldehyde', AUTHOR => 1 }, # furfural
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    eval { is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
}
