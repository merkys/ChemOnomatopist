#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %SMILES_cases = (
    # From BBv2 P-66.6.1.1.1
    'CCCCC=O' => 'pentanal',
    'O=CCCCC=O' => 'pentanedial',

    'O=CCC(C=O)CCC=O' => 'butane-1,2,4-tricarbaldehyde', # From BBv2 P-66.6.1.1.2
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
