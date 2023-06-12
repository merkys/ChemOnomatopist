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

    'CN(N=NN(C)C)C=O' => '1,4,4-trimethyltetraaz-2-ene-1-carbaldehyde', # From BBv2 P-66.6.1.1.3
);

plan tests => scalar( keys %SMILES_cases );

for my $case (sort keys %SMILES_cases) {
    is ChemOnomatopist::get_name( $case ), $SMILES_cases{$case};
}
