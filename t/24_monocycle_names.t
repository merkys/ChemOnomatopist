#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Chain::Circular;
use Test::More;

my @cases = grep { !/=/ } keys %ChemOnomatopist::Chain::Circular::names;

plan tests => scalar @cases;

for my $SMILES (@cases) {
    my @vertices = map { { symbol => $_ } } split '', $SMILES;
    my $g = Graph::Undirected->new( refvertexed => 1 );
    $g->add_cycle( @vertices );
    my $cycle = ChemOnomatopist::Chain::Circular->new( $g, @vertices );
    is $cycle->backbone_SMILES, $SMILES;
}
