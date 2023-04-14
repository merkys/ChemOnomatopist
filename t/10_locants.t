#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Chain;
use ChemOnomatopist::ChainHalf;
use Graph::Undirected;
use Test::More;

sub chain
{
    my( $graph, $other_center, @vertices ) = @_;
    return ChemOnomatopist::ChainHalf->new( $graph, $other_center, @vertices );
}

plan tests => 10;

my @atoms = map { { symbol => 'C', number => $_ } } 0..99;

# A graph representing octane
my $graph = Graph::Undirected->new( refvertexed => 1 );
$graph->add_path( map { $atoms[$_] } ( 1..8 ) );

my $A = chain( $graph, map { $atoms[$_] } 5, reverse 1..4 );
my $B = chain( $graph, map { $atoms[$_] } 4, 5..8 );

is $A->locant_positions_forward, 0;
is $A->locant_positions_backward, 0;

is( ChemOnomatopist::Chain->new( $A, $B )->locant_positions, 0 );
is( ChemOnomatopist::Chain->new( $B, $A )->locant_positions, 0 );

# The following additions transform the graph to 2,4,7-trimethyloctane
$graph->add_edge( map { $atoms[$_] } ( 2, 12 ) );
$graph->add_edge( map { $atoms[$_] } ( 4, 42 ) );
$graph->add_edge( map { $atoms[$_] } ( 7, 72 ) );

is $A->locant_positions_forward, 12; # == 5+7
is $A->locant_positions_backward, 6; # == 2+4
is $B->locant_positions_forward,  7; # == 7
is $B->locant_positions_backward, 2; # == 2

is( ChemOnomatopist::Chain->new( $A, $B )->locant_positions, 13 );
is( ChemOnomatopist::Chain->new( $B, $A )->locant_positions, 14 );
