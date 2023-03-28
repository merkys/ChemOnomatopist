#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Chain;
use ChemOnomatopist::ChainHalf;
use Graph::Undirected;
use Test::More;

sub chain($@)
{
    my( $graph, @vertices ) = @_;
    return ChemOnomatopist::ChainHalf->new( $graph, $vertices[0], @vertices );
}

plan tests => 8;

my @atoms = map { { symbol => 'C', number => $_ } } 0..99;

# A graph representing octane
my $graph = Graph::Undirected->new( refvertexed => 1 );
$graph->add_path( map { $atoms[$_] } ( 1..8 ) );

my $A = chain( $graph, map { $atoms[$_] } reverse 1..4 );
my $B = chain( $graph, map { $atoms[$_] } 5..8 );

is $A->locant_positions_forward, 0;
is $A->locant_positions_backward, 0;

is( ChemOnomatopist::Chain->new( $A, $B )->locant_positions, 0 );
is( ChemOnomatopist::Chain->new( $B, $A )->locant_positions, 0 );

# The following additions transform the graph to 2,4,7-trimethyloctane
$graph->add_edge( map { $atoms[$_] } ( 2, 12 ) );
$graph->add_edge( map { $atoms[$_] } ( 4, 42 ) );
$graph->add_edge( map { $atoms[$_] } ( 7, 72 ) );

is $A->locant_positions_forward, 10;
is $A->locant_positions_backward, 8;

is( ChemOnomatopist::Chain->new( $A, $B )->locant_positions, 14 );
is( ChemOnomatopist::Chain->new( $B, $A )->locant_positions, 13 );
