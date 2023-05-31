#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Chain::FromHalves;
use ChemOnomatopist::ChainHalf;
use Graph::Undirected;
use Test::More;

sub chain
{
    my( $graph, $other_center, @vertices ) = @_;
    return ChemOnomatopist::ChainHalf->new( $graph, $other_center, @vertices );
}

plan tests => 4;

my @atoms = map { { symbol => 'C', number => $_ } } 0..99;

# A graph representing octane
my $graph = Graph::Undirected->new( refvertexed => 1 );
$graph->add_path( map { $atoms[$_] } ( 1..8 ) );

my( $A, $B );

$A = chain( $graph, map { $atoms[$_] } ( 5, reverse 1..4 ) );
$B = chain( $graph, map { $atoms[$_] } ( 4, 5..8 ) );

is( ChemOnomatopist::Chain::FromHalves->new( $A, $B )->branch_positions, 0 );
is( ChemOnomatopist::Chain::FromHalves->new( $B, $A )->branch_positions, 0 );

# The following additions transform the graph to 2,4,7-trimethyloctane
$graph->add_edge( map { $atoms[$_] } ( 2, 12 ) );
$graph->add_edge( map { $atoms[$_] } ( 4, 42 ) );
$graph->add_edge( map { $atoms[$_] } ( 7, 72 ) );

$A = chain( $graph, map { $atoms[$_] } ( 5, reverse 1..4 ) );
$B = chain( $graph, map { $atoms[$_] } ( 4, 5..8 ) );

is join( ',', ChemOnomatopist::Chain::FromHalves->new( $A, $B )->branch_positions ), '1,3,6';
is join( ',', ChemOnomatopist::Chain::FromHalves->new( $B, $A )->branch_positions ), '1,4,6';
