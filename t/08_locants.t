#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::ChainHalf;
use Graph::Undirected;
use Test::More;

sub chain($@)
{
    my( $graph, @vertices ) = @_;
    return ChemOnomatopist::ChainHalf->new( $graph, $vertices[1], @vertices );
}

plan tests => 5;

my $graph;
my @paths;
my @chain;

my @atoms = map { { symbol => 'C', number => $_ } } 0..38;

$graph = Graph::Undirected->new( refvertexed => 1 );
$graph->add_path( map { $atoms[$_] } ( 1, 0, 11..17 ) );
$graph->add_path( map { $atoms[$_] } ( 1, 22..27 ) );
$graph->add_path( map { $atoms[$_] } ( 1, 32..37 ) );
$graph->add_edge( map { $atoms[$_] } ( 11, 18 ) );
$graph->add_edge( map { $atoms[$_] } ( 25, 28 ) );
$graph->add_edge( map { $atoms[$_] } ( 32, 38 ) );

@paths = ChemOnomatopist::rule_lowest_numbered_locants_new( $graph, chain( $graph, map { $atoms[$_] } ( 0, 11..17 ) ),
                                                                    chain( $graph, map { $atoms[$_] } ( 0, 1, 22..27 ) ),
                                                                    chain( $graph, map { $atoms[$_] } ( 0, 1, 32..37 ) ) );
is scalar( @paths ), 2;
is join( ',', map { $_->{number} } $paths[0]->vertices ), '0,11,12,13,14,15,16,17';
is join( ',', map { $_->{number} } $paths[1]->vertices ), '0,1,22,23,24,25,26,27';

@chain = ChemOnomatopist::select_main_chain_new( $graph );
is join( ',', map { $_->{number} } @chain ), '27,26,25,24,23,22,1,0,11,12,13,14,15,16,17';

# Figure 7 from Urbonaitė, 2022.
# In the image, however, one branch is held as having priority over another, while in fact they are equal.
$graph = Graph::Undirected->new;
$graph->add_path( 'A'..'G' );
$graph->add_path( 'B', 'H' );
$graph->add_path( 'D', 'I'..'K' );
$graph->add_path( 'J', 'L' );
$graph->add_path( 'E', 'M' );

@paths = ChemOnomatopist::rule_lowest_numbered_locants_new( $graph, chain( $graph, reverse 'A'..'D' ),
                                                                    chain( $graph, 'D'..'G' ),
                                                                    chain( $graph, 'D', 'I'..'K' ) );
is scalar( @paths ), 0;
