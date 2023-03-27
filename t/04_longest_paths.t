#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::ChainHalf;
use ChemOnomatopist::Util::Graph qw(
    graph_longest_paths
    graph_longest_paths_from_vertex
);
use Graph::Undirected;
use Test::More;

sub chainE($@)
{
    my( $graph, @vertices ) = @_;
    return ChemOnomatopist::ChainHalf->new( $graph, $vertices[0], @vertices );
}

sub chainO($@)
{
    my( $graph, @vertices ) = @_;
    return ChemOnomatopist::ChainHalf->new( $graph, $vertices[1], @vertices );
}

plan tests => 12;

my $graph;
my @paths;

$graph = Graph::Undirected->new;
for (1..10) {
    $graph->add_edge( 0, $_ );
}
is( scalar graph_longest_paths_from_vertex( $graph, 0 ), 10 );
is( scalar graph_longest_paths( $graph ), 45 );

$graph->add_edge( 1, 11 );
is( scalar graph_longest_paths_from_vertex( $graph, 0 ), 1 );
is( scalar graph_longest_paths( $graph ), 9 );

# Elongated X-shaped graph with an odd-numbered longest path
$graph = Graph::Undirected->new;
$graph->add_path( 'A'..'C' );
$graph->add_edge( 'A', 'A1' );
$graph->add_edge( 'A', 'A2' );
$graph->add_edge( 'C', 'C1' );
$graph->add_edge( 'C', 'C2' );

is scalar graph_longest_paths( $graph ), 4;
@paths = ChemOnomatopist::rule_greatest_number_of_side_chains_new( $graph,
                                                                   chainO( $graph, 'B', 'A', 'A1' ),
                                                                   chainO( $graph, 'B', 'C', 'C1' ) );
is scalar @paths, 2;

@paths = ChemOnomatopist::rule_greatest_number_of_side_chains_new( $graph,
                                                                   chainO( $graph, 'B', 'A', 'A1' ),
                                                                   chainO( $graph, 'B', 'A', 'A2' ),
                                                                   chainO( $graph, 'B', 'C', 'C1' ) );
is scalar @paths, 2;

# Elongated X-shaped graph with an even-numbered longest path
$graph = Graph::Undirected->new;
$graph->add_path( 'A'..'D' );
$graph->add_edge( 'A', 'A1' );
$graph->add_edge( 'A', 'A2' );
$graph->add_edge( 'D', 'D1' );
$graph->add_edge( 'D', 'D2' );

is scalar graph_longest_paths( $graph ), 4;
@paths = ChemOnomatopist::rule_greatest_number_of_side_chains_new( $graph,
                                                                   chainE( $graph, 'B', 'A', 'A1' ),
                                                                   chainE( $graph, 'C', 'D', 'D1' ) ),
is scalar @paths, 2;

@paths = ChemOnomatopist::rule_greatest_number_of_side_chains_new( $graph,
                                                                   chainE( $graph, 'B', 'A', 'A1' ),
                                                                   chainE( $graph, 'B', 'A', 'A2' ),
                                                                   chainE( $graph, 'C', 'D', 'D1' ) );
is scalar @paths, 2;

# Elongated Y-shaped graph with an even-numbered longest path
$graph = Graph::Undirected->new;
$graph->add_path( 'AA1'..'AA5' );
$graph->add_path( 'AB1'..'AB5' );
$graph->add_path( 'BA1'..'BA5' );
$graph->add_edge( 'A', 'AA1' );
$graph->add_edge( 'A', 'AB1' );
$graph->add_edge( 'B', 'BA1' );
$graph->add_edge( 'A', 'B' );

@paths = ChemOnomatopist::rule_greatest_number_of_side_chains_new( $graph,
                                                                   chainE( $graph, 'A', 'AA1'..'AA5' ),
                                                                   chainE( $graph, 'A', 'AB1'..'AB5' ),
                                                                   chainE( $graph, 'B', 'BA1'..'BA5' ) );
is scalar @paths, 2;

$graph->add_path( 'AA2', 'branch' );
@paths = ChemOnomatopist::rule_greatest_number_of_side_chains_new( $graph,
                                                                   chainE( $graph, 'A', 'AA1'..'AA5' ),
                                                                   chainE( $graph, 'A', 'AB1'..'AB5' ),
                                                                   chainE( $graph, 'B', 'BA1'..'BA5' ) );
is scalar @paths, 2;
