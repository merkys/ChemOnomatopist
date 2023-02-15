#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist::Util::Graph qw(
    graph_longest_paths
    graph_longest_paths_from_vertex
);
use Graph::Undirected;
use Test::More;

plan tests => 4;

my $graph = Graph::Undirected->new;
for (1..10) {
    $graph->add_edge( 0, $_ );
}
is( scalar graph_longest_paths_from_vertex( $graph, 0 ), 10 );
is( scalar graph_longest_paths( $graph ), 90 );

$graph->add_edge( 1, 11 );
is( scalar graph_longest_paths_from_vertex( $graph, 0 ), 1 );
is( scalar graph_longest_paths( $graph ), 9 );
