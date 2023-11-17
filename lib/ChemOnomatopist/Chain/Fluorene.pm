package ChemOnomatopist::Chain::Fluorene;

# ABSTRACT: Fluorene or its derivative
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Util::Graph qw( merge_graphs );
use Graph::Undirected;

use parent ChemOnomatopist::Chain::Circular::;

sub ideal_graph($)
{
    my( $class ) = @_;

    my @graphs;
    for (0..1) {
        my $graph = Graph::Undirected->new( refvertexed => 1 );
        $graph->add_cycle( map { { symbol => 'C', number => $_-1 } } 1..6 );
        push @graphs, $graph;
    }
    my $graph = merge_graphs( @graphs );

    # Pick an edge from each graph
    my( $A ) = $graphs[0]->edges;
    my( $B ) = $graphs[1]->edges;

    # Join a pair of atoms with an edge
    $graph->add_edge( $A->[0], $B->[0] );

    # Add a longer arc between other two atoms
    $graph->add_path( $A->[1], { symbol => 'C' }, $B->[1] );

    return $graph;
}

1;
