package ChemOnomatopist::Chain::Phenanthrene;

use strict;
use warnings;

# ABSTRACT: Phenanthrene or its derivative
# VERSION

use parent ChemOnomatopist::Chain::Polyaphene::;

use ChemOnomatopist::Util::Graph qw(
    merge_graphs
    subgraph
);
use Graph::Undirected;
use Set::Object qw( set );

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my $subgraph = subgraph( $graph, map { $_->vertices } @cycles );

    # Constructing the connectivity graph for cycles
    my $connectivity_graph = Graph::Undirected->new( refvertexed => 1 );
    for my $vertex ($subgraph->vertices) {
        next unless $subgraph->degree( $vertex ) == 3;
        $connectivity_graph->add_edge( grep { set( $_->vertices )->has( $vertex ) } @cycles );
    }

    # Detecting the common ring which has three edges having degree 3 at their ends
    my $common_ring;
    for my $cycle (@cycles) {
        next unless scalar( grep { $subgraph->degree( $_->[0] ) == 3 &&
                                   $subgraph->degree( $_->[1] ) == 3 }
                                 subgraph( $graph, $cycle->vertices )->edges ) == 3;
        $common_ring = $cycle;
    }

    # Finding the correct order of cycles, flipping if needed
    my $common_ring_pos = @cycles % 2 ? (@cycles - 1) / 2 : @cycles / 2 - 1;
    my( $start ) = grep { $connectivity_graph->degree( $_ ) == 1 } @cycles;
    my @cycles_in_order = Graph::Traversal::DFS->new( $connectivity_graph,
                                                      start => $start )->dfs;
    if( !@cycles % 2 && $cycles_in_order[$common_ring_pos] != $common_ring ) {
        @cycles_in_order = reverse @cycles_in_order;
    }

    # Finding the atom in the common ring which will get the lowest number
    $subgraph = subgraph( $graph, $common_ring->vertices );
    $subgraph->delete_edges( map { @$_ }
                             map { subgraph( $graph, $_->vertices )->edges }
                                 ( $cycles_in_order[$common_ring_pos-1],
                                   $cycles_in_order[$common_ring_pos+1] ) );
    my( $short_edge ) = grep { $subgraph->degree( $_->[0] ) == 1 &&
                               $subgraph->degree( $_->[1] ) == 1 }
                               $subgraph->edges;
    my( $junction ) = (set( $cycles_in_order[$common_ring_pos-1]->vertices ) *
                       set( @$short_edge ))->members;

    # Finding the candidates of the starting atom
    $subgraph = subgraph( $graph, $cycles_in_order[0]->vertices );
    $subgraph->delete_vertices(   $cycles_in_order[1]->vertices );
    my @candidates = grep { $subgraph->degree( $_ ) == 1 } $subgraph->vertices;

    # Finding the first and the last atom in the enumeration order
    $subgraph = subgraph( $graph, map { $_->vertices } @cycles );
    my $shortest_paths = $subgraph->single_source_shortest_paths( $junction );
    my $min_length;
    my $first;
    for my $vertex (@candidates) {
        my $length = 0;
        my $v = $vertex;
        while(   $shortest_paths->has_vertex_attribute( $v, 'p' ) ) {
            $v = $shortest_paths->get_vertex_attribute( $v, 'p' );
            $length++;
        }
        if( !defined $min_length || $min_length > $length ) {
            $min_length = $length;
            $first = $vertex;
        }
    }

    my( $last ) = grep { $subgraph->degree( $_ ) == 3 }
                       $subgraph->neighbours( $first );

    # Deleting chords and connection between first and last atoms
    $subgraph->delete_edges( map { (set( $cycles_in_order[$_  ]->vertices ) *
                                    set( $cycles_in_order[$_+1]->vertices ))->members }
                                 0..$#cycles-1 );
    $subgraph->delete_edge( $first, $last );

    my @vertices = Graph::Traversal::DFS->new( $subgraph, start => $last )->dfs;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

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
    $graph->add_path( $A->[1], { symbol => 'C' }, { symbol => 'C' }, $B->[1] );

    return $graph;
}

sub prefix()
{
    my( $self ) = @_;
    return 'phenanthrene';
}

sub suffix() { return $_[0]->prefix }

1;
