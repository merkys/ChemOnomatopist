package ChemOnomatopist::Chain::Fluorene;

# ABSTRACT: Fluorene or its derivative
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Name;
use ChemOnomatopist::Util::Graph qw(
    graph_without_edge_attributes
    merge_graphs
);
use Graph::Nauty qw( are_isomorphic );
use Graph::Traversal::DFS;
use Graph::Undirected;
use List::Util qw( first );
use Set::Object qw( set );

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my $subgraph = $graph->subgraph( map { $_->vertices } @cycles );

    my( $cyclopentane, @benzenes ) = sort { $a->length <=> $b->length } @cycles;
    my( $apex ) = (set( $cyclopentane->vertices ) -
                   set( $benzenes[0]->vertices  ) -
                   set( $benzenes[1]->vertices  ))->members;
    $subgraph->delete_vertices( $apex, $subgraph->neighbours( $apex ) );
    my $start = first { $subgraph->degree( $_ ) == 1 } $subgraph->vertices;

    $subgraph = $graph->subgraph( map { $_->vertices } @cycles );
    $subgraph->delete_edge( (set( $cyclopentane->vertices ) * set( $benzenes[0]->vertices ))->members );
    $subgraph->delete_edge( (set( $cyclopentane->vertices ) * set( $benzenes[1]->vertices ))->members );
    $subgraph->delete_edge( $start, (set( $subgraph->neighbours( $start ) ) *
                                     set( $cyclopentane->vertices ))->members );
    my @vertices = Graph::Traversal::DFS->new( $subgraph, start => $start )->dfs;

    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub has_form($$)
{
    my( $class, $graph ) = @_;

    return '' unless $graph->vertices == 13;

    return are_isomorphic( graph_without_edge_attributes( $graph ),
                           $class->ideal_graph,
                           sub { ChemOnomatopist::element( $_[0] ) } );
}

sub ideal_graph($)
{
    my( $class ) = @_;

    my @graphs;
    for (0..1) {
        my $graph = Graph::Undirected->new( refvertexed => 1 );
        $graph->add_cycle( map { { symbol => 'C' } } 1..6 );
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

sub number_of_rings() { 3 }

sub prefix() { return ChemOnomatopist::Name->new( 'fluorene' ) }
sub suffix() { return ChemOnomatopist::Name->new( 'fluorene' ) }

1;
