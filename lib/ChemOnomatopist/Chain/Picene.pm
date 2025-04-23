package ChemOnomatopist::Chain::Picene;

# ABSTRACT: Picene chain
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Name;
use ChemOnomatopist::Util;
use ChemOnomatopist::Util::Graph qw(
    graph_without_edge_attributes
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
    # Isolating a path of vertices with degree of 3
    my $subgraph_d3 = $subgraph->subgraph( grep { $subgraph->degree( $_ ) == 3 }
                                                $subgraph->vertices );
    my $d3_start = first { $subgraph_d3->degree( $_ ) == 1 } $subgraph_d3->vertices;
    my @d3_path = reverse Graph::Traversal::DFS->new( $subgraph_d3, start => $d3_start )->dfs;

    my $first = first { !set( @d3_path )->has( $_ ) }
                      $subgraph->neighbours( $d3_path[1] );
    $subgraph->delete_edge( $first, $d3_path[1] );

    # Taking pairs of vertices along the path, deleting edges between them
    while( @d3_path ) {
        $subgraph->delete_edge( shift @d3_path, shift @d3_path );
    }

    my @vertices = reverse Graph::Traversal::DFS->new( $subgraph, start => $first )->dfs;
    return bless { graph => $graph,
                   vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;
    my @candidates = ( $self );

    my @vertices = reverse $self->vertices;
    for (1..6) {
        push @vertices, shift @vertices;
    }
    push @candidates, bless { graph => $self->graph,
                              vertices => \@vertices,
                              candidate_for => $self };

    return @candidates;
}

sub has_form($$)
{
    my( $class, $graph ) = @_;
    my @vertices = $graph->vertices;

    return '' unless @vertices == 22;
    return are_isomorphic( graph_without_edge_attributes( $graph ),
                           $class->ideal_graph,
                           sub { ChemOnomatopist::Util::element( $_[0] ) } );
}

sub ideal_graph($)
{
    my( $class ) = @_;

    my $graph = Graph::Undirected->new( refvertexed => 1 );
    my @vertices = map { { symbol => 'C' } } 1..8;
    $graph->add_path( @vertices );
    $graph->add_cycle( @vertices[0..1], map { { symbol => 'C' } } 1..4 );
    $graph->add_cycle( @vertices[6..7], map { { symbol => 'C' } } 1..4 );
    for (0..2) {
        $graph->add_cycle( @vertices[0+(2*$_)..3+(2*$_)], map { { symbol => 'C' } } 1..2 );
    }

    return $graph;
}

sub number_of_rings() { 5 }

sub prefix() { ChemOnomatopist::Name->new( 'picene' ) }
sub suffix() { $_[0]->prefix }

1;
