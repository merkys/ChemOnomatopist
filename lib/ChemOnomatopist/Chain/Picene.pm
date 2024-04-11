package ChemOnomatopist::Chain::Picene;

# ABSTRACT: Picene chain
# VERSION

use strict;
use warnings;

use Graph::Traversal::DFS;
use Graph::Undirected;

use parent ChemOnomatopist::Chain::Circular::;

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my $subgraph = $graph->subgraph( map { $_->vertices } @cycles );
    # Isolating a path of vertices with degree of 3
    my $subgraph_d3 = $subgraph->subgraph( grep { $subgraph->degree( $_ ) == 3 }
                                                $subgraph->vertices );
    my $d3_start = first { $subgraph_d3->degree( $_ ) == 1 } $subgraph_d3->vertices;
    my @d3_path = reverse Graph::Traversal::DFS->new( $subgraph_d3, start => $d3_start )->dfs;

    while( @d3_path ) {
        $subgraph->delete_edge( shift @d3_path, shift @d3_path );
    }
}

sub candidates()
{
    my( $self ) = @_;
    my @candidates = ( $self );
    # TODO: Add one more candidate (there are two)

    return @candidates;
}

sub ideal_graph($)
{
    my( $class ) = @_;

    my $graph = Graph::Undirected->new( refvertexed => 1 );
    my @vertices = map { { symbol => 'C' } } 1..8;
    $graph->add_path( @vertices );
    $graph->add_cycle( @vertices[0..1], map { { symbol => 'C' } } 1..4;
    $graph->add_cycle( @vertices[6..7], map { { symbol => 'C' } } 1..4;
    for (0..2) {
        $graph->add_cycle( @vertices[0+(2*$_)..3+(2*$_)], map { { symbol => 'C' } } 1..2 );
    }

    return $graph;
}

sub prefix() { 'picene' }
sub suffix() { $_[0]->prefix }

1;
