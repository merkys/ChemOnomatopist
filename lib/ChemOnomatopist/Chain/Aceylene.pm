package ChemOnomatopist::Chain::Aceylene;

# ABSTRACT: Ace...ylene chain, as per BBv3 P-25.1.2.7
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Chain::Phenanthrene;
use ChemOnomatopist::Chain::Polyacene;
use ChemOnomatopist::Util::Graph qw(
    graph_without_edge_attributes
);
use Graph::Nauty qw( are_isomorphic );
use Graph::Undirected;
use List::Util qw( first );

sub has_form($$)
{
    my( $class, $graph ) = @_;

    return 1 if are_isomorphic( graph_without_edge_attributes( $graph ),
                                $class->ideal_graph_acenaphthylene,
                                sub { ChemOnomatopist::element( $_[0] ) } );
    return 1 if are_isomorphic( graph_without_edge_attributes( $graph ),
                                $class->ideal_graph_aceanthrylene,
                                sub { ChemOnomatopist::element( $_[0] ) } );
    return 1 if are_isomorphic( graph_without_edge_attributes( $graph ),
                                $class->ideal_graph_acephenanthrylene,
                                sub { ChemOnomatopist::element( $_[0] ) } );
    return '';
}

sub ideal_graph_acenaphthylene()
{
    my( $class ) = @_;

    my $graph = Graph::Undirected->new( refvertexed => 1 );
    my @vertices = map { { symbol => 'C' } } 1..11;
    $graph->add_cycle( @vertices );
    my $center = { symbol => 'C' };
    $graph->add_path( $vertices[2], $center, $vertices[6] );
    $graph->add_path( $vertices[2], $center, $vertices[10] );

    return $graph;
}

sub ideal_graph_aceanthrylene()
{
    my( $class ) = @_;

    my $graph = ChemOnomatopist::Chain::Polyacene::ideal_graph( 14 );
    my $d3 = first { $graph->degree( $_ ) == 3 } $graph->vertices;
    my @d2 = grep  { $graph->degree( $_ ) == 2 } $graph->neighbours( $d3 );
    $graph->add_path( $d2[0], { symbol => 'C' }, { symbol => 'C' }, $d2[1] );

    return $graph;
}

sub ideal_graph_acephenanthrylene()
{
    my( $class ) = @_;

    my $graph = ChemOnomatopist::Chain::Phenanthrene::ideal_graph();
    my $d3_subgraph = subgraph( $graph, grep { $graph->degree( $_ ) == 3 } $graph->vertices );
    my $d3 = first { $subgraph->degree( $_ ) == 1 } $d3_subgraph->vertices;
    my @d2 = grep { $graph->degree( $_ ) == 2 } $graph->neighbours( $d3 );
    $graph->add_path( $d2[0], { symbol => 'C' }, { symbol => 'C' }, $d2[1] );

    return $graph;
}

# TODO: Override locants() to give the internal vertex a number

1;
