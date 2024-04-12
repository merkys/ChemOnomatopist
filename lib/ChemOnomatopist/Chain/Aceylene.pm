package ChemOnomatopist::Chain::Aceylene;

# ABSTRACT: Ace...ylene chain, as per BBv3 P-25.1.2.7
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Chain::Polyacene;
use Graph::Undirected;
use List::Util qw( first );

sub has_form($$)
{
    my( $class, $graph ) = @_;
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

# TODO: Override locants() to give the internal vertex a number

1;
