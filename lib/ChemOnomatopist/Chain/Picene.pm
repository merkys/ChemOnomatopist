package ChemOnomatopist::Chain::Picene;

# ABSTRACT: Picene chain
# VERSION

use strict;
use warnings;

use Graph::Undirected;

use parent ChemOnomatopist::Chain::Circular::;

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my $subgraph = $graph->subgraph( map { $_->vertices } @cycles );
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
