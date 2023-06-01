package ChemOnomatopist::ChainHalf;

use strict;
use warnings;

use ChemOnomatopist::Chain;

use parent ChemOnomatopist::Chain::;

# ABSTRACT: Half of a longest chain
# VERSION

sub new
{
    my( $class, $graph, $other_center, @vertices ) = @_;
    my $self = { vertices => \@vertices, graph => $graph, other_center => $other_center, cache => {} };
    return bless $self, $class;
}

# Accessors

# Groups are used to check which halves of a chain can be combined together.
# If a graph contains single center, all halves will share the center.
sub group()
{
    my( $self ) = @_;
    return $self->{vertices}[1 - defined $self->{other_center}];
}

sub _disconnected_chain_graph()
{
    my( $self ) = @_;

    return $self->{_disconnected_chain_graph} if $self->{_disconnected_chain_graph};

    my $graph = $self->graph->copy;
    my @vertices = $self->vertices;

    if( $self->{other_center} ) {
        # Cut the edge to the other center
        $graph->delete_edge( $vertices[0], $self->{other_center} );
    } else {
        # Cut the edges to the other candidates
        for ($graph->neighbours( $vertices[0] )) {
            $graph->delete_edge( $vertices[0], $_ );
        }
    }
    $graph->delete_path( @vertices );

    $self->{_disconnected_chain_graph} = $graph;
    return $graph;
}

1;
