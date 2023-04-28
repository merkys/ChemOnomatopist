package ChemOnomatopist::Chain::VertexArray;

use strict;
use warnings;

# ABSTRACT: Chain built upon an array of vertices
# VERSION

use ChemOnomatopist::ChainHalf; # FIXME: Not sure why it is needed

use parent ChemOnomatopist::ChainHalf::;

use List::Util qw( sum0 );

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

# sub branch_positions() # TODO: Maybe need to add 1 to all returned positions?

sub _disconnected_chain_graph()
{
    my( $self ) = @_;

    my $graph = $self->{graph}->copy;
    $graph->delete_path( $self->vertices );

    return $graph;
}

1;
