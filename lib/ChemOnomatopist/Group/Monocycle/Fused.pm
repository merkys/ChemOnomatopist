package ChemOnomatopist::Group::Monocycle::Fused;

use strict;
use warnings;

# ABSTRACT: Monocyclic group which is fused to some other cycle
# VERSION

use parent ChemOnomatopist::Group::Monocycle::;

use ChemOnomatopist;
use List::Util qw( all );

sub new
{
    my( $class, $graph, $system, @vertices ) = @_;

    return bless { graph => $graph, system => \$system, vertices => \@vertices }, $class;
}

sub is_aromatic()
{
    my( $self ) = @_;
    return 1 if $self->SUPER::is_aromatic;

    # Check if the whole system can be deemed as aromatic
    my @system_vertices;
    for ($self->{system}->cycles) {
        my @vertices = $_->vertices;
        pop @vertices;
        push @system_vertices, @vertices;
    }
    return ChemOnomatopist::Chain::Circular->new( $self->graph, @system_vertices )->is_aromatic;
}

1;
