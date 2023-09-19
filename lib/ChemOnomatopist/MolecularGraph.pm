package ChemOnomatopist::MolecularGraph;

use strict;
use warnings;

# ABSTRACT: Graph extension for molecular graphs
# VERSION

use parent Graph::Undirected::;

use ChemOnomatopist::Util::Graph;
use List::Util qw( all );

sub new
{
    my $class = shift;
    # use Carp; confess "object method new() for ChemOnomatopist::MolecularGraph is not implemented yet\n" if ref $class;

    if( ref $class ) {
        return bless $class->SUPER::new;
    } elsif( @_ == 1 && $_[0]->isa( Graph::Undirected:: ) ) {
        return bless $_[0], $class;
    } else {
        return bless Graph::Undirected->new( @_, refvertexed => 1 ), $class;
    }
}

# Graph::copy() does not copy edge attributes
sub copy()
{
    my( $self ) = @_;

    my $copy = $self->SUPER::copy;
    for my $edge ($self->edges) {
        next unless $self->has_edge_attributes( @$edge );
        $copy->set_edge_attributes( @$edge, $self->get_edge_attributes( @$edge ) );
    }
    return bless $copy;
}

sub subgraph()
{
    my( $self, @vertices ) = @_;
    if( all { ref $_ eq 'ARRAY' } @vertices ) {
        return bless $self->SUPER::subgraph( @vertices );
    } else {
        return ChemOnomatopist::Util::Graph::subgraph( $self, @vertices );
    }
}

1;
