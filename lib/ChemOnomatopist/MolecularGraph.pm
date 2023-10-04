package ChemOnomatopist::MolecularGraph;

use strict;
use warnings;

# ABSTRACT: Graph extension for molecular graphs
# VERSION

use parent Graph::Undirected::;

use ChemOnomatopist::Util::Graph;
use List::Util qw( all );
use Set::Object qw( set );

sub new
{
    my $class = shift;

    if( ref $class ) {
        return bless $class->SUPER::new( refvertexed => 1 );
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
    $copy->set_graph_attributes( $self->get_graph_attributes );

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

sub add_group($)
{
    my( $self, $group ) = @_;
    $self->set_graph_attribute( 'groups', [] ) unless $self->has_graph_attribute( 'groups' );

    push @{$self->get_graph_attribute( 'groups' )}, $group;
}

sub delete_group($)
{
    my( $self, $group ) = @_;
    return unless $self->has_graph_attribute( 'groups' );

    @{$self->get_graph_attribute( 'groups' )} =
        grep { $_ != $group } @{$self->get_graph_attribute( 'groups' )};
}

sub groups(@)
{
    my( $self, @vertices ) = @_;
    return () unless $self->has_graph_attribute( 'groups' );

    if( @vertices ) {
        my $vertices = set( @vertices );
        return grep { $vertices <= set( $_->vertices ) }
                    @{$self->get_graph_attribute( 'groups' )};
    } else {
        return @{$self->get_graph_attribute( 'groups' )};
    }
}

1;
