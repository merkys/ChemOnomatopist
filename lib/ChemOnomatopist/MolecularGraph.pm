package ChemOnomatopist::MolecularGraph;

# ABSTRACT: Graph extension for molecular graphs
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Util::Graph;
use Graph::MoreUtils qw( graph_replace );
use Graph::Undirected;
use List::Util qw( all any );
use Set::Object qw( set );

use parent Graph::Undirected::;

sub new
{
    my $class = shift;

    my $self;
    if( ref $class ) {
        $self = bless $class->SUPER::new( refvertexed => 1 );
    } elsif( @_ == 1 && $_[0]->isa( Graph::Undirected:: ) ) {
        $self = bless $_[0], $class;
    } else {
        $self = bless Graph::Undirected->new( @_, refvertexed => 1 ), $class;
    }

    return $self;
}

# copy() is overridden as Graph::copy() does not copy edge attributes
sub copy()
{
    my( $self ) = @_;

    my $copy = $self->SUPER::copy;

    for my $edge ($self->edges) {
        next unless $self->has_edge_attributes( @$edge );
        $copy->set_edge_attributes( @$edge, $self->get_edge_attributes( @$edge ) );
    }
    for my $attribute ($self->get_graph_attribute_names) {
        next if $attribute =~ /^_/; # Skip internal attributes
        $copy->set_graph_attribute( $attribute, $self->get_graph_attribute( $attribute ) );
    }

    return $copy;
}

sub replace()
{
    my( $self, $new, @old ) = @_;
    return graph_replace( $self, $new, @old );
}

# TODO: Add edge attributes
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

sub has_negative_charge()
{
    my( $self ) = @_;
    return any { $_->{charge} && $_->{charge} < 0 } $self->vertices;
}

sub has_positive_charge()
{
    my( $self ) = @_;
    return any { $_->{charge} && $_->{charge} > 0 } $self->vertices;
}

sub is_anion()      {  $_[0]->has_negative_charge && !$_[0]->has_positive_charge }
sub is_cation()     { !$_[0]->has_negative_charge &&  $_[0]->has_positive_charge }
sub is_zwitterion() {  $_[0]->has_negative_charge &&  $_[0]->has_positive_charge }

1;
