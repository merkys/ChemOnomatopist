package ChemOnomatopist::MolecularGraph;

# ABSTRACT: Graph extension for molecular graphs
# VERSION

use strict;
use warnings;

use ChemOnomatopist::DigraphComparator;
use ChemOnomatopist::Util::Graph;
use Chemistry::OpenSMILES qw(
    is_chiral
    is_chiral_tetrahedral
    mirror
);
use Graph::MoreUtils qw( graph_replace );
use Graph::Undirected;
use List::Util qw( all any first );
use Scalar::Util qw( blessed );
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

    # Detect chiralities that cannot be handled yet
    if( any { is_chiral $_ && !is_chiral_tetrahedral $_ } $self->atoms ) {
        die 'cannot handle chirality other than tetrahedral' . "\n";
    }
    if( any { is_chiral_tetrahedral $_ && exists $_->{chirality_neighbours} &&
              @{$_->{chirality_neighbours}} != 4 } $self->atoms ) {
        die 'cannot handle chiral tetrahedral centers with other than ' .
            '4 neighbours' . "\n";
    }

    # Reorder chiral centers
    for my $atom (grep { is_chiral $_ } $self->atoms) {
        if( !exists $atom->{chirality_neighbours} ) {
            delete  $atom->{chirality};
            next;
        }

        my @chirality_neighbours = @{$atom->{chirality_neighbours}};
        my %order = map { ( $chirality_neighbours[$_] => $_ ) }
                        0..$#chirality_neighbours;
        my @order_now = sort { ChemOnomatopist::DigraphComparator->new( $self, $atom, $a, $b )->compare } @chirality_neighbours;
        if( any { !ChemOnomatopist::DigraphComparator->new( $self, $atom,  $order_now[$_], $order_now[$_+1] )->compare }
                0..@order_now-2 ) {
            # Unimportant chiral center
            delete $atom->{chirality};
            delete $atom->{chirality_neighbours};
            next;
        }

        $atom->{chirality_neighbours} = \@order_now;
        if( join( '', Chemistry::OpenSMILES::Writer::_permutation_order( map { $order{$_} } @order_now ) ) ne '0123' ) {
            mirror $atom;
        }
    }

    return $self;
}

sub atoms() { $_[0]->vertices }
sub bonds() { $_[0]->edges }

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
    graph_replace( $self, $new, @old );

    # Adjust chiral neighbours
    for my $vertex (grep { !blessed $_ && is_chiral $_ } $self->vertices) {
        next unless exists $vertex->{chirality_neighbours};
        my $common = set( @old ) * set( @{$vertex->{chirality_neighbours}} );
        next unless $common->size;

        if( $common->size == 1 ) {
            my $old = first { 1 } @$common;
            @{$vertex->{chirality_neighbours}} =
                map { $_ == $old ? $new : $_ } @{$vertex->{chirality_neighbours}};
        } else {
            # Multiple changes in chirality, cannot retain the center
            delete $vertex->{chirality};
            delete $vertex->{chirality_neighbours};
        }
    }

    return $self;
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
    return any { $_->{charge} && $_->{charge} < 0 } $self->atoms;
}

sub has_positive_charge()
{
    my( $self ) = @_;
    return any { $_->{charge} && $_->{charge} > 0 } $self->atoms;
}

sub is_anion()      {  $_[0]->has_negative_charge && !$_[0]->has_positive_charge }
sub is_cation()     { !$_[0]->has_negative_charge &&  $_[0]->has_positive_charge }
sub is_zwitterion() {  $_[0]->has_negative_charge &&  $_[0]->has_positive_charge }

1;
