package ChemOnomatopist::ChainHalf;

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Util::Graph qw(
    tree_branch_positions
);
use Graph::Traversal::DFS;
use List::Util qw( sum0 );

# ABSTRACT: Half of a longest chain
# VERSION

sub new
{
    my( $class, $graph, $number_of_centers, $other_center, @vertices ) = @_;
    my $self = { vertices => \@vertices, graph => $graph, number_of_centers => $number_of_centers, other_center => $other_center };
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

sub vertices()
{
    my( $self ) = @_;
    return @{$self->{vertices}};
}

# Properties

sub branch_positions()
{
    my( $self ) = @_;
    return tree_branch_positions( $self->{graph}, $self->vertices );
}

sub length()
{
    my( $self ) = @_;
    return scalar $self->vertices;
}

sub locant_names()
{
    my( $self ) = @_;

    my $graph = $self->{graph}->copy;
    $graph->delete_path( $self->vertices );

    my @locants;
    for my $vertex ($self->vertices) {
        my @current_locants;
        for my $neighbour ($graph->neighbours( $vertex )) {
            $graph->delete_edge( $vertex, $neighbour );
            push @current_locants, ChemOnomatopist::get_sidechain_name( $graph, $neighbour );
        }
        push @locants, \@current_locants;
    }

    return @locants;
}

sub locant_positions_forward()
{
    my( $self ) = @_;
    my @locants = $self->branch_positions;
    return ($self->length + (defined $self->{other_center})) * @locants + sum0 @locants;
}

sub locant_positions_backward()
{
    my( $self ) = @_;
    my @locants = $self->branch_positions;
    return $self->length * @locants - sum0 @locants;
}

sub number_of_branches_in_sidechains()
{
    my( $self ) = @_;

    # Make a copy with all atoms from candidate chains removed.
    my $copy = $self->{graph}->copy;
    $copy->delete_vertices( $self->vertices );

    return sum0 map  { $_ > 2 ? $_ - 2 : 0 }
                map  { $copy->degree( $_ ) }
                map  { Graph::Traversal::DFS->new( $copy, start => $_ )->dfs }
                grep { $copy->has_vertex( $_ ) }
                map  { $self->{graph}->neighbours( $_ ) }
                     $self->vertices;
}

sub number_of_carbons()
{
    my( $self ) = @_;

    # Make a copy with all atoms from candidate chains removed.
    my $copy = $self->{graph}->copy;
    $copy->delete_vertices( $self->vertices );

    my $C = # grep { is_element( $_, 'C' ) } # FIXME: Will not work in tests, have to enable later.
            map  { Graph::Traversal::DFS->new( $copy, start => $_ )->dfs }
            grep { $copy->has_vertex( $_ ) }
            map  { $self->{graph}->neighbours( $_ ) }
                 $self->vertices;

    return $C;
}

sub number_of_branches()
{
    my( $self ) = @_;
    return scalar $self->branch_positions;
}

1;
