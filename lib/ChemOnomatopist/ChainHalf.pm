package ChemOnomatopist::ChainHalf;

use strict;
use warnings;

use ChemOnomatopist::Util::Graph qw(
    tree_branch_positions
);
use Graph::Traversal::DFS;
use List::Util qw( sum0 );

# ABSTRACT: Half of a longest chain
# VERSION

sub new
{
    my( $class, $graph, $group, @vertices ) = @_;
    my $self = { vertices => \@vertices, graph => $graph, group => $group };
    return bless $self, $class;
}

# Accessors

sub group()
{
    my( $self ) = @_;
    return $self->{group};
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

sub locant_positions_forward()
{
    my( $self ) = @_;
    my @locants = $self->branch_positions;
    return $self->length * @locants + sum0 @locants;
}

sub locant_positions_backward()
{
    my( $self ) = @_;
    my @locants = $self->branch_positions;
    return ($self->length + 1) * @locants - sum0 @locants;
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
