package ChemOnomatopist::ChainHalf;

use strict;
use warnings;

use List::Util qw( sum0 );

# ABSTRACT: Half of a longest chain
# VERSION

sub new
{
    my( $class, $graph, $group, @vertices ) = @_;
    my $self = { vertices => \@vertices, graph => $graph, group => $group };
    return bless $self, $class;
}

sub branch_positions()
{
    my( $self ) = @_;
    return tree_branch_positions( $self->{graph}, @{$self->{vertices}} );
}

sub group()
{
    my( $self ) = @_;
    return $self->{group};
}

sub length()
{
    my( $self ) = @_;
    return scalar @{$self->{vertices}};
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

sub number_of_branches()
{
    my( $self ) = @_;
    return scalar $self->branch_positions;
}

1;
