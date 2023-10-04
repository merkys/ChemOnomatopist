package ChemOnomatopist::Group::Hydrazide;

use strict;
use warnings;

# ABSTRACT: Hydrazide group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub needs_heteroatom_locants() { return '' }
sub needs_heteroatom_names() { return '' }
sub needs_substituent_locants()
{
    my( $self ) = @_;
    return $self->number_of_branches > 1 && $self->number_of_branches < $self->max_valence;
}

sub prefix()
{
    my( $self ) = @_;
    return 'hydrazidyl';
}

sub suffix()
{
    my( $self ) = @_;
    return 'hydrazide';
}

1;
