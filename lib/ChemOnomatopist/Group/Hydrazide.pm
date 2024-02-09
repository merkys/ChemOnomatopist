package ChemOnomatopist::Group::Hydrazide;

# ABSTRACT: Hydrazide group
# VERSION

use strict;
use warnings;

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

sub locants(@)
{
    my $self = shift;
    return map { $_ == 0 ? "N'" : $_ == 1 ? 'N' : $_ - 1 } @_;
}

sub prefix() { 'hydrazidyl' }
sub suffix() { 'hydrazide' }

1;
