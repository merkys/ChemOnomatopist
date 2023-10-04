package ChemOnomatopist::Group::Hydrazine;

use strict;
use warnings;

# ABSTRACT: Hydrazine group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;

    my @chains = ( $self, ChemOnomatopist::Group::Hydrazine->new( $self->graph, reverse $self->vertices ) );
    $chains[1]->{candidate_for} = $self;

    return @chains;
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
    return 'hydrazine';
}

sub suffix()
{
    my( $self ) = @_;
    return 'hydrazine';
}

1;
