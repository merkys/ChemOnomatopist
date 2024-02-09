package ChemOnomatopist::Group::Diazene;

# ABSTRACT: Diazene group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;

    my @chains = ( $self, ChemOnomatopist::Group::Diazene->new( $self->graph, reverse $self->vertices ) );
    $chains[1]->{candidate_for} = $self;

    return @chains;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }
sub needs_substituent_locants()
{
    my( $self ) = @_;
    return $self->number_of_branches > 1 && $self->number_of_branches < $self->max_valence;
}
sub needs_suffix_locant() { $_[0]->number_of_branches != 2 }

sub prefix() { 'diazenyl' }
sub suffix() { 'diazene' }

1;
