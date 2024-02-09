package ChemOnomatopist::Group::Hydrazine;

# ABSTRACT: Hydrazine group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

use Chemistry::OpenSMILES qw( is_double_bond );

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

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }
sub needs_substituent_locants()
{
    my( $self ) = @_;
    return $self->number_of_branches > 1 && $self->number_of_branches < $self->max_valence;
}

# Twice-substituted diazenes do not need suffix locant
sub needs_suffix_locant() { $_[0]->max_valence != 2 || $_[0]->number_of_branches != 2 }

sub prefix()
{
    my( $self ) = @_;
    if( is_double_bond( $self->{graph}, $self->vertices ) ) {
        return 'diazenyl';
    } else {
        return 'hydrazinyl';
    }
}

sub suffix()
{
    my( $self ) = @_;
    if( is_double_bond( $self->{graph}, $self->vertices ) ) {
        return 'diazene';
    } else {
        return 'hydrazine';
    }
}

1;
